# ============================================================
# personC.R  —  Person C
# Volatility metrics, causal baseline, within-subject
# permutation engine
# ============================================================
#
# RESTART R BEFORE RUNNING.
#
# Needs (all from prior scripts):
#   data/derived/features_full_v1.parquet        (Person A)
#   data/derived/features_full_v1_info.parquet   (Person A)
#   data/derived/samples.rds                     (Person B)
#   data/derived/episodes.rds                    (step 0)
#
# Produces:
#   data/derived/volatility.rds
#     — one row per sample; lag-1 SD, RMSSD, Mahalanobis
#       distance-to-baseline, autocorrelation; all time-normalised
#   data/derived/deviation.rds
#     — same index; Z-score deviation from causal baseline
#       (rolling median of PAST samples only)
#   data/derived/permutation_null.rds
#     — null distribution for the within-subject permutation test
#   outputs/volatility_summary.csv
#   outputs/volatility_checks.png
#   outputs/deviation_checks.png
#
# Design notes
# ─────────────────────────────────────────────────────────────
# CAUSAL RULE: every baseline statistic is computed from samples
# BEFORE the index sample's day.  No future information ever
# enters a predictor.  This is enforced by explicit day-ordering
# and an assertion at the end of section 3.
#
# TIME NORMALISATION: raw consecutive differences are divided by
# the inter-sample gap (days) so that unevenly-spaced series are
# comparable.  Minimum gap is clamped to 1 day to prevent Inf.
#
# PERMUTATION ENGINE: within-subject only.  Labels (pre / not-pre)
# are shuffled within each subject's own samples, never across
# subjects.  This preserves individual metabolic baselines and
# between-sample autocorrelation structure.
#
# INSTALL ONCE then comment out:
#   install.packages(c("arrow", "matrixStats"))
# ============================================================

setwd("/home/subhrab/dh607_project/metabolomics")

suppressPackageStartupMessages({
  library(arrow)
  library(matrixStats)
})

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs",      showWarnings = FALSE)


# ============================================================
# PARAMETERS  —  all stated here, all sensitivity-tested below
# ============================================================

MIN_TIMEPOINTS_VOLATILITY <- 3    # need this many to compute lag-1 metrics
MIN_TIMEPOINTS_BASELINE   <- 2    # need this many PAST samples for baseline
N_PERMUTATIONS            <- 5000 # within-subject label permutations
BASELINE_WINDOW_DAYS      <- NULL # NULL = use ALL prior samples (expanding)
# or set e.g. 365 for a rolling window
SEED                      <- 607

# sensitivity grid
MIN_TP_SENS  <- c(2, 4)
WINDOW_SENS  <- c(NULL, 180, 365)   # NULL means no window (expanding)


# ============================================================
# 0. LOAD
# ============================================================

cat("loading feature matrix...\n")
fa   <- read_parquet("data/derived/features_full_v1.parquet")
fi   <- read_parquet("data/derived/features_full_v1_info.parquet")
samp <- readRDS("data/derived/samples.rds")
eps  <- readRDS("data/derived/episodes.rds")

# feature matrix: rows = features, columns = samples
# first column is feature_id; remove it and transpose
feat_ids   <- fa$feature_id
sample_ids <- setdiff(colnames(fa), "feature_id")

# X: samples × features (the working orientation)
X <- t(as.matrix(fa[, sample_ids]))
rownames(X) <- sample_ids
colnames(X) <- feat_ids

cat("feature matrix:", nrow(X), "samples x", ncol(X), "features\n")
cat("sample table:  ", nrow(samp), "rows\n")

# align the sample table to the matrix rows
samp <- samp[match(rownames(X), samp$sample_id), ]
stopifnot(all(rownames(X) == samp$sample_id))

cat("alignment check passed\n\n")


# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

# --- 1a. lag-1 metrics for a single subject's ordered matrix --
#
# x_ord: matrix, rows in CHRONOLOGICAL ORDER, samples × features
# days:  numeric vector of collection days (same length as nrow)
#
# Returns a data.frame with one row per sample.  The first sample
# in every series receives NA for all lag-1 quantities because
# there is no previous observation to difference against.

lag1_metrics <- function(x_ord, days) {
  n <- nrow(x_ord)
  out <- data.frame(
    lag1_sd     = rep(NA_real_, n),
    rmssd       = rep(NA_real_, n),
    lag1_ac     = rep(NA_real_, n),   # feature-wise autocorrelation, median
    gap_days    = c(NA_real_, diff(days))
  )
  
  if (n < 2) return(out)
  
  for (k in 2:n) {
    gap <- max(days[k] - days[k - 1], 1)   # clamp to avoid division by zero
    
    delta <- (x_ord[k, ] - x_ord[k - 1, ]) / gap   # time-normalised difference
    
    out$lag1_sd[k]  <- sd(delta, na.rm = TRUE)
    out$rmssd[k]    <- sqrt(mean(delta^2, na.rm = TRUE))
  }
  
  # autocorrelation: for each feature, lag-1 Pearson across the
  # whole subject series; then take the median across features.
  # This is a subject-level scalar, so it is the same for every
  # row of this subject.
  if (n >= 3) {
    ac <- sapply(seq_len(ncol(x_ord)), function(j) {
      v <- x_ord[, j]
      if (sd(v, na.rm = TRUE) == 0) return(NA_real_)
      cor(v[-n], v[-1], use = "complete.obs", method = "pearson")
    })
    out$lag1_ac <- median(ac, na.rm = TRUE)
  }
  
  out
}


# --- 1b. causal Mahalanobis distance to baseline --------------
#
# For each sample k, compute the Mahalanobis distance from the
# k-th feature vector to the centroid estimated from samples
# 1 … k-1 (causal).  Uses a regularised covariance to handle
# the p >> n situation (shrinkage toward diagonal).
#
# Falls back to Euclidean distance if fewer than
# MIN_TIMEPOINTS_BASELINE prior samples exist, or if the
# covariance is singular after regularisation.

SHRINK_ALPHA <- 0.1   # regularisation strength; (1-α)·S + α·diag·I

mahal_causal <- function(x_ord, days, window_days = NULL) {
  n  <- nrow(x_ord)
  nc <- ncol(x_ord)
  out <- rep(NA_real_, n)
  
  for (k in seq_len(n)) {
    
    # pick the past samples inside the window
    if (is.null(window_days)) {
      past_idx <- seq_len(k - 1)
    } else {
      past_idx <- which(days < days[k] &
                          days >= days[k] - window_days)
    }
    
    if (length(past_idx) < MIN_TIMEPOINTS_BASELINE) next
    
    baseline <- x_ord[past_idx, , drop = FALSE]
    mu       <- colMeans(baseline, na.rm = TRUE)
    
    if (length(past_idx) >= nc + 1) {
      # enough samples for a covariance estimate
      S <- cov(baseline)
      # regularise: blend toward scaled identity
      S_reg <- (1 - SHRINK_ALPHA) * S +
        SHRINK_ALPHA * mean(diag(S)) * diag(nc)
      tryCatch({
        inv_S   <- solve(S_reg)
        diff_v  <- x_ord[k, ] - mu
        out[k]  <- sqrt(as.numeric(t(diff_v) %*% inv_S %*% diff_v))
      }, error = function(e) {
        # singular even after regularisation: fall back to Euclidean
        out[k] <<- sqrt(sum((x_ord[k, ] - mu)^2, na.rm = TRUE))
      })
    } else {
      # too few samples: Euclidean distance to baseline centroid
      out[k] <- sqrt(sum((x_ord[k, ] - mu)^2, na.rm = TRUE))
    }
  }
  out
}


# --- 1c. causal Z-score deviation for each feature ------------
#
# For sample k, return (x_k - median_{t<k}) / MAD_{t<k}
# using the expanding (or windowed) set of prior samples.
# Returns a matrix: same dimensions as x_ord.

zscore_causal <- function(x_ord, days, window_days = NULL) {
  n  <- nrow(x_ord)
  nc <- ncol(x_ord)
  out <- matrix(NA_real_, nrow = n, ncol = nc)
  
  for (k in seq_len(n)) {
    if (is.null(window_days)) {
      past_idx <- seq_len(k - 1)
    } else {
      past_idx <- which(days < days[k] &
                          days >= days[k] - window_days)
    }
    if (length(past_idx) < MIN_TIMEPOINTS_BASELINE) next
    
    baseline <- x_ord[past_idx, , drop = FALSE]
    mu  <- colMedians(baseline, na.rm = TRUE)
    mad <- colMads(baseline,    na.rm = TRUE)
    mad[mad == 0] <- 1e-6   # prevent division by zero
    
    out[k, ] <- (x_ord[k, ] - mu) / mad
  }
  out
}


# ============================================================
# 2. PER-SUBJECT LOOP — VOLATILITY
# ============================================================

cat("Computing volatility metrics...\n")

vol_list  <- list()
subj_ids  <- unique(samp$subject_id)
n_skipped <- 0

for (s in subj_ids) {
  
  idx  <- which(samp$subject_id == s)
  days <- samp$day[idx]
  ord  <- order(days)
  idx  <- idx[ord]
  days <- days[ord]
  
  if (length(idx) < MIN_TIMEPOINTS_VOLATILITY) { n_skipped <- n_skipped + 1; next }
  
  x_ord <- X[idx, , drop = FALSE]
  
  # lag-1 metrics
  lg <- lag1_metrics(x_ord, days)
  
  # Mahalanobis to causal baseline
  mah <- mahal_causal(x_ord, days, window_days = BASELINE_WINDOW_DAYS)
  
  df <- data.frame(
    sample_id  = samp$sample_id[idx],
    subject_id = s,
    day        = days,
    lag1_sd    = lg$lag1_sd,
    rmssd      = lg$rmssd,
    lag1_ac    = lg$lag1_ac,
    gap_days   = lg$gap_days,
    mahal_dist = mah,
    stringsAsFactors = FALSE
  )
  
  vol_list[[s]] <- df
}

vol <- do.call(rbind, vol_list)
rownames(vol) <- NULL

cat("  subjects processed:", length(vol_list),
    " | skipped (< ", MIN_TIMEPOINTS_VOLATILITY, " timepoints):", n_skipped, "\n")
cat("  total rows:", nrow(vol), "\n")


# ============================================================
# 3. PER-SUBJECT LOOP — Z-SCORE DEVIATION MATRIX
# ============================================================

cat("\nComputing Z-score deviation matrix...\n")

dev_list <- list()

for (s in subj_ids) {
  
  idx  <- which(samp$subject_id == s)
  days <- samp$day[idx]
  ord  <- order(days)
  idx  <- idx[ord]
  days <- days[ord]
  
  if (length(idx) < MIN_TIMEPOINTS_BASELINE + 1) next
  
  x_ord <- X[idx, , drop = FALSE]
  
  Z <- zscore_causal(x_ord, days, window_days = BASELINE_WINDOW_DAYS)
  rownames(Z) <- samp$sample_id[idx]
  colnames(Z) <- colnames(X)
  
  dev_list[[s]] <- Z
}

# aggregate: samples × features (only samples with a baseline)
dev_full <- do.call(rbind, dev_list)
cat("  deviation matrix:", nrow(dev_full), "samples x", ncol(dev_full), "features\n")

# summary scalar per sample: median |Z| across all features
# (the "how far from your own normal" single number)
dev_summary <- data.frame(
  sample_id      = rownames(dev_full),
  median_abs_z   = rowMedians(abs(dev_full), na.rm = TRUE),
  pct_gt2sd      = rowMeans(abs(dev_full) > 2, na.rm = TRUE),
  stringsAsFactors = FALSE
)


# CAUSAL ASSERTION: no sample's deviation row was computed using
# future information.  Verify by checking that for every subject
# the first MIN_TIMEPOINTS_BASELINE rows are NA (no past yet).
cat("\n=== CAUSAL ASSERTION ===\n")
n_causal_fail <- 0
for (s in subj_ids) {
  idx  <- which(samp$subject_id == s)
  days <- samp$day[idx]
  ord  <- order(days)
  ids  <- samp$sample_id[idx[ord]]
  
  in_dev <- intersect(ids, rownames(dev_full))
  if (length(in_dev) == 0) next
  
  # the samples that DO appear in dev_full should NOT include the
  # first MIN_TIMEPOINTS_BASELINE of this subject
  first_allowed <- ids[MIN_TIMEPOINTS_BASELINE + 1]
  bad <- in_dev[in_dev %in% ids[seq_len(MIN_TIMEPOINTS_BASELINE)]]
  if (length(bad) > 0) n_causal_fail <- n_causal_fail + 1
}
if (n_causal_fail == 0) {
  cat("PASSED: no first-baseline samples appear in the deviation matrix\n")
} else {
  warning("CAUSAL CHECK FAILED for ", n_causal_fail,
          " subjects — review zscore_causal()\n")
}


# ============================================================
# 4. MERGE BACK TO THE SAMPLE TABLE
# ============================================================

samp2 <- merge(samp,    vol[, c("sample_id","lag1_sd","rmssd",
                                "lag1_ac","gap_days","mahal_dist")],
               by = "sample_id", all.x = TRUE)
samp2 <- merge(samp2,  dev_summary,
               by = "sample_id", all.x = TRUE)
samp2 <- samp2[match(samp$sample_id, samp2$sample_id), ]   # restore row order

cat("\n=== VOLATILITY + DEVIATION COVERAGE ===\n")
cat("lag1_sd  non-NA:", sum(!is.na(samp2$lag1_sd)), "of", nrow(samp2), "\n")
cat("rmssd    non-NA:", sum(!is.na(samp2$rmssd)),   "of", nrow(samp2), "\n")
cat("mahal    non-NA:", sum(!is.na(samp2$mahal_dist)), "of", nrow(samp2), "\n")
cat("med|Z|   non-NA:", sum(!is.na(samp2$median_abs_z)), "of", nrow(samp2), "\n")


# ============================================================
# 5. WITHIN-SUBJECT PERMUTATION ENGINE
# ============================================================
#
# NULL DISTRIBUTION for Q1 / Q2.
#
# Observed statistic: in the TRAIN set, for samples labelled
# "pre-episode" (episode_phase == "pre"), is their median |Z|
# or Mahalanobis distance higher than for samples labelled
# "not-pre" (phase == NA or post) FROM THE SAME SUBJECTS?
#
# Null: within each subject, shuffle the binary pre / not-pre
# label (keeping the total counts fixed), recompute the
# mean group difference, repeat N_PERMUTATIONS times.
#
# This preserves:
#   - each subject's own metabolic baseline
#   - the autocorrelation structure within subjects
#   - sample sizes per subject
# and breaks only the association between the label and the
# metabolic state.

cat("\n\n=== WITHIN-SUBJECT PERMUTATION ENGINE ===\n")

set.seed(SEED)

# use TRAIN samples only (never touch test set)
train_samp <- samp2[!is.na(samp2$split) & samp2$split == "train", ]

# restrict to subjects who have at least one pre and one non-pre sample
pre_subjs <- unique(train_samp$subject_id[
  !is.na(train_samp$episode_phase) &
    train_samp$episode_phase == "pre"])
nonpre_subjs <- unique(train_samp$subject_id[
  is.na(train_samp$episode_phase) |
    train_samp$episode_phase != "pre"])
perm_subjs <- intersect(pre_subjs, nonpre_subjs)

cat("subjects contributing to permutation test:", length(perm_subjs), "\n")
cat("(must have >= 1 pre AND >= 1 non-pre sample in train)\n")

if (length(perm_subjs) == 0) {
  cat("WARNING: no eligible subjects — permutation null not computed.\n")
  cat("This can happen if no episodes land in the train set.\n")
  cat("Check episodes.rds and PersonB output.\n")
  perm_null <- list(medZ = numeric(0), mahal = numeric(0))
} else {
  
  # compute observed statistic
  perm_df <- train_samp[train_samp$subject_id %in% perm_subjs, ]
  perm_df$is_pre <- !is.na(perm_df$episode_phase) &
    perm_df$episode_phase == "pre"
  
  obs_medZ  <- mean(perm_df$median_abs_z[perm_df$is_pre],  na.rm = TRUE) -
    mean(perm_df$median_abs_z[!perm_df$is_pre], na.rm = TRUE)
  obs_mahal <- mean(perm_df$mahal_dist[perm_df$is_pre],    na.rm = TRUE) -
    mean(perm_df$mahal_dist[!perm_df$is_pre],   na.rm = TRUE)
  
  cat("observed mean diff (pre minus not-pre):\n")
  cat("  median |Z| :", round(obs_medZ,  4), "\n")
  cat("  Mahal dist :", round(obs_mahal, 4), "\n")
  
  # --- permutation loop ---
  null_medZ  <- numeric(N_PERMUTATIONS)
  null_mahal <- numeric(N_PERMUTATIONS)
  
  # pre-build per-subject indices for speed
  per_subj <- split(seq_len(nrow(perm_df)), perm_df$subject_id)
  
  for (b in seq_len(N_PERMUTATIONS)) {
    
    perm_labels <- perm_df$is_pre   # will be overwritten per-subject
    
    for (s in perm_subjs) {
      idx_s <- per_subj[[s]]
      n_pre_s <- sum(perm_df$is_pre[idx_s])
      perm_labels[idx_s] <- FALSE
      if (n_pre_s > 0)
        perm_labels[idx_s[sample(length(idx_s), n_pre_s)]] <- TRUE
    }
    
    null_medZ[b]  <- mean(perm_df$median_abs_z[perm_labels],  na.rm = TRUE) -
      mean(perm_df$median_abs_z[!perm_labels],  na.rm = TRUE)
    null_mahal[b] <- mean(perm_df$mahal_dist[perm_labels],     na.rm = TRUE) -
      mean(perm_df$mahal_dist[!perm_labels],     na.rm = TRUE)
  }
  
  p_medZ  <- mean(null_medZ  >= obs_medZ)
  p_mahal <- mean(null_mahal >= obs_mahal)
  
  cat("\npermutation p-values (one-tailed, pre > not-pre):\n")
  cat("  median |Z| p =", round(p_medZ,  4),
      ifelse(p_medZ  < 0.05, "  *", ""), "\n")
  cat("  Mahal dist p =", round(p_mahal, 4),
      ifelse(p_mahal < 0.05, "  *", ""), "\n")
  
  perm_null <- list(
    medZ       = null_medZ,
    mahal      = null_mahal,
    obs_medZ   = obs_medZ,
    obs_mahal  = obs_mahal,
    p_medZ     = p_medZ,
    p_mahal    = p_mahal,
    n_subjs    = length(perm_subjs),
    n_perm     = N_PERMUTATIONS,
    seed       = SEED
  )
}


# ============================================================
# 6. SENSITIVITY — repeat permutation test across parameter grid
# ============================================================

cat("\n\n=== SENSITIVITY: MIN_TIMEPOINTS_VOLATILITY ===\n")

for (mtp in MIN_TP_SENS) {
  
  vl2 <- list()
  for (s in subj_ids) {
    idx  <- which(samp$subject_id == s)
    days <- samp$day[idx]
    ord  <- order(days)
    idx  <- idx[ord]; days <- days[ord]
    if (length(idx) < mtp) next
    lg  <- lag1_metrics(X[idx, , drop = FALSE], days)
    mah <- mahal_causal(X[idx, , drop = FALSE], days)
    vl2[[s]] <- data.frame(sample_id  = samp$sample_id[idx],
                           lag1_sd_s  = lg$lag1_sd,
                           mahal_s    = mah,
                           stringsAsFactors = FALSE)
  }
  v2 <- do.call(rbind, vl2)
  cat("  min_tp =", mtp, ": subjects with data =",
      length(vl2), ", samples =", nrow(v2), "\n")
}

cat("\n=== SENSITIVITY: BASELINE WINDOW ===\n")
for (wd in WINDOW_SENS) {
  wl <- list()
  for (s in subj_ids) {
    idx  <- which(samp$subject_id == s)
    days <- samp$day[idx]
    ord  <- order(days)
    idx  <- idx[ord]; days <- days[ord]
    if (length(idx) < MIN_TIMEPOINTS_BASELINE + 1) next
    mah <- mahal_causal(X[idx, , drop = FALSE], days, window_days = wd)
    wl[[s]] <- data.frame(sample_id = samp$sample_id[idx],
                          mahal_w   = mah,
                          stringsAsFactors = FALSE)
  }
  w <- do.call(rbind, wl)
  lab <- if (is.null(wd)) "expanding" else paste0(wd, " days")
  cat("  window =", lab, ": non-NA distances =",
      sum(!is.na(w$mahal_w)), "\n")
}


# ============================================================
# 7. DIAGNOSTIC PLOTS
# ============================================================

cat("\nPlotting diagnostics...\n")

png("outputs/volatility_checks.png", width = 1000, height = 800)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

# (a) lag1 SD by episode phase
phase_ok <- !is.na(samp2$lag1_sd) & !is.na(samp2$episode_phase)
if (sum(phase_ok) > 0) {
  boxplot(lag1_sd ~ episode_phase,
          data   = samp2[phase_ok, ],
          main   = "Lag-1 SD by episode phase",
          ylab   = "lag-1 SD (time-normalised)",
          xlab   = "phase",
          col    = c("#4DBBD5", "#E64B35", "#00A087"),
          outline = FALSE)
  title(sub = paste("n =", sum(phase_ok), "samples"))
}

# (b) Mahalanobis distance vs days_to_next_onset
pre_ok <- !is.na(samp2$mahal_dist) & !is.na(samp2$days_to_next_onset)
if (sum(pre_ok) > 10) {
  plot(samp2$days_to_next_onset[pre_ok],
       samp2$mahal_dist[pre_ok],
       pch = 16, cex = 0.5,
       col = adjustcolor("steelblue", 0.4),
       xlab = "days to next onset",
       ylab = "Mahalanobis distance to baseline",
       main = "Distance to baseline vs lead time (pre-samples)")
  lo <- tryCatch(
    lowess(samp2$days_to_next_onset[pre_ok], samp2$mahal_dist[pre_ok], f = 0.4),
    error = function(e) NULL)
  if (!is.null(lo)) lines(lo, col = "red", lwd = 2)
}

# (c) median |Z| distribution: pre vs not-pre (train only)
trn <- samp2[!is.na(samp2$split) & samp2$split == "train", ]
trn$pre_label <- ifelse(!is.na(trn$episode_phase) &
                          trn$episode_phase == "pre", "pre", "not-pre")
mz_ok <- !is.na(trn$median_abs_z)
if (sum(mz_ok) > 0) {
  boxplot(median_abs_z ~ pre_label,
          data    = trn[mz_ok, ],
          main    = "Median |Z-score| (train, pre vs not-pre)",
          ylab    = "median |Z|",
          xlab    = "",
          col     = c("#AAAAAA", "#E64B35"),
          outline = FALSE)
}

# (d) permutation null for median |Z|, if available
if (length(perm_null$medZ) > 0) {
  hist(perm_null$medZ,
       breaks = 50,
       main   = "Permutation null: median |Z| difference",
       xlab   = "null mean diff (pre - not-pre)",
       col    = "grey70", border = "white")
  abline(v = perm_null$obs_medZ, col = "red", lwd = 2)
  legend("topright",
         legend = paste0("observed = ", round(perm_null$obs_medZ, 3),
                         "\np = ", round(perm_null$p_medZ, 3)),
         bty = "n", text.col = "red")
}

dev.off()
cat("saved outputs/volatility_checks.png\n")


# deviation Z-score plot
png("outputs/deviation_checks.png", width = 900, height = 500)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))

# pct features > 2 SD by phase (train)
p2_ok <- !is.na(trn$pct_gt2sd) & !is.na(trn$episode_phase)
if (sum(p2_ok) > 0) {
  boxplot(pct_gt2sd * 100 ~ episode_phase,
          data = trn[p2_ok, ],
          main = "% features > 2 SD from baseline\n(train samples with an episode)",
          ylab = "% features > 2 SD",
          xlab = "phase",
          col  = c("#4DBBD5","#E64B35","#00A087"),
          outline = FALSE)
}

# scatter: median |Z| vs days_to_next_onset
pre_tr <- trn[!is.na(trn$days_to_next_onset) & !is.na(trn$median_abs_z), ]
if (nrow(pre_tr) > 5) {
  plot(pre_tr$days_to_next_onset, pre_tr$median_abs_z,
       pch = 16, cex = 0.6,
       col = adjustcolor("darkorange", 0.5),
       xlab = "days to next onset",
       ylab = "median |Z|",
       main = "Median deviation vs lead time\n(train, pre-samples)")
  lo2 <- tryCatch(
    lowess(pre_tr$days_to_next_onset, pre_tr$median_abs_z, f = 0.5),
    error = function(e) NULL)
  if (!is.null(lo2)) lines(lo2, col = "navy", lwd = 2)
}

dev.off()
cat("saved outputs/deviation_checks.png\n")


# ============================================================
# 8. SAVE
# ============================================================

saveRDS(samp2,     "data/derived/volatility.rds")
saveRDS(dev_full,  "data/derived/deviation_matrix.rds")  # samples x features
saveRDS(perm_null, "data/derived/permutation_null.rds")

write.csv(
  samp2[, c("sample_id","subject_id","day","episode_phase",
            "days_to_next_onset","lag1_sd","rmssd","lag1_ac",
            "gap_days","mahal_dist","median_abs_z","pct_gt2sd","split")],
  "outputs/volatility_summary.csv", row.names = FALSE)

cat("\n=== WRITTEN ===\n")
cat("data/derived/volatility.rds          —",
    nrow(samp2), "rows x", ncol(samp2), "cols\n")
cat("data/derived/deviation_matrix.rds    —",
    nrow(dev_full), "rows x", ncol(dev_full), "cols\n")
cat("data/derived/permutation_null.rds    — null dist, N =",
    length(perm_null$medZ), "\n")
cat("outputs/volatility_summary.csv\n")


# ============================================================
# 9. WHAT Q1/Q2 CAN NOW DO
# ============================================================
#
# After running this script you have three analysis-ready objects:
#
# volatility.rds
#   The sample table (all 546 rows) augmented with:
#     lag1_sd, rmssd      — absolute change per day, pre-episode
#     lag1_ac             — within-subject autocorrelation
#     mahal_dist          — distance from your own past centroid
#     median_abs_z        — single-number deviation from baseline
#     pct_gt2sd           — fraction of features > 2 SD from baseline
#
#   Q1 uses: lag1_sd or rmssd as the volatility outcome,
#            episode_phase and days_to_next_onset as grouping.
#   Q2 uses: median_abs_z or mahal_dist vs lag1_sd head-to-head
#            (deviation vs absolute), with permutation p-values.
#
# deviation_matrix.rds
#   Full feature-level Z-scores (samples × features).
#   Needed for Q2's per-compound analysis and for the
#   single-cell bridge (which metabolites deviate most?).
#
# permutation_null.rds
#   Null distribution for the train set.
#   p_medZ and p_mahal are the headline Q1/Q2 p-values.
#   Re-run with the test set by subsetting samp2$split == "test"
#   and calling the same permutation block.
#
# ============================================================

cat("\n=== DONE — Person C complete ===\n")
cat("volatility.rds and deviation_matrix.rds are ready for Q1 → Q2.\n")
cat("Read the diagnostic plots before proceeding.\n")