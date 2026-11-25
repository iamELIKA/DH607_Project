# ============================================================
# Q2_sensitivity.R
# Does DEV performance improve with better baseline quality?
# (n prior timepoints as proxy for baseline reliability)
# ============================================================
#
# RESTART R BEFORE RUNNING.
#
# Needs:
#   data/derived/volatility.rds
#   data/derived/deviation_matrix.rds
#   data/derived/features_full_v1.parquet
#   data/derived/features_full_v1_info.parquet
#   outputs/Q2v2_results.csv            (for reference AUC)
#
# Produces:
#   outputs/Q2_sensitivity_results.csv
#   outputs/Q2_sensitivity_figure.png
#   outputs/Q2_sensitivity_numbers.txt
#
# LOGIC
# ─────────────────────────────────────────────────────────────
# The DEV LOSO-CV AUC was 0.302 — below chance.
# Two explanations:
#   (A) NOISE: too few prior timepoints → unreliable baseline
#              → Z-scores are dominated by estimation error
#   (B) BIOLOGY: deviation genuinely doesn't generalise across
#                subjects regardless of baseline quality
#
# TEST: compute DEV LOSO-CV AUC separately for subjects grouped
# by the number of prior timepoints available at each sample.
#
# "n_prior" = number of samples from that subject BEFORE the
# index sample (causal, same as Person C used for baseline).
#
# Groups:
#   weak:     n_prior = 2–3   (minimum viable baseline)
#   moderate: n_prior = 4–5
#   strong:   n_prior = 6+    (best baseline)
#
# PREDICTION:
#   If (A): AUC(strong) >> AUC(weak), trending upward
#   If (B): AUC flat across groups, all near 0.5 or below
#
# Also run ABS model on same groups as control — ABS should
# be FLAT across groups (absolute levels don't depend on
# prior timepoints), confirming the grouping isn't confounded.
# ============================================================

setwd("/home/subhrab/dh607_project/metabolomics")

suppressPackageStartupMessages({
  library(arrow)
  library(glmnet)
})

dir.create("outputs", showWarnings = FALSE)
set.seed(607)

MAD_FLOOR_PCTILE <- 5
MIN_NONZERO_MAD  <- 0.01
TOP_SCREEN       <- 50
N_BOOT           <- 1000   # fewer for speed; enough for CI


# ============================================================
# 0. LOAD  (identical to Q2_v2.R)
# ============================================================

cat("loading...\n")
samp <- readRDS("data/derived/volatility.rds")
dev  <- readRDS("data/derived/deviation_matrix.rds")

fa         <- read_parquet("data/derived/features_full_v1.parquet")
feat_ids   <- fa$feature_id
sample_ids <- setdiff(colnames(fa), "feature_id")
ABS_full   <- t(as.matrix(fa[, sample_ids]))
rownames(ABS_full) <- sample_ids
colnames(ABS_full) <- feat_ids

fi <- read_parquet("data/derived/features_full_v1_info.parquet")


# ============================================================
# 1. ALIGN + LABEL  (identical to Q2_v2.R)
# ============================================================

common  <- Reduce(intersect, list(rownames(dev),
                                  rownames(ABS_full),
                                  samp$sample_id))
s2      <- samp[match(common, samp$sample_id), ]
DEV_all <- dev[match(common, rownames(dev)),        , drop = FALSE]
ABS_all <- ABS_full[match(common, rownames(ABS_full)), , drop = FALSE]

pre_mask <- !is.na(s2$days_to_next_onset) &
  s2$days_to_next_onset <= 400 &
  !is.na(s2$episode_phase) &
  s2$episode_phase == "pre" &
  s2$diagnosis != "nonIBD"

stable_mask <- is.na(s2$days_to_next_onset) &
  is.na(s2$episode_phase)

s2$outcome <- NA_integer_
s2$outcome[pre_mask]    <- 1L
s2$outcome[stable_mask] <- 0L

lab     <- !is.na(s2$outcome)
s2      <- s2[lab, ]
DEV_lab <- DEV_all[lab, , drop = FALSE]
ABS_lab <- ABS_all[lab, , drop = FALSE]
y       <- s2$outcome

tr <- s2$split == "train" & !is.na(s2$split)
cat("train: pre =", sum(y[tr] == 1),
    ", stable =", sum(y[tr] == 0), "\n")


# ============================================================
# 2. COMPUTE N_PRIOR PER SAMPLE
#
# n_prior = number of samples from the SAME SUBJECT collected
# BEFORE this sample's day (causal, same rule as Person C).
# This is the number of timepoints used to build the baseline
# that generated this sample's Z-score.
# ============================================================

cat("\ncomputing n_prior per sample...\n")

s2$n_prior <- NA_integer_

for (subj in unique(s2$subject_id)) {
  idx  <- which(s2$subject_id == subj)
  days <- s2$day[idx]
  ord  <- order(days)
  idx  <- idx[ord]
  
  # for each sample, count how many EARLIER samples exist
  for (k in seq_along(idx)) {
    s2$n_prior[idx[k]] <- k - 1L
  }
}

cat("n_prior distribution (all labelled samples):\n")
print(table(s2$n_prior))
cat("\nn_prior distribution (train only):\n")
print(table(s2$n_prior[tr]))


# ============================================================
# 3. DEFINE BASELINE QUALITY GROUPS
# ============================================================

# DEV scores require n_prior >= 2 (Person C minimum)
# Group boundaries chosen to give roughly equal group sizes
# while maintaining biological meaning

s2$baseline_group <- NA_character_
s2$baseline_group[s2$n_prior >= 2 & s2$n_prior <= 3] <- "weak (2–3)"
s2$baseline_group[s2$n_prior >= 4 & s2$n_prior <= 5] <- "moderate (4–5)"
s2$baseline_group[s2$n_prior >= 6]                    <- "strong (6+)"

cat("\nbaseline group x outcome (train):\n")
print(table(s2$baseline_group[tr], y[tr],
            dnn = c("group", "outcome")))


# ============================================================
# 4. ANNOTATED FEATURES + MAD CLIPPING
#    (identical pipeline to Q2_v2.R — must be consistent)
# ============================================================

ann_ids <- fi$feature_id[fi$has_hmdb &
                           fi$feature_id %in% colnames(DEV_lab) &
                           fi$feature_id %in% colnames(ABS_lab)]

DEV_ann <- DEV_lab[, ann_ids, drop = FALSE]
ABS_ann <- ABS_lab[, ann_ids, drop = FALSE]

# MAD clipping on TRAIN
dev_tr    <- DEV_ann[tr, , drop = FALSE]
feat_mads <- apply(dev_tr, 2, function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  mad(x, constant = 1)
})
nonzero_mads <- feat_mads[!is.na(feat_mads) & feat_mads > 0]
mad_floor    <- max(quantile(nonzero_mads, MAD_FLOOR_PCTILE / 100),
                    MIN_NONZERO_MAD)
keep_feats   <- ann_ids[!is.na(feat_mads) & feat_mads >= mad_floor]

DEV_clean <- DEV_ann[, keep_feats, drop = FALSE]
ABS_clean <- ABS_ann[, keep_feats, drop = FALSE]

cat("\nfeatures after MAD clip:", length(keep_feats), "\n")


# ============================================================
# 5. HELPERS
# ============================================================

auc_roc <- function(labels, scores) {
  labels <- as.integer(labels)
  ok  <- !is.na(scores) & !is.na(labels)
  if (sum(ok) < 4) return(NA_real_)
  labels <- labels[ok]; scores <- scores[ok]
  if (length(unique(labels)) < 2) return(NA_real_)
  ord <- order(scores, decreasing = TRUE)
  lab <- labels[ord]
  n1  <- sum(lab == 1); n0 <- sum(lab == 0)
  tp  <- cumsum(lab == 1); fp <- cumsum(lab == 0)
  sum(diff(fp / n0) * (tp[-1] + tp[-length(tp)]) / 2 / n1)
}

screen_top <- function(mat, y_tr, top_n) {
  aucs <- sapply(seq_len(ncol(mat)), function(j) {
    v <- mat[, j]
    a <- auc_roc(y_tr, v)
    if (is.na(a)) return(NA_real_)
    max(a, 1 - a)
  })
  names(aucs) <- colnames(mat)
  aucs_sorted <- sort(aucs[!is.na(aucs)], decreasing = TRUE)
  names(head(aucs_sorted, top_n))
}

loso_ridge <- function(X, y, subject_ids) {
  # returns vector of out-of-fold predicted probabilities
  subjs <- unique(subject_ids)
  preds <- rep(NA_real_, length(y))
  
  for (s in subjs) {
    ho  <- which(subject_ids == s)
    trn <- setdiff(seq_along(y), ho)
    if (length(unique(y[trn])) < 2) next
    
    X_tr <- X[trn, , drop = FALSE]
    X_ho <- X[ho,  , drop = FALSE]
    y_tr <- y[trn]
    
    # impute NAs with train column mean
    col_mn <- colMeans(X_tr, na.rm = TRUE)
    for (j in seq_len(ncol(X_tr))) {
      X_tr[is.na(X_tr[, j]), j] <- col_mn[j]
      X_ho[is.na(X_ho[, j]), j] <- col_mn[j]
    }
    
    # scale on train fold
    mu  <- colMeans(X_tr)
    sig <- apply(X_tr, 2, sd); sig[sig == 0] <- 1
    X_tr_sc <- scale(X_tr, center = mu, scale = sig)
    X_ho_sc <- scale(X_ho, center = mu, scale = sig)
    
    cv_fit <- tryCatch(
      cv.glmnet(X_tr_sc, y_tr, family = "binomial",
                alpha = 0, nfolds = 5, type.measure = "auc"),
      error = function(e) NULL)
    if (is.null(cv_fit)) next
    
    preds[ho] <- as.numeric(
      predict(cv_fit, newx = X_ho_sc,
              s = "lambda.min", type = "response"))
  }
  preds
}


# ============================================================
# 6. SENSITIVITY ANALYSIS
#
# For each baseline group, run LOSO-CV for DEV and ABS.
# Use ALL train samples for feature screening (train-global),
# then evaluate only on the samples in the group.
# This separates feature selection from evaluation.
# ============================================================

cat("\n=== SENSITIVITY ANALYSIS ===\n")

# global feature screening on train (same as Q2_v2)
top_feats <- screen_top(DEV_clean[tr, ], y[tr], TOP_SCREEN)
top_feats <- union(top_feats,
                   screen_top(ABS_clean[tr, ], y[tr], TOP_SCREEN))
cat("features used:", length(top_feats), "\n\n")

DEV_top <- DEV_clean[, top_feats, drop = FALSE]
ABS_top <- ABS_clean[, top_feats, drop = FALSE]

groups    <- c("weak (2–3)", "moderate (4–5)", "strong (6+)", "all")
sens_res  <- data.frame()

for (grp in groups) {
  
  if (grp == "all") {
    grp_mask <- tr
  } else {
    grp_mask <- tr & !is.na(s2$baseline_group) &
      s2$baseline_group == grp
  }
  
  n_pre    <- sum(y[grp_mask] == 1, na.rm = TRUE)
  n_stable <- sum(y[grp_mask] == 0, na.rm = TRUE)
  n_subj   <- length(unique(s2$subject_id[grp_mask]))
  
  cat(sprintf("group: %-16s  pre=%2d  stable=%3d  subj=%2d\n",
              grp, n_pre, n_stable, n_subj))
  
  if (n_pre < 3 || n_stable < 3 || n_subj < 4) {
    cat("  SKIP: too few samples\n")
    sens_res <- rbind(sens_res, data.frame(
      group = grp, n_pre = n_pre, n_stable = n_stable,
      n_subj = n_subj,
      auc_dev = NA, auc_abs = NA, diff_dev_abs = NA,
      stringsAsFactors = FALSE))
    next
  }
  
  subj_grp <- s2$subject_id[grp_mask]
  y_grp    <- y[grp_mask]
  
  # DEV
  preds_dev <- loso_ridge(DEV_top[grp_mask, ], y_grp, subj_grp)
  auc_dev   <- auc_roc(y_grp, preds_dev)
  
  # ABS
  preds_abs <- loso_ridge(ABS_top[grp_mask, ], y_grp, subj_grp)
  auc_abs   <- auc_roc(y_grp, preds_abs)
  
  cat(sprintf("  DEV AUC = %s   ABS AUC = %s   diff = %s\n",
              ifelse(is.na(auc_dev), "—", sprintf("%.3f", auc_dev)),
              ifelse(is.na(auc_abs), "—", sprintf("%.3f", auc_abs)),
              ifelse(is.na(auc_dev) | is.na(auc_abs), "—",
                     sprintf("%+.3f", auc_dev - auc_abs))))
  
  sens_res <- rbind(sens_res, data.frame(
    group        = grp,
    n_pre        = n_pre,
    n_stable     = n_stable,
    n_subj       = n_subj,
    auc_dev      = round(auc_dev, 3),
    auc_abs      = round(auc_abs, 3),
    diff_dev_abs = round(auc_dev - auc_abs, 3),
    stringsAsFactors = FALSE))
}

cat("\n=== SENSITIVITY RESULTS TABLE ===\n")
print(sens_res, row.names = FALSE)


# ============================================================
# 7. BOOTSTRAP CIs FOR THE STRONG GROUP
#    (most important comparison — best baseline quality)
# ============================================================

strong_mask <- tr & !is.na(s2$baseline_group) &
  s2$baseline_group == "strong (6+)"
n_strong    <- sum(strong_mask)

if (n_strong >= 10 &&
    sum(y[strong_mask] == 1) >= 3 &&
    sum(y[strong_mask] == 0) >= 3) {
  
  cat("\n=== Bootstrap CI for strong group (N =", N_BOOT, ") ===\n")
  
  subj_strong <- s2$subject_id[strong_mask]
  y_strong    <- y[strong_mask]
  subjs_u     <- unique(subj_strong)
  
  # get full LOSO predictions first
  p_dev_str <- loso_ridge(DEV_top[strong_mask, ], y_strong, subj_strong)
  p_abs_str <- loso_ridge(ABS_top[strong_mask, ], y_strong, subj_strong)
  
  boot_diffs <- replicate(N_BOOT, {
    boot_s <- sample(subjs_u, length(subjs_u), replace = TRUE)
    idx    <- unlist(lapply(boot_s, function(s)
      which(subj_strong == s)))
    ok <- !is.na(p_dev_str[idx]) & !is.na(p_abs_str[idx])
    if (sum(ok) < 5 || length(unique(y_strong[idx][ok])) < 2)
      return(NA_real_)
    auc_roc(y_strong[idx][ok], p_dev_str[idx][ok]) -
      auc_roc(y_strong[idx][ok], p_abs_str[idx][ok])
  })
  
  ci <- quantile(boot_diffs, c(0.025, 0.975), na.rm = TRUE)
  obs_diff <- auc_roc(y_strong, p_dev_str) -
    auc_roc(y_strong, p_abs_str)
  
  cat(sprintf("Strong group DEV−ABS: obs=%+.3f  95%%CI [%+.3f, %+.3f]\n",
              obs_diff, ci[1], ci[2]))
  
  if (ci[1] > 0) {
    cat("INTERPRETATION: DEV significantly BETTER in strong-baseline subjects\n")
    cat("→ The Q2 failure was NOISE — more timepoints fixes it\n")
  } else if (ci[2] < 0) {
    cat("INTERPRETATION: ABS still better even with strong baseline\n")
    cat("→ The Q2 failure is BIOLOGY — deviation doesn't generalise\n")
  } else {
    cat("INTERPRETATION: CI spans 0 in strong group\n")
    cat("→ MIXED: some improvement with better baseline but not definitive\n")
    cat("→ Report the trend (if positive) and note underpowered\n")
  }
  
  strong_ci <- data.frame(
    group = "strong (6+)",
    obs_diff = round(obs_diff, 3),
    ci_lo    = round(ci[1], 3),
    ci_hi    = round(ci[2], 3),
    stringsAsFactors = FALSE)
  
} else {
  cat("\nStrong group too small for bootstrap CI\n")
  strong_ci <- NULL
}


# ============================================================
# 8. FIGURE
# ============================================================

png("outputs/Q2_sensitivity_figure.png", width = 1050, height = 500)
par(mfrow = c(1, 3), mar = c(7, 5, 4, 2))

groups_plot <- groups[groups != "all"]
auc_dev_vec <- sens_res$auc_dev[match(groups_plot, sens_res$group)]
auc_abs_vec <- sens_res$auc_abs[match(groups_plot, sens_res$group)]
n_vec       <- sens_res$n_pre[match(groups_plot, sens_res$group)]

# --- Panel A: DEV AUC by group (the key plot) ---------------
ylim_a <- c(0.2, 1.0)
bp <- barplot(auc_dev_vec,
              names.arg = paste0(groups_plot, "\nn_pre=", n_vec),
              col  = c("#FFCCCC", "#FF6666", "#E64B35"),
              ylim = ylim_a, las = 2,
              ylab = "LOSO-CV AUC",
              main = "DEV model AUC\nby baseline quality",
              cex.names = 0.8)
abline(h = 0.5, lty = 2, col = "grey50")
abline(h = sens_res$auc_dev[sens_res$group == "all"],
       lty = 3, col = "darkred", lwd = 1.5)

# reference line: ABS overall AUC from Q2v2
abline(h = 0.638, lty = 2, col = "#4DBBD5", lwd = 1.5)
legend("topright",
       legend = c("chance (0.5)", "DEV overall", "ABS reference (0.638)"),
       lty    = c(2, 3, 2),
       col    = c("grey50", "darkred", "#4DBBD5"),
       lwd    = c(1, 1.5, 1.5),
       bty    = "n", cex = 0.8)

# --- Panel B: ABS AUC by group (should be flat = control) ---
barplot(auc_abs_vec,
        names.arg = paste0(groups_plot, "\nn_pre=", n_vec),
        col  = c("#CCE5FF", "#66B2FF", "#4DBBD5"),
        ylim = ylim_a, las = 2,
        ylab = "LOSO-CV AUC",
        main = "ABS model AUC\nby baseline quality\n(should be flat = control)",
        cex.names = 0.8)
abline(h = 0.5, lty = 2, col = "grey50")
abline(h = sens_res$auc_abs[sens_res$group == "all"],
       lty = 3, col = "navy", lwd = 1.5)

# --- Panel C: DEV − ABS difference by group -----------------
diff_vec <- auc_dev_vec - auc_abs_vec
bar_col3 <- ifelse(is.na(diff_vec), "grey80",
                   ifelse(diff_vec > 0, "#E64B35", "#4DBBD5"))

bp3 <- barplot(diff_vec,
               names.arg = paste0(groups_plot, "\nn_pre=", n_vec),
               col  = bar_col3,
               ylim = c(min(diff_vec, na.rm = TRUE) - 0.1,
                        max(diff_vec, na.rm = TRUE) + 0.15),
               las  = 2,
               ylab = "AUC(DEV) − AUC(ABS)",
               main = "DEV advantage over ABS\nby baseline quality",
               cex.names = 0.8)
abline(h = 0, lty = 2, col = "black")

# add CI for strong group if computed
if (!is.null(strong_ci)) {
  strong_x <- bp3[length(bp3)]
  arrows(strong_x, strong_ci$ci_lo,
         strong_x, strong_ci$ci_hi,
         angle = 90, code = 3, length = 0.07, lwd = 2,
         col = "black")
  text(strong_x,
       strong_ci$ci_hi + 0.03,
       paste0("CI [", strong_ci$ci_lo, ",", strong_ci$ci_hi, "]"),
       cex = 0.7, col = "black")
}

dev.off()
cat("saved outputs/Q2_sensitivity_figure.png\n")


# ============================================================
# 9. SAVE + METHODS NUMBERS
# ============================================================

write.csv(sens_res, "outputs/Q2_sensitivity_results.csv",
          row.names = FALSE)

sink("outputs/Q2_sensitivity_numbers.txt")
cat("Q2 SENSITIVITY ANALYSIS\n")
cat("generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("QUESTION BEING TESTED\n")
cat("  Does DEV LOSO-CV AUC improve as n_prior (baseline quality) increases?\n")
cat("  If yes → Q2 failure was noise, more timepoints would fix it\n")
cat("  If no  → Q2 failure is biology, deviation doesn't generalise\n\n")

cat("RESULTS\n")
print(sens_res, row.names = FALSE)

cat("\nINTERPRETATION\n")
# auto-interpret the trend
dev_aucs <- sens_res$auc_dev[sens_res$group != "all"]
dev_aucs <- dev_aucs[!is.na(dev_aucs)]
if (length(dev_aucs) >= 2) {
  trend <- dev_aucs[length(dev_aucs)] - dev_aucs[1]
  if (trend > 0.1) {
    cat("  POSITIVE TREND: DEV AUC increases with baseline quality\n")
    cat("  Magnitude:", round(trend, 3), "AUC units\n")
    cat("  Conclusion: noise is likely the dominant explanation\n")
    cat("  More timepoints per subject would likely improve DEV\n")
  } else if (trend < -0.05) {
    cat("  NO POSITIVE TREND: DEV AUC does not improve with more timepoints\n")
    cat("  Conclusion: biology is the dominant explanation\n")
    cat("  Deviation genuinely doesn't generalise across subjects\n")
  } else {
    cat("  FLAT/UNCLEAR TREND: magnitude =", round(trend, 3), "\n")
    cat("  Cannot distinguish noise vs biology at current sample size\n")
    cat("  Report as: insufficient data to resolve the question\n")
  }
}

if (!is.null(strong_ci)) {
  cat("\nSTRONG GROUP BOOTSTRAP CI\n")
  cat("  DEV−ABS observed:", strong_ci$obs_diff, "\n")
  cat("  95% CI: [", strong_ci$ci_lo, ",", strong_ci$ci_hi, "]\n")
}
sink()

cat("wrote outputs/Q2_sensitivity_numbers.txt\n")
cat("\n=== DONE ===\n")
cat("Read outputs/Q2_sensitivity_numbers.txt for the interpretation.\n")
cat("Then share the results table here.\n")

setwd("/home/subhrab/dh607_project/metabolomics")

suppressPackageStartupMessages({
  library(arrow)
  library(glmnet)
})

dir.create("outputs", showWarnings = FALSE)
set.seed(607)

# ============================================================
# PARAMETERS
# ============================================================

MAD_FLOOR_PCTILE <- 5     # exclude features below this MAD percentile
MIN_NONZERO_MAD  <- 0.01  # hard minimum MAD (absolute floor)
N_BOOT           <- 2000  # bootstrap replicates for AUC difference CI
TOP_SCREEN       <- 50    # top features by univariate AUC to pass to ridge
# (screens 405 → 50 to reduce noise, still generous)


# ============================================================
# 0. LOAD
# ============================================================

cat("loading data...\n")
samp <- readRDS("data/derived/volatility.rds")
dev  <- readRDS("data/derived/deviation_matrix.rds")

fa          <- read_parquet("data/derived/features_full_v1.parquet")
feat_ids    <- fa$feature_id
sample_ids  <- setdiff(colnames(fa), "feature_id")
ABS_full    <- t(as.matrix(fa[, sample_ids]))
rownames(ABS_full) <- sample_ids
colnames(ABS_full) <- feat_ids

fi <- read_parquet("data/derived/features_full_v1_info.parquet")

cat("deviation matrix:", nrow(dev), "x", ncol(dev), "\n")
cat("absolute matrix: ", nrow(ABS_full), "x", ncol(ABS_full), "\n\n")


# ============================================================
# 1. ALIGN — same samples in both matrices
# ============================================================

common <- Reduce(intersect, list(rownames(dev),
                                 rownames(ABS_full),
                                 samp$sample_id))
cat("samples with BOTH deviation and absolute:", length(common), "\n")

s2      <- samp[match(common, samp$sample_id), ]
DEV_all <- dev[match(common, rownames(dev)),       , drop = FALSE]
ABS_all <- ABS_full[match(common, rownames(ABS_full)), , drop = FALSE]

stopifnot(all(rownames(DEV_all) == s2$sample_id))
stopifnot(all(rownames(ABS_all) == s2$sample_id))


# ============================================================
# 2. LABEL ASSIGNMENT (strict, same as Q2_v1 fix)
# ============================================================

pre_mask <- !is.na(s2$days_to_next_onset) &
  s2$days_to_next_onset <= 400 &
  !is.na(s2$episode_phase) &
  s2$episode_phase == "pre" &
  s2$diagnosis != "nonIBD"

stable_mask <- is.na(s2$days_to_next_onset) &
  is.na(s2$episode_phase)

s2$outcome <- NA_integer_
s2$outcome[pre_mask]    <- 1L
s2$outcome[stable_mask] <- 0L

cat("outcome counts:\n")
cat("  pre-onset (1):", sum(s2$outcome == 1, na.rm = TRUE), "\n")
cat("  stable    (0):", sum(s2$outcome == 0, na.rm = TRUE), "\n")
cat("  excluded:     ", sum(is.na(s2$outcome)), "\n\n")

# restrict to labelled samples
lab      <- !is.na(s2$outcome)
s2       <- s2[lab, ]
DEV_lab  <- DEV_all[lab, , drop = FALSE]
ABS_lab  <- ABS_all[lab, , drop = FALSE]
y        <- s2$outcome

# TRAIN ONLY for all modelling
tr       <- s2$split == "train" & !is.na(s2$split)
te       <- s2$split == "test"  & !is.na(s2$split)

cat("train: pre =", sum(y[tr] == 1), ", stable =", sum(y[tr] == 0), "\n")
cat("test:  pre =", sum(y[te] == 1), ", stable =", sum(y[te] == 0), "\n\n")


# ============================================================
# 3. ANNOTATED FEATURES ONLY
# ============================================================

ann_ids <- fi$feature_id[fi$has_hmdb &
                           fi$feature_id %in% colnames(DEV_lab) &
                           fi$feature_id %in% colnames(ABS_lab)]
cat("annotated features in both matrices:", length(ann_ids), "\n")

DEV_ann <- DEV_lab[, ann_ids, drop = FALSE]
ABS_ann <- ABS_lab[, ann_ids, drop = FALSE]


# ============================================================
# 4. FIX 2 — MAD CLIPPING
#
# For each annotated feature, compute the median within-subject
# MAD of its DEVIATION Z-scores (across the training subjects).
# Features whose MAD is below the floor are excluded — they
# are features where every subject has near-identical deviation,
# meaning the Z-score is dominated by noise/precision artifacts.
# ============================================================

cat("\n=== FIX 2: MAD clipping ===\n")

# compute per-feature median absolute deviation of Z-scores
# across all labelled TRAINING samples
dev_tr <- DEV_ann[tr, , drop = FALSE]

# for each feature: MAD of its Z-score values across train samples
feat_mads <- apply(dev_tr, 2, function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3) return(NA_real_)
  mad(x, constant = 1)   # constant=1: not scaled to normal
})

# floor = max(5th percentile of non-zero MADs, MIN_NONZERO_MAD)
nonzero_mads <- feat_mads[!is.na(feat_mads) & feat_mads > 0]
mad_floor    <- max(quantile(nonzero_mads, MAD_FLOOR_PCTILE / 100),
                    MIN_NONZERO_MAD)

cat("MAD floor:", round(mad_floor, 4), "\n")
cat("features below MAD floor:",
    sum(!is.na(feat_mads) & feat_mads < mad_floor), "\n")

keep_feats <- ann_ids[!is.na(feat_mads) & feat_mads >= mad_floor]
cat("features retained after MAD clipping:", length(keep_feats), "\n\n")

DEV_clean <- DEV_ann[, keep_feats, drop = FALSE]
ABS_clean <- ABS_ann[, keep_feats, drop = FALSE]

# sanity check: no extreme Z-scores
max_z <- max(abs(DEV_clean[tr, ]), na.rm = TRUE)
cat("max |Z-score| after clipping (train):", round(max_z, 1), "\n")
if (max_z > 50) {
  cat("WARNING: still extreme Z-scores present.",
      "Consider increasing MAD_FLOOR_PCTILE.\n")
} else {
  cat("Z-score range looks reasonable.\n\n")
}


# ============================================================
# 5. UNIVARIATE SCREENING — top TOP_SCREEN features
#
# Screen on TRAIN only. Rank by univariate AUC (both ABS and
# DEV separately). Take union of top TOP_SCREEN from each.
# This reduces 405 → ~100 features before ridge, reducing
# noise while keeping the best candidates from both sides.
# ============================================================

cat("=== Feature screening (train, univariate AUC) ===\n")

# trapezoidal AUC helper
auc_roc <- function(labels, scores) {
  labels <- as.integer(labels)
  if (length(unique(labels[!is.na(labels)])) < 2) return(NA_real_)
  ok  <- !is.na(scores) & !is.na(labels)
  labels <- labels[ok]; scores <- scores[ok]
  ord <- order(scores, decreasing = TRUE)
  lab <- labels[ord]
  n1  <- sum(lab == 1); n0 <- sum(lab == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  tp  <- cumsum(lab == 1); fp <- cumsum(lab == 0)
  sum(diff(fp / n0) * (tp[-1] + tp[-length(tp)]) / 2 / n1)
}

screen <- function(mat, y_tr) {
  aucs <- sapply(seq_len(ncol(mat)), function(j) {
    v <- mat[, j]
    a <- auc_roc(y_tr, v)
    if (is.na(a)) return(NA_real_)
    max(a, 1 - a)   # flip so always >= 0.5
  })
  names(aucs) <- colnames(mat)
  sort(aucs, decreasing = TRUE)
}

auc_dev <- screen(DEV_clean[tr, ], y[tr])
auc_abs <- screen(ABS_clean[tr, ], y[tr])

top_dev <- names(head(auc_dev[!is.na(auc_dev)], TOP_SCREEN))
top_abs <- names(head(auc_abs[!is.na(auc_abs)], TOP_SCREEN))
top_both <- union(top_dev, top_abs)

cat("top", TOP_SCREEN, "DEV features: max AUC =",
    round(auc_dev[1], 3), ", min =",
    round(auc_dev[TOP_SCREEN], 3), "\n")
cat("top", TOP_SCREEN, "ABS features: max AUC =",
    round(auc_abs[1], 3), ", min =",
    round(auc_abs[TOP_SCREEN], 3), "\n")
cat("union of top features passed to ridge:", length(top_both), "\n\n")

DEV_top <- DEV_clean[, top_both, drop = FALSE]
ABS_top <- ABS_clean[, top_both, drop = FALSE]


# ============================================================
# 6. LOSO-CV WITH RIDGE LOGISTIC REGRESSION
#
# For each subject in the TRAIN set:
#   - hold out all their samples
#   - fit ridge on remaining train samples (inner CV for lambda)
#   - predict held-out samples
# Collect all out-of-fold predictions → compute AUC once.
#
# Three models: ABS, DEV, ABS+DEV combined.
# ============================================================

cat("=== LOSO-CV (train set only) ===\n")

loso_cv <- function(X, y, subject_ids, label) {
  
  subjs   <- unique(subject_ids)
  n_subjs <- length(subjs)
  preds   <- rep(NA_real_, length(y))
  
  for (i in seq_along(subjs)) {
    s       <- subjs[i]
    ho_idx  <- which(subject_ids == s)   # held-out
    tr_idx  <- setdiff(seq_along(y), ho_idx)
    
    # need both classes in training fold
    if (length(unique(y[tr_idx])) < 2) next
    
    X_tr <- X[tr_idx, , drop = FALSE]
    y_tr <- y[tr_idx]
    X_ho <- X[ho_idx, , drop = FALSE]
    
    # handle NAs: impute with column mean of training fold
    col_means <- colMeans(X_tr, na.rm = TRUE)
    for (j in seq_len(ncol(X_tr))) {
      na_tr <- is.na(X_tr[, j])
      na_ho <- is.na(X_ho[, j])
      if (any(na_tr)) X_tr[na_tr, j] <- col_means[j]
      if (any(na_ho)) X_ho[na_ho, j] <- col_means[j]
    }
    
    # scale on training fold
    mu  <- colMeans(X_tr)
    sig <- apply(X_tr, 2, sd); sig[sig == 0] <- 1
    X_tr_sc <- scale(X_tr, center = mu, scale = sig)
    X_ho_sc <- scale(X_ho, center = mu, scale = sig)
    
    # ridge (alpha=0), lambda by 5-fold inner CV
    cv_fit <- tryCatch(
      cv.glmnet(X_tr_sc, y_tr, family = "binomial",
                alpha = 0, nfolds = 5,
                type.measure = "auc"),
      error = function(e) NULL)
    
    if (is.null(cv_fit)) next
    
    preds[ho_idx] <- as.numeric(
      predict(cv_fit, newx = X_ho_sc,
              s = "lambda.min", type = "response"))
  }
  
  auc_val <- auc_roc(y, preds)
  cat(sprintf("  %-12s  LOSO-CV AUC = %.3f  (n_subjs=%d, n_pred=%d)\n",
              label, auc_val, n_subjs, sum(!is.na(preds))))
  list(preds = preds, auc = auc_val)
}

# run on train set only
y_tr      <- y[tr]
subj_tr   <- s2$subject_id[tr]

cat("running ABS model...\n")
res_abs  <- loso_cv(ABS_top[tr, ], y_tr, subj_tr, "ABS")

cat("running DEV model...\n")
res_dev  <- loso_cv(DEV_top[tr, ], y_tr, subj_tr, "DEV")

cat("running BOTH model...\n")
res_both <- loso_cv(cbind(ABS_top[tr, ], DEV_top[tr, ]),
                    y_tr, subj_tr, "ABS+DEV")


# ============================================================
# 7. BOOTSTRAP CI FOR AUC DIFFERENCE
#
# Subject-level bootstrap: resample subjects with replacement,
# recompute AUC for each model, compute difference.
# ============================================================

cat("\n=== Bootstrap CI for AUC difference (N =", N_BOOT, ") ===\n")

boot_diff <- function(preds_a, preds_b, y, subject_ids, n_boot) {
  subjs <- unique(subject_ids)
  replicate(n_boot, {
    boot_s  <- sample(subjs, length(subjs), replace = TRUE)
    idx     <- unlist(lapply(boot_s, function(s) which(subject_ids == s)))
    ok      <- !is.na(preds_a[idx]) & !is.na(preds_b[idx])
    if (sum(ok) < 5 || length(unique(y[idx][ok])) < 2)
      return(NA_real_)
    auc_roc(y[idx][ok], preds_a[idx][ok]) -
      auc_roc(y[idx][ok], preds_b[idx][ok])
  })
}

# DEV − ABS: positive = DEV wins
diff_dev_abs <- boot_diff(res_dev$preds, res_abs$preds,
                          y_tr, subj_tr, N_BOOT)
ci_dev_abs   <- quantile(diff_dev_abs, c(0.025, 0.975), na.rm = TRUE)

# BOTH − ABS: positive = combining adds value
diff_both_abs <- boot_diff(res_both$preds, res_abs$preds,
                           y_tr, subj_tr, N_BOOT)
ci_both_abs   <- quantile(diff_both_abs, c(0.025, 0.975), na.rm = TRUE)

cat(sprintf("\nDEV − ABS:  diff = %+.3f  95%% CI [%+.3f, %+.3f]\n",
            res_dev$auc - res_abs$auc,
            ci_dev_abs[1], ci_dev_abs[2]))
cat(sprintf("BOTH − ABS: diff = %+.3f  95%% CI [%+.3f, %+.3f]\n",
            res_both$auc - res_abs$auc,
            ci_both_abs[1], ci_both_abs[2]))

# interpretation
cat("\n--- Q2 answer ---\n")
if (ci_dev_abs[1] > 0) {
  cat("DEV significantly BETTER than ABS (CI excludes 0 from above)\n")
} else if (ci_dev_abs[2] < 0) {
  cat("ABS significantly BETTER than DEV (CI excludes 0 from below)\n")
} else {
  cat("DEV and ABS are NOT significantly different (CI spans 0)\n")
  cat("Both carry similar predictive information at this sample size.\n")
}
if (ci_both_abs[1] > 0) {
  cat("COMBINING DEV+ABS adds significant value over ABS alone.\n")
} else {
  cat("Combining DEV+ABS does NOT significantly improve over ABS alone.\n")
}


# ============================================================
# 8. TEST SET CONFIRMATION (descriptive, one look only)
# ============================================================

cat("\n=== TEST SET (descriptive only) ===\n")

y_te    <- y[te]
subj_te <- s2$subject_id[te]

# fit final model on ALL train, predict test
fit_and_predict_test <- function(X_tr, y_tr, X_te, label) {
  # impute + scale on train
  col_means <- colMeans(X_tr, na.rm = TRUE)
  col_sds   <- apply(X_tr, 2, sd, na.rm = TRUE); col_sds[col_sds == 0] <- 1
  for (j in seq_len(ncol(X_tr))) {
    X_tr[is.na(X_tr[, j]), j] <- col_means[j]
    X_te[is.na(X_te[, j]), j] <- col_means[j]
  }
  X_tr_sc <- scale(X_tr, center = col_means, scale = col_sds)
  X_te_sc <- scale(X_te, center = col_means, scale = col_sds)
  
  cv_fit <- tryCatch(
    cv.glmnet(X_tr_sc, y_tr, family = "binomial",
              alpha = 0, nfolds = 5, type.measure = "auc"),
    error = function(e) NULL)
  if (is.null(cv_fit)) return(NA_real_)
  
  p_te <- as.numeric(predict(cv_fit, newx = X_te_sc,
                             s = "lambda.min", type = "response"))
  auc  <- auc_roc(y_te, p_te)
  cat(sprintf("  %-12s  test AUC = %s\n", label,
              ifelse(is.na(auc), "—", sprintf("%.3f", auc))))
  
  list(preds = p_te, auc = auc, fit = cv_fit)
}

test_abs  <- fit_and_predict_test(ABS_top[tr,], y_tr, ABS_top[te,], "ABS")
test_dev  <- fit_and_predict_test(DEV_top[tr,], y_tr, DEV_top[te,], "DEV")
test_both <- fit_and_predict_test(
  cbind(ABS_top[tr,], DEV_top[tr,]), y_tr,
  cbind(ABS_top[te,], DEV_top[te,]), "ABS+DEV")


# ============================================================
# 9. FEATURE IMPORTANCE FROM FINAL RIDGE MODEL
#
# Use the ABS and DEV full-train ridge fits.
# Features with largest |coefficient| are the important ones.
# These are stable (ridge shrinks but doesn't zero out).
# ============================================================

cat("\n=== Feature importance (ridge coefficients, full train) ===\n")

get_importance <- function(fit_obj, feat_names, fi_table, label) {
  if (is.null(fit_obj) || is.null(fit_obj$fit)) return(NULL)
  coefs <- coef(fit_obj$fit, s = "lambda.min")
  coef_vec <- as.numeric(coefs)[-1]   # drop intercept
  names(coef_vec) <- feat_names
  
  # merge with metabolite names
  meta <- fi_table[match(feat_names, fi_table$feature_id),
                   c("feature_id", "metabolite", "hmdb_clean", "kegg")]
  df <- data.frame(
    feature_id = feat_names,
    metabolite = meta$metabolite,
    hmdb       = meta$hmdb_clean,
    kegg       = meta$kegg,
    coefficient = coef_vec,
    abs_coef    = abs(coef_vec),
    model       = label,
    stringsAsFactors = FALSE
  )
  df <- df[order(-df$abs_coef), ]
  cat("\n", label, "— top 10 features by |coefficient|:\n")
  print(head(df[, c("metabolite","coefficient","abs_coef")], 10),
        row.names = FALSE, digits = 3)
  df
}

imp_abs  <- get_importance(test_abs,  top_both, fi, "ABS")
imp_dev  <- get_importance(test_dev,  top_both, fi, "DEV")


# ============================================================
# 10. RESULTS TABLE
# ============================================================

results <- data.frame(
  model         = c("ABS", "DEV", "ABS+DEV"),
  n_train_subj  = length(unique(subj_tr)),
  n_train_samp  = sum(tr),
  loso_cv_auc   = round(c(res_abs$auc, res_dev$auc, res_both$auc), 3),
  test_auc      = round(c(test_abs$auc, test_dev$auc, test_both$auc), 3),
  stringsAsFactors = FALSE
)

results$diff_vs_abs     <- round(results$loso_cv_auc - results$loso_cv_auc[1], 3)
results$diff_ci_lo      <- c(NA, round(ci_dev_abs[1], 3), round(ci_both_abs[1], 3))
results$diff_ci_hi      <- c(NA, round(ci_dev_abs[2], 3), round(ci_both_abs[2], 3))

cat("\n=== FINAL RESULTS TABLE ===\n")
print(results, row.names = FALSE)


# ============================================================
# 11. FIGURES
# ============================================================

png("outputs/Q2v2_main_figure.png", width = 1100, height = 520)
par(mfrow = c(1, 3), mar = c(6, 5, 4, 2))

# --- Panel A: LOSO-CV AUC comparison -----------------------
auc_vals <- c(res_abs$auc, res_dev$auc, res_both$auc)
bar_cols <- c("#4DBBD5", "#E64B35", "#00A087")
bp <- barplot(auc_vals,
              names.arg = c("ABS\n(log2 levels)",
                            "DEV\n(Z-score deviation)",
                            "ABS+DEV\n(combined)"),
              col  = bar_cols, ylim = c(0.4, 1.0),
              ylab = "LOSO-CV AUC (train set)",
              main = "Q2: Deviation vs Absolute\n(LOSO-CV, ridge logistic, same samples)",
              las  = 1, cex.names = 0.85)
abline(h = 0.5, lty = 2, col = "grey50")

# test AUC as points
test_aucs <- c(test_abs$auc, test_dev$auc, test_both$auc)
ok_te <- !is.na(test_aucs)
if (any(ok_te))
  points(bp[ok_te], test_aucs[ok_te],
         pch = 23, bg = "white", col = bar_cols[ok_te], cex = 1.5)

legend("topright",
       legend = c("LOSO-CV AUC", "test AUC (◇)"),
       pch    = c(NA, 23), col = c(NA, "black"),
       fill   = c("grey70", NA), border = c("black", NA),
       bty = "n", cex = 0.85)


# --- Panel B: bootstrap distribution of AUC difference -----
hist(diff_dev_abs[!is.na(diff_dev_abs)],
     breaks = 50,
     col    = "grey70", border = "white",
     main   = "Bootstrap: AUC(DEV) − AUC(ABS)\n(subject-level, N=2000)",
     xlab   = "AUC difference (DEV − ABS)",
     ylab   = "frequency")
abline(v = res_dev$auc - res_abs$auc, col = "#E64B35", lwd = 2)
abline(v = 0, lty = 2, col = "black")
abline(v = ci_dev_abs, lty = 3, col = "#E64B35")
legend("topright",
       legend = c(sprintf("obs = %+.3f", res_dev$auc - res_abs$auc),
                  sprintf("95%% CI [%+.3f, %+.3f]",
                          ci_dev_abs[1], ci_dev_abs[2])),
       bty = "n", text.col = "#E64B35", cex = 0.85)


# --- Panel C: top features by |ridge coefficient| (DEV) ----
if (!is.null(imp_dev)) {
  top15_dev <- head(imp_dev, 15)
  top15_dev <- top15_dev[order(top15_dev$coefficient), ]
  bar_c2    <- ifelse(top15_dev$coefficient > 0, "#E64B35", "#4DBBD5")
  par(mar = c(5, 13, 4, 2))
  barplot(top15_dev$coefficient,
          horiz     = TRUE,
          names.arg = ifelse(is.na(top15_dev$metabolite) |
                               top15_dev$metabolite == "",
                             top15_dev$feature_id,
                             top15_dev$metabolite),
          las    = 1, col = bar_c2, cex.names = 0.7,
          xlab   = "ridge coefficient (DEV model)",
          main   = "Top 15 features\n(DEV model, full train fit)")
  abline(v = 0, lty = 2)
}

dev.off()
cat("saved outputs/Q2v2_main_figure.png\n")


# ============================================================
# 12. SAVE
# ============================================================

write.csv(results, "outputs/Q2v2_results.csv", row.names = FALSE)

feat_imp_full <- rbind(
  if (!is.null(imp_abs)) cbind(imp_abs, model = "ABS") else NULL,
  if (!is.null(imp_dev)) cbind(imp_dev, model = "DEV") else NULL
)
if (!is.null(feat_imp_full))
  write.csv(feat_imp_full, "outputs/Q2v2_feature_importance.csv",
            row.names = FALSE)

cat("wrote outputs/Q2v2_results.csv\n")
cat("wrote outputs/Q2v2_feature_importance.csv\n")


# ============================================================
# 13. METHODS NUMBERS
# ============================================================

sink("outputs/Q2v2_methods_numbers.txt")
cat("Q2 v2 METHODS NUMBERS\n")
cat("generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("DESIGN\n")
cat("  method: LOSO-CV with ridge logistic regression (alpha=0)\n")
cat("  same samples for DEV and ABS:", sum(tr), "train,", sum(te), "test\n")
cat("  annotated features (pre-MAD clip):", length(ann_ids), "\n")
cat("  MAD floor (5th pctile):", round(mad_floor, 4), "\n")
cat("  features after MAD clip:", length(keep_feats), "\n")
cat("  features passed to ridge (top screen):", length(top_both), "\n\n")

cat("SAMPLE COUNTS\n")
cat("  train: pre =", sum(y[tr] == 1),
    ", stable =", sum(y[tr] == 0),
    ", subjects =", length(unique(subj_tr)), "\n")
cat("  test:  pre =", sum(y[te] == 1),
    ", stable =", sum(y[te] == 0),
    ", subjects =", length(unique(subj_te)), "\n\n")

cat("LOSO-CV AUC\n")
cat("  ABS:     ", round(res_abs$auc,  3), "\n")
cat("  DEV:     ", round(res_dev$auc,  3), "\n")
cat("  ABS+DEV: ", round(res_both$auc, 3), "\n\n")

cat("AUC DIFFERENCE (DEV − ABS)\n")
cat("  observed:", round(res_dev$auc - res_abs$auc, 3), "\n")
cat("  95% CI:  [", round(ci_dev_abs[1], 3), ",",
    round(ci_dev_abs[2], 3), "]\n")
cat("  CI spans 0:", ci_dev_abs[1] < 0 & ci_dev_abs[2] > 0, "\n\n")

cat("TEST AUC\n")
cat("  ABS:     ",
    ifelse(is.na(test_abs$auc),  "—", round(test_abs$auc,  3)), "\n")
cat("  DEV:     ",
    ifelse(is.na(test_dev$auc),  "—", round(test_dev$auc,  3)), "\n")
cat("  ABS+DEV: ",
    ifelse(is.na(test_both$auc), "—", round(test_both$auc, 3)), "\n\n")

if (!is.null(imp_dev) && nrow(imp_dev) > 0) {
  cat("TOP DEV FEATURES (by |ridge coefficient|)\n")
  for (i in seq_len(min(10, nrow(imp_dev)))) {
    cat(sprintf("  %2d. %-30s coef=%+.4f\n",
                i,
                ifelse(is.na(imp_dev$metabolite[i]) |
                         imp_dev$metabolite[i] == "",
                       imp_dev$feature_id[i],
                       imp_dev$metabolite[i]),
                imp_dev$coefficient[i]))
  }
}
sink()

cat("wrote outputs/Q2v2_methods_numbers.txt\n")
cat("\n=== DONE — Q2 v2 complete ===\n")
cat("Key result: outputs/Q2v2_results.csv\n")
cat("Feature list for single-cell bridge: outputs/Q2v2_feature_importance.csv\n")