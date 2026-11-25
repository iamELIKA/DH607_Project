# ============================================================
# Q1_volatility.R
# Does gut chemistry become more erratic before dysbiosis onset?
# ============================================================
#
# RESTART R BEFORE RUNNING.
#
# Needs:
#   data/derived/volatility.rds          (Person C)
#   data/derived/permutation_null.rds    (Person C)
#   data/derived/episodes.rds            (step 0)
#   output/leadtime_bins.csv             (step 0)
#
# Produces:
#   outputs/Q1_results.csv               — headline numbers, one row per bin
#   outputs/Q1_permutation_bins.rds      — full null dist per bin
#   outputs/Q1_main_figure.png           — publication figure
#   outputs/Q1_supplementary.png         — sensitivity / diagnostics
#   outputs/Q1_methods_numbers.txt       — every number cited in the methods
#
# Analysis logic
# ─────────────────────────────────────────────────────────────
# PRIMARY OUTCOME:  median |Z-score| (median_abs_z)
#   = how far the metabolome is from that subject's own past
#     median, expressed in units of that subject's own MAD.
#   This is the Q2 "deviation" metric; Q1 asks whether it
#   rises before onset, Q2 will compare it to absolute levels.
#
# SECONDARY OUTCOMES (reported alongside):
#   lag1_sd   — absolute change per day between visits
#   rmssd     — root-mean-square successive difference per day
#   mahal_dist — Mahalanobis distance to causal baseline
#
# REFERENCE GROUP: "stable" samples from the SAME SUBJECTS
#   who have an episode, taken at lead time > 400 days or from
#   subjects who never had an episode.  Using control samples
#   from subjects with episodes preserves the within-person
#   comparison and avoids confounding by IBD status.
#
# BINNING: the five lead-time bins from step_0
#   [0,45], (45,75], (75,120], (120,200], (200,400]
#   The question is whether any bin differs from stable.
#
# TEST: within-subject permutation, done PER BIN.
#   For each bin, subjects contributing a pre-sample to that bin
#   are identified.  Within each such subject, the bin label is
#   shuffled (keeping the number of "pre" samples fixed).
#   This preserves individual baselines and autocorrelation.
#   One-tailed: H_a is pre > stable.
#
# MULTIPLE COMPARISONS: Benjamini-Hochberg across the 5 bins,
#   reported but NOT used as a gate (5 bins, low power).
#
# TRAIN ONLY for all inferential analysis.
# TEST SET held out; a brief confirmation table is printed but
# not interpreted — save it for the write-up.
# ============================================================

setwd("/home/subhrab/dh607_project/metabolomics")

dir.create("outputs", showWarnings = FALSE)

set.seed(607)
N_PERM <- 5000

# ---- lead-time bin edges (must match step_0) ---------------
BIN_EDGES   <- c(0, 45, 75, 120, 200, 400)
BIN_LABELS  <- c("0–45 d", "45–75 d", "75–120 d", "120–200 d", "200–400 d")

# minimum pre-samples per bin to report a p-value
MIN_PRE_SAMPLES <- 5
MIN_PRE_SUBJECTS <- 3


# ============================================================
# 0. LOAD
# ============================================================

samp  <- readRDS("data/derived/volatility.rds")
eps   <- readRDS("data/derived/episodes.rds")

cat("sample table:", nrow(samp), "rows x", ncol(samp), "cols\n")
cat("episodes:    ", nrow(eps),  "\n")

# confirm required columns exist
need <- c("sample_id","subject_id","day","split","episode_phase",
          "days_to_next_onset","median_abs_z","lag1_sd","rmssd","mahal_dist")
missing_cols <- setdiff(need, colnames(samp))
if (length(missing_cols) > 0)
  stop("Missing columns from volatility.rds: ", paste(missing_cols, collapse = ", "))

cat("all required columns present\n\n")


# ============================================================
# 1. DEFINE ANALYSIS SETS
# ============================================================

# --- 1a. TRAIN set only for inference -----------------------
trn <- samp[!is.na(samp$split) & samp$split == "train", ]
cat("train samples:", nrow(trn), "\n")

# --- 1b. Assign lead-time bins to pre-samples ---------------
trn$leadtime_bin <- NA_character_
pre_mask <- !is.na(trn$days_to_next_onset) &
  trn$days_to_next_onset >= BIN_EDGES[1] &
  trn$days_to_next_onset <= BIN_EDGES[length(BIN_EDGES)]

trn$leadtime_bin[pre_mask] <- as.character(
  cut(trn$days_to_next_onset[pre_mask],
      breaks = BIN_EDGES,
      labels = BIN_LABELS,
      include.lowest = TRUE, right = TRUE))

cat("\npre-samples by lead-time bin (train):\n")
print(table(trn$leadtime_bin, useNA = "ifany"))

# --- 1c. STABLE reference: same subjects, no upcoming episode
#         within 400 days (or subjects with no episodes at all)
trn$is_stable <- (is.na(trn$days_to_next_onset) |
                    trn$days_to_next_onset > 400) &
  (is.na(trn$episode_phase) |
     trn$episode_phase == "post")   # post is also acceptable ref
# remove if you prefer purer

cat("\nstable reference samples (train):", sum(trn$is_stable), "\n")
cat("from", length(unique(trn$subject_id[trn$is_stable])), "subjects\n")


# ============================================================
# 2. DESCRIPTIVE TABLE — mean ± SD per bin + stable
# ============================================================

outcomes <- c("median_abs_z", "lag1_sd", "rmssd", "mahal_dist")
outcome_labels <- c("Median |Z|", "Lag-1 SD/day", "RMSSD/day",
                    "Mahalanobis dist")

groups <- c(BIN_LABELS, "stable")

desc <- data.frame()
for (g in groups) {
  if (g == "stable") {
    rows <- trn[trn$is_stable, ]
  } else {
    rows <- trn[!is.na(trn$leadtime_bin) & trn$leadtime_bin == g, ]
  }
  for (oc in outcomes) {
    v <- rows[[oc]]
    v <- v[!is.na(v)]
    desc <- rbind(desc, data.frame(
      group    = g,
      outcome  = oc,
      n        = length(v),
      n_subj   = length(unique(rows$subject_id[!is.na(rows[[oc]])])),
      mean     = if (length(v) > 0) mean(v)   else NA,
      sd       = if (length(v) > 0) sd(v)     else NA,
      median   = if (length(v) > 0) median(v) else NA,
      q25      = if (length(v) > 0) quantile(v, 0.25) else NA,
      q75      = if (length(v) > 0) quantile(v, 0.75) else NA,
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n=== DESCRIPTIVE TABLE (train) ===\n")
print(desc[desc$outcome == "median_abs_z",
           c("group","n","n_subj","mean","sd","median","q25","q75")],
      row.names = FALSE, digits = 3)


# ============================================================
# 3. WITHIN-SUBJECT PERMUTATION TEST — PER BIN
# ============================================================
#
# For each lead-time bin:
#   - identify subjects with >= 1 sample in this bin AND >= 1
#     stable sample
#   - observed statistic: mean(metric[in bin]) - mean(metric[stable])
#     computed within each eligible subject, then averaged
#   - permute: within each subject, randomly reassign which of
#     (bin + stable) samples get the "bin" label; recompute stat
#   - p-value: proportion of null >= observed (one-tailed)

cat("\n\n=== PER-BIN PERMUTATION TESTS ===\n")

perm_results <- list()
all_nulls    <- list()

for (bn in BIN_LABELS) {
  
  cat("\n--- bin:", bn, "---\n")
  
  # subjects who have at least one sample in this bin AND one stable
  bin_subjs    <- unique(trn$subject_id[!is.na(trn$leadtime_bin) &
                                          trn$leadtime_bin == bn])
  stable_subjs <- unique(trn$subject_id[trn$is_stable])
  elig_subjs   <- intersect(bin_subjs, stable_subjs)
  
  cat("  eligible subjects (bin AND stable):", length(elig_subjs), "\n")
  
  n_bin_samp <- sum(!is.na(trn$leadtime_bin) & trn$leadtime_bin == bn &
                      trn$subject_id %in% elig_subjs)
  cat("  bin samples from eligible subjects:", n_bin_samp, "\n")
  
  if (length(elig_subjs) < MIN_PRE_SUBJECTS ||
      n_bin_samp < MIN_PRE_SAMPLES) {
    cat("  SKIP: below minimum thresholds\n")
    perm_results[[bn]] <- data.frame(
      bin = bn, n_pre = n_bin_samp, n_subj = length(elig_subjs),
      obs_medZ = NA, p_medZ = NA,
      obs_lag1sd = NA, p_lag1sd = NA,
      obs_mahal = NA, p_mahal = NA,
      stringsAsFactors = FALSE)
    next
  }
  
  # working data: bin samples + stable samples, eligible subjects only
  wd <- trn[trn$subject_id %in% elig_subjs &
              ((!is.na(trn$leadtime_bin) & trn$leadtime_bin == bn) |
                 trn$is_stable), ]
  wd$in_bin <- !is.na(wd$leadtime_bin) & wd$leadtime_bin == bn
  
  # observed within-subject mean differences
  obs_stat <- function(data, label_col, metrics) {
    sapply(metrics, function(oc) {
      diffs <- sapply(elig_subjs, function(s) {
        d <- data[data$subject_id == s, ]
        v_bin    <- d[[oc]][d[[label_col]]]
        v_stable <- d[[oc]][!d[[label_col]]]
        v_bin    <- v_bin[!is.na(v_bin)]
        v_stable <- v_stable[!is.na(v_stable)]
        if (length(v_bin) == 0 || length(v_stable) == 0) return(NA_real_)
        mean(v_bin) - mean(v_stable)
      })
      mean(diffs, na.rm = TRUE)
    })
  }
  
  metrics <- c("median_abs_z", "lag1_sd", "mahal_dist")
  obs <- obs_stat(wd, "in_bin", metrics)
  names(obs) <- metrics
  cat("  observed diffs (bin - stable):\n")
  for (m in metrics) cat("   ", m, ":", round(obs[m], 4), "\n")
  
  # permutation loop
  per_subj_idx <- split(seq_len(nrow(wd)), wd$subject_id)
  n_bin_per_subj <- sapply(elig_subjs, function(s)
    sum(wd$in_bin[wd$subject_id == s]))
  
  null_mat <- matrix(NA_real_, nrow = N_PERM, ncol = length(metrics))
  colnames(null_mat) <- metrics
  
  for (b in seq_len(N_PERM)) {
    perm_label <- wd$in_bin
    for (s in elig_subjs) {
      idx_s <- per_subj_idx[[s]]
      perm_label[idx_s] <- FALSE
      n_s <- n_bin_per_subj[s]
      if (n_s > 0 && n_s < length(idx_s))
        perm_label[idx_s[sample(length(idx_s), n_s)]] <- TRUE
    }
    null_mat[b, ] <- sapply(metrics, function(oc) {
      diffs <- sapply(elig_subjs, function(s) {
        d  <- wd[wd$subject_id == s, ]
        pl <- perm_label[wd$subject_id == s]
        v1 <- d[[oc]][pl];  v1 <- v1[!is.na(v1)]
        v0 <- d[[oc]][!pl]; v0 <- v0[!is.na(v0)]
        if (length(v1) == 0 || length(v0) == 0) return(NA_real_)
        mean(v1) - mean(v0)
      })
      mean(diffs, na.rm = TRUE)
    })
  }
  
  p_vals <- colMeans(null_mat >= matrix(obs, nrow = N_PERM,
                                        ncol = length(metrics),
                                        byrow = TRUE), na.rm = TRUE)
  
  cat("  p-values (one-tailed, bin > stable):\n")
  for (m in metrics) cat("   ", m, ": p =", round(p_vals[m], 4), "\n")
  
  perm_results[[bn]] <- data.frame(
    bin        = bn,
    n_pre      = n_bin_samp,
    n_subj     = length(elig_subjs),
    obs_medZ   = obs["median_abs_z"],
    p_medZ     = p_vals["median_abs_z"],
    obs_lag1sd = obs["lag1_sd"],
    p_lag1sd   = p_vals["lag1_sd"],
    obs_mahal  = obs["mahal_dist"],
    p_mahal    = p_vals["mahal_dist"],
    stringsAsFactors = FALSE, row.names = NULL
  )
  all_nulls[[bn]] <- null_mat
}

results <- do.call(rbind, perm_results)

# Benjamini-Hochberg across the 5 bins
for (pc in c("p_medZ", "p_lag1sd", "p_mahal")) {
  results[[paste0(pc, "_bh")]] <- p.adjust(results[[pc]], method = "BH")
}

cat("\n\n=== RESULTS TABLE ===\n")
print(results[, c("bin","n_pre","n_subj",
                  "obs_medZ","p_medZ","p_medZ_bh",
                  "obs_lag1sd","p_lag1sd",
                  "obs_mahal","p_mahal")],
      row.names = FALSE, digits = 3)


# ============================================================
# 4. EFFECT SIZE — Cohen's d vs stable
# ============================================================
#
# Pooled Cohen's d (not within-subject), reported alongside
# the permutation p-value.  Gives readers a sense of magnitude.

cat("\n=== EFFECT SIZES (Cohen's d, median_abs_z vs stable) ===\n")

stable_z <- trn$median_abs_z[trn$is_stable & !is.na(trn$median_abs_z)]
sd_stable <- sd(stable_z)

results$cohens_d_medZ <- NA_real_
for (i in seq_len(nrow(results))) {
  bn <- results$bin[i]
  v  <- trn$median_abs_z[!is.na(trn$leadtime_bin) &
                           trn$leadtime_bin == bn &
                           !is.na(trn$median_abs_z)]
  if (length(v) < 2 || is.na(sd_stable) || sd_stable == 0) next
  sd_pool <- sqrt((sd(v)^2 + sd_stable^2) / 2)
  results$cohens_d_medZ[i] <- (mean(v) - mean(stable_z)) / sd_pool
}
print(results[, c("bin","obs_medZ","cohens_d_medZ","p_medZ","p_medZ_bh")],
      row.names = FALSE, digits = 3)


# ============================================================
# 5. TEST SET CONFIRMATION (look but don't interpret)
# ============================================================

tst <- samp[!is.na(samp$split) & samp$split == "test", ]
tst$leadtime_bin <- NA_character_
pre_mask_t <- !is.na(tst$days_to_next_onset) &
  tst$days_to_next_onset >= BIN_EDGES[1] &
  tst$days_to_next_onset <= BIN_EDGES[length(BIN_EDGES)]
tst$leadtime_bin[pre_mask_t] <- as.character(
  cut(tst$days_to_next_onset[pre_mask_t],
      breaks = BIN_EDGES, labels = BIN_LABELS,
      include.lowest = TRUE, right = TRUE))
tst$is_stable <- (is.na(tst$days_to_next_onset) |
                    tst$days_to_next_onset > 400) &
  (is.na(tst$episode_phase) | tst$episode_phase == "post")

cat("\n\n=== TEST SET: descriptive only (DO NOT use to select) ===\n")
test_desc <- data.frame()
for (g in c(BIN_LABELS, "stable")) {
  rows <- if (g == "stable") tst[tst$is_stable, ] else
    tst[!is.na(tst$leadtime_bin) & tst$leadtime_bin == g, ]
  v <- rows$median_abs_z; v <- v[!is.na(v)]
  test_desc <- rbind(test_desc, data.frame(
    group = g, n = length(v),
    mean  = if (length(v) > 0) round(mean(v), 3) else NA,
    sd    = if (length(v) > 0) round(sd(v), 3)   else NA,
    stringsAsFactors = FALSE))
}
print(test_desc, row.names = FALSE)


# ============================================================
# 6. MAIN FIGURE
# ============================================================
#
# Panel A: median |Z| by lead-time bin, with stable reference.
#          Points are subject-level means; box shows distribution.
# Panel B: permutation nulls for the two most informative bins
#          (closest to onset: 0-45 d and 45-75 d).

png("outputs/Q1_main_figure.png", width = 1100, height = 550)
par(mfrow = c(1, 2), mar = c(6, 5, 4, 2))

# --- Panel A: boxplot per group ---
groups_all <- c(BIN_LABELS, "stable")
plot_data <- lapply(groups_all, function(g) {
  if (g == "stable")
    trn$median_abs_z[trn$is_stable & !is.na(trn$median_abs_z)]
  else
    trn$median_abs_z[!is.na(trn$leadtime_bin) &
                       trn$leadtime_bin == g &
                       !is.na(trn$median_abs_z)]
})
names(plot_data) <- groups_all

n_per <- sapply(plot_data, length)
cols  <- c(rep("#E64B35", length(BIN_LABELS)), "#AAAAAA")

bp <- boxplot(plot_data,
              col     = cols,
              outline = FALSE,
              las     = 2,
              ylab    = "Median |Z-score from personal baseline|",
              main    = "Q1: metabolome deviation before dysbiosis onset\n(train set)",
              xaxt    = "n",
              cex.lab = 0.9)

# x-axis labels with n
axis(1, at = seq_along(groups_all),
     labels = paste0(groups_all, "\nn=", n_per),
     las = 2, cex.axis = 0.75, tick = FALSE)

# add significance markers from BH-corrected p-values
for (i in seq_len(nrow(results))) {
  p <- results$p_medZ_bh[i]
  if (!is.na(p)) {
    sym <- if (p < 0.001) "***" else if (p < 0.01) "**" else
      if (p < 0.05) "*" else if (p < 0.1) "." else ""
    if (sym != "") {
      y_top <- max(unlist(plot_data[i]), na.rm = TRUE) * 1.05
      text(i, y_top, sym, cex = 1.2, col = "red")
    }
  }
}

# reference line at stable median
abline(h = median(plot_data[["stable"]], na.rm = TRUE),
       col = "grey50", lty = 2, lwd = 1.5)
legend("topright", legend = c("pre-onset bins", "stable reference",
                              "stable median"),
       fill = c("#E64B35", "#AAAAAA", NA),
       lty  = c(NA, NA, 2), lwd = c(NA, NA, 1.5),
       col  = c(NA, NA, "grey50"), bty = "n", cex = 0.8,
       border = c("black","black", NA))


# --- Panel B: permutation nulls for the nearest two bins ---
plot_bins <- BIN_LABELS[1:min(2, length(BIN_LABELS))]
null_cols  <- c("#E64B35", "#4DBBD5")

xlims <- range(unlist(lapply(plot_bins, function(bn)
  if (!is.null(all_nulls[[bn]])) all_nulls[[bn]][, "median_abs_z"])), na.rm = TRUE)
xlims <- c(min(xlims, -0.3), max(xlims, 0.3))

plot(NULL, xlim = xlims,
     ylim = c(0, N_PERM * 0.12),
     xlab = "null mean difference (bin − stable)",
     ylab = "frequency",
     main = "Permutation null distributions\n(closest lead-time bins)")

for (bi in seq_along(plot_bins)) {
  bn <- plot_bins[bi]
  if (is.null(all_nulls[[bn]])) next
  null_v <- all_nulls[[bn]][, "median_abs_z"]
  h <- hist(null_v, breaks = 40, plot = FALSE)
  for (j in seq_along(h$mids)) {
    rect(h$breaks[j], 0, h$breaks[j + 1], h$counts[j],
         col = adjustcolor(null_cols[bi], 0.4),
         border = NA)
  }
  obs_v <- results$obs_medZ[results$bin == bn]
  if (!is.na(obs_v)) {
    abline(v = obs_v, col = null_cols[bi], lwd = 2, lty = 1)
    p_v <- results$p_medZ[results$bin == bn]
    text(obs_v, N_PERM * 0.11,
         labels = paste0(bn, "\nobs=", round(obs_v, 3),
                         "\np=",  round(p_v, 3)),
         col = null_cols[bi], cex = 0.75, adj = 0)
  }
}
abline(v = 0, col = "black", lty = 2)
legend("topleft", legend = plot_bins, fill = adjustcolor(null_cols, 0.4),
       border = NA, bty = "n", cex = 0.8)

dev.off()
cat("saved outputs/Q1_main_figure.png\n")


# ============================================================
# 7. SUPPLEMENTARY FIGURE — all bins + secondary outcomes
# ============================================================

png("outputs/Q1_supplementary.png", width = 1200, height = 900)
par(mfrow = c(2, 3), mar = c(7, 5, 4, 2))

# (a-c) boxplots for lag1_sd, rmssd, mahal_dist
for (oc in c("lag1_sd", "rmssd", "mahal_dist")) {
  pd <- lapply(groups_all, function(g) {
    if (g == "stable") trn[[oc]][trn$is_stable & !is.na(trn[[oc]])]
    else trn[[oc]][!is.na(trn$leadtime_bin) &
                     trn$leadtime_bin == g & !is.na(trn[[oc]])]
  })
  names(pd) <- groups_all
  ns <- sapply(pd, length)
  boxplot(pd, col = cols, outline = FALSE, las = 2,
          ylab = oc,
          main = paste("Q1:", oc, "(train)"),
          xaxt = "n")
  axis(1, at = seq_along(groups_all),
       labels = paste0(groups_all, "\nn=", ns),
       las = 2, cex.axis = 0.7, tick = FALSE)
  abline(h = median(pd[["stable"]], na.rm = TRUE),
         col = "grey50", lty = 2)
}

# (d) effect sizes across bins
bn_ok <- !is.na(results$cohens_d_medZ)
if (sum(bn_ok) > 0) {
  barplot(results$cohens_d_medZ[bn_ok],
          names.arg = results$bin[bn_ok],
          col = "#E64B35", las = 2,
          ylab = "Cohen's d vs stable",
          main = "Effect size (median |Z|)",
          cex.names = 0.8)
  abline(h = 0, lty = 2)
  abline(h = c(-0.2, 0.2), lty = 3, col = "grey60")
  legend("topright", legend = "small effect = 0.2",
         lty = 3, col = "grey60", bty = "n", cex = 0.8)
}

# (e) n samples contributing to each bin
barplot(c(sapply(BIN_LABELS, function(bn)
  sum(!is.na(trn$leadtime_bin) & trn$leadtime_bin == bn)),
  sum(trn$is_stable)),
  names.arg = c(BIN_LABELS, "stable"),
  col = c(rep("#E64B35", length(BIN_LABELS)), "#AAAAAA"),
  las = 2, ylab = "n samples",
  main = "Sample counts per group\n(train)",
  cex.names = 0.75)

# (f) p-value profile across bins
p_ok <- !is.na(results$p_medZ)
if (sum(p_ok) > 1) {
  plot(seq_len(sum(p_ok)), results$p_medZ[p_ok],
       type = "b", pch = 16, col = "#E64B35",
       xlab = "lead-time bin (furthest → closest)",
       ylab = "permutation p-value (one-tailed)",
       main = "p-value profile across bins\n(median |Z|, train)",
       xaxt = "n", ylim = c(0, 1))
  axis(1, at = seq_len(sum(p_ok)),
       labels = results$bin[p_ok], las = 2, cex.axis = 0.75)
  abline(h = 0.05, lty = 2, col = "grey50")
  lines(seq_len(sum(p_ok)), results$p_medZ_bh[p_ok],
        type = "b", pch = 17, col = "#4DBBD5", lty = 2)
  legend("topright",
         legend = c("nominal p", "BH-adjusted"),
         col = c("#E64B35","#4DBBD5"),
         pch = c(16, 17), lty = c(1, 2), bty = "n", cex = 0.85)
}

dev.off()
cat("saved outputs/Q1_supplementary.png\n")


# ============================================================
# 8. SAVE RESULTS
# ============================================================

saveRDS(all_nulls, "data/derived/Q1_permutation_bins.rds")
write.csv(results, "outputs/Q1_results.csv", row.names = FALSE)
write.csv(desc,    "outputs/Q1_descriptive.csv", row.names = FALSE)

cat("\nwrote outputs/Q1_results.csv\n")
cat("wrote outputs/Q1_descriptive.csv\n")


# ============================================================
# 9. METHODS NUMBERS  — everything cited in the text
# ============================================================

sink("outputs/Q1_methods_numbers.txt")
cat("Q1 METHODS NUMBERS\n")
cat("generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("TRAIN SET\n")
cat("  total train samples:", nrow(trn), "\n")
cat("  stable reference samples:", sum(trn$is_stable), "\n")
cat("  stable reference subjects:",
    length(unique(trn$subject_id[trn$is_stable])), "\n\n")

cat("PRE-ONSET SAMPLES BY BIN\n")
for (bn in BIN_LABELS) {
  ns <- sum(!is.na(trn$leadtime_bin) & trn$leadtime_bin == bn)
  nsubj <- length(unique(trn$subject_id[
    !is.na(trn$leadtime_bin) & trn$leadtime_bin == bn]))
  cat(sprintf("  %-12s : %3d samples, %2d subjects\n", bn, ns, nsubj))
}

cat("\nPERMUTATION TEST (N =", N_PERM, ")\n")
cat("  outcome: median |Z-score from personal baseline|\n")
cat("  direction: one-tailed (pre > stable)\n")
cat("  correction: Benjamini-Hochberg across 5 bins\n\n")

cat("RESULTS\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("  bin %-12s : n=%2d, d=%.3f, p=%.3f, p_BH=%.3f\n",
              results$bin[i], results$n_pre[i],
              ifelse(is.na(results$cohens_d_medZ[i]), NA,
                     results$cohens_d_medZ[i]),
              ifelse(is.na(results$p_medZ[i]), NA, results$p_medZ[i]),
              ifelse(is.na(results$p_medZ_bh[i]), NA,
                     results$p_medZ_bh[i])))
}
sink()

cat("wrote outputs/Q1_methods_numbers.txt\n")


# ============================================================
# 10. INTERPRETATION GUIDE
# ============================================================
#
# Reading the results:
#
# If the 0–45 d or 45–75 d bins show p < 0.05 (permutation)
# AND Cohen's d >= 0.2 AND the BH-corrected p < 0.1:
#   → "gut chemistry becomes more erratic in the ~1–2 months
#     preceding dysbiosis onset"
#   → This feeds directly into Q2 (is deviation a better
#     predictor than absolute levels?) and Q3 (at which
#     lead time is the signal strongest?).
#
# If only the 60–80 d region is elevated (as the scatter plot
# in Person C suggested):
#   → "signal is not monotone — it peaks ~2 months before
#     onset, not in the final weeks"
#   → This is a real finding, not a negative result.
#     It is consistent with a prodromal phase rather than
#     an acute precursor.
#
# If no bin reaches p < 0.05:
#   → Q1 is a negative result WITH the correct N (~40 onsets,
#     40 pre-samples spread across 5 bins means ~8/bin).
#   → Report the effect sizes and CIs, not just p-values.
#   → The scatter plot shape from Person C is still reportable
#     as a visual observation.
#
# The permutation test is correctly powered for the data you
# have.  Do NOT reframe a negative result as "insufficient
# power" without reporting the observed effect sizes.
#
# ============================================================

cat("\n=== DONE — Q1 complete ===\n")
cat("Main finding is in outputs/Q1_results.csv\n")
cat("Read outputs/Q1_methods_numbers.txt before writing up.\n")