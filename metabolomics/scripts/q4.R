setwd("/home/subhrab/dh607_project/metabolomics")
suppressPackageStartupMessages(library(arrow))

samp <- readRDS("data/derived/volatility.rds")
dev  <- readRDS("data/derived/deviation_matrix.rds")
fa   <- read_parquet("data/derived/features_full_v1.parquet")
fi   <- read_parquet("data/derived/features_full_v1_info.parquet")

feat_ids   <- fa$feature_id
sample_ids <- setdiff(colnames(fa), "feature_id")
ABS_full   <- t(as.matrix(fa[, sample_ids]))
rownames(ABS_full) <- sample_ids
colnames(ABS_full) <- feat_ids

# Rebuild s2 cleanly from scratch (do NOT reuse a filtered s2 from Q2)
common <- Reduce(intersect, list(rownames(dev), rownames(ABS_full), samp$sample_id))
s2     <- samp[match(common, samp$sample_id), ]
D_all  <- dev[match(common, rownames(dev)), , drop = FALSE]
A_all  <- ABS_full[match(common, rownames(ABS_full)), , drop = FALSE]

cat("=== S2 INTEGRITY AUDIT ===\n")
cat("1. Total aligned rows:", nrow(s2), "(expected ~542)\n")
cat("2. Unique sample_ids:  ", length(unique(s2$sample_id)), "\n")
cat("3. Row order exact match (D & A):",
    all(rownames(D_all) == s2$sample_id) && all(rownames(A_all) == s2$sample_id), "\n")
cat("4. Missing subject_id or day:",
    sum(is.na(s2$subject_id) | is.na(s2$day)), "\n")
cat("5. nonIBD rows with episode_phase or days_to_next_onset:\n")
print(table(
  has_phase = !is.na(s2$episode_phase[s2$diagnosis == "nonIBD"]),
  has_onset = !is.na(s2$days_to_next_onset[s2$diagnosis == "nonIBD"])
))
cat("\n6. Cross-tab of episode_phase vs. non-NA days_to_next_onset:\n")
print(table(
  phase     = ifelse(is.na(s2$episode_phase), "<NA>", s2$episode_phase),
  has_onset = !is.na(s2$days_to_next_onset)
))
##############################################################################################################################
# ============================================================
# Q4_partA_setpoint.R
# Part A: Does the IBD metabolome return to its personal
#         pre-flare / baseline set point after a flare ("post")?
# ============================================================

setwd("/home/subhrab/dh607_project/metabolomics")

suppressPackageStartupMessages(library(arrow))

dir.create("outputs", showWarnings = FALSE)
set.seed(607)

# ============================================================
# 1. LOAD & ALIGN
# ============================================================

samp <- readRDS("data/derived/volatility.rds")
dev  <- readRDS("data/derived/deviation_matrix.rds")
fa   <- read_parquet("data/derived/features_full_v1.parquet")
fi   <- read_parquet("data/derived/features_full_v1_info.parquet")

feat_ids   <- fa$feature_id
sample_ids <- setdiff(colnames(fa), "feature_id")
ABS_full   <- t(as.matrix(fa[, sample_ids]))
rownames(ABS_full) <- sample_ids
colnames(ABS_full) <- feat_ids

common <- Reduce(intersect, list(rownames(dev), rownames(ABS_full), samp$sample_id))
s2     <- samp[match(common, samp$sample_id), ]
D_all  <- dev[match(common, rownames(dev)), , drop = FALSE]
A_all  <- ABS_full[match(common, rownames(ABS_full)), , drop = FALSE]

stopifnot(all(rownames(D_all) == s2$sample_id))
stopifnot(all(rownames(A_all) == s2$sample_id))

# Restrict strictly to IBD subjects (CD & UC)
ibd   <- s2$diagnosis %in% c("CD", "UC")
s_ibd <- s2[ibd, ]
A_ibd <- A_all[ibd, , drop = FALSE]

# Use annotated HMDB features in log2 space
ann_ids <- fi$feature_id[fi$has_hmdb & fi$feature_id %in% colnames(A_ibd)]
A_ann   <- A_ibd[, ann_ids, drop = FALSE]

cat("IBD samples:", nrow(s_ibd), "across", length(unique(s_ibd$subject_id)), "subjects\n")
cat("Annotated features used:", length(ann_ids), "\n\n")

# ============================================================
# 2. COMPUTE DISTANCE TO PERSONAL SET POINT
#    Set point = per-subject median across baseline/pre-flare
#    visits (episode_phase is NA or "pre").
#    Metric = root-mean-square log2 shift (RMSD) from set point.
# ============================================================

s_ibd$dist_to_setpoint <- NA_real_

for (subj in unique(s_ibd$subject_id)) {
  idx_all <- which(s_ibd$subject_id == subj)
  idx_ref <- idx_all[is.na(s_ibd$episode_phase[idx_all]) |
                       s_ibd$episode_phase[idx_all] == "pre"]
  if (length(idx_ref) < 2) next
  
  ref_vec <- apply(A_ann[idx_ref, , drop = FALSE], 2, median, na.rm = TRUE)
  
  for (i in idx_all) {
    v  <- A_ann[i, ]
    ok <- !is.na(v) & !is.na(ref_vec)
    if (sum(ok) >= 10) {
      s_ibd$dist_to_setpoint[i] <- sqrt(mean((v[ok] - ref_vec[ok])^2))
    }
  }
}

s_ibd$phase_clean <- ifelse(is.na(s_ibd$episode_phase), "stable", s_ibd$episode_phase)
s_ibd$phase_clean <- factor(s_ibd$phase_clean,
                            levels = c("stable", "pre", "during", "post"))

# ============================================================
# 3. SUMMARY BY PHASE & PAIRED WITHIN-SUBJECT TESTS
# ============================================================

cat("=== PART A: DISTANCE TO PERSONAL SET POINT BY PHASE ===\n")
phase_summary <- do.call(rbind, lapply(levels(s_ibd$phase_clean), function(ph) {
  x <- s_ibd$dist_to_setpoint[s_ibd$phase_clean == ph]
  x <- x[!is.na(x)]
  data.frame(
    phase      = ph,
    n_samples  = length(x),
    n_subjects = length(unique(s_ibd$subject_id[s_ibd$phase_clean == ph &
                                                  !is.na(s_ibd$dist_to_setpoint)])),
    median_rms = round(median(x), 4),
    iqr_rms    = round(IQR(x), 4),
    mean_rms   = round(mean(x), 4),
    stringsAsFactors = FALSE
  )
}))
print(phase_summary, row.names = FALSE)

# Helper for paired within-subject Wilcoxon comparison between two phases
paired_phase_test <- function(df, ph1, ph2) {
  s1 <- unique(df$subject_id[df$phase_clean == ph1 & !is.na(df$dist_to_setpoint)])
  s2_u <- unique(df$subject_id[df$phase_clean == ph2 & !is.na(df$dist_to_setpoint)])
  both <- intersect(s1, s2_u)
  if (length(both) < 3) {
    cat(sprintf("\nPaired %s vs %s: too few overlapping subjects (n=%d)\n",
                ph1, ph2, length(both)))
    return(NULL)
  }
  v1 <- vapply(both, function(s)
    mean(df$dist_to_setpoint[df$subject_id == s & df$phase_clean == ph1], na.rm = TRUE),
    numeric(1))
  v2 <- vapply(both, function(s)
    mean(df$dist_to_setpoint[df$subject_id == s & df$phase_clean == ph2], na.rm = TRUE),
    numeric(1))
  wt <- suppressWarnings(wilcox.test(v1, v2, paired = TRUE))
  cat(sprintf("\nPaired %s vs %s (n=%d subjects): median %s=%.4f, %s=%.4f, diff=%+.4f, p=%.4f\n",
              ph1, ph2, length(both), ph1, median(v1), ph2, median(v2),
              median(v2 - v1), wt$p.value))
  invisible(list(subjects = both, v1 = v1, v2 = v2, p = wt$p.value))
}

cat("\n=== PAIRED WITHIN-SUBJECT COMPARISONS ===\n")
res_dur_post  <- paired_phase_test(s_ibd, "during", "post")
res_pre_post  <- paired_phase_test(s_ibd, "pre",    "post")
res_stab_post <- paired_phase_test(s_ibd, "stable", "post")

# ============================================================
# 4. FIGURE: PHASE DISTRIBUTION + SUBJECT TRAJECTORIES
# ============================================================

png("outputs/Q4_partA_setpoint.png", width = 1000, height = 480)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Panel A: Boxplot by phase
phase_list <- split(s_ibd$dist_to_setpoint, s_ibd$phase_clean)
phase_list <- lapply(phase_list, function(x) x[!is.na(x)])
ns <- vapply(phase_list, length, integer(1))

boxplot(phase_list,
        col     = c("#AAAAAA", "olivedrab", "#E64B35", "#4DBBD5"),
        outline = FALSE,
        ylab    = "RMS log2 distance from personal set point",
        main    = "IBD Metabolome Distance\nto Personal Set Point by Phase",
        xaxt    = "n")
axis(1, at = 1:4, labels = paste0(names(phase_list), "\n(n=", ns, ")"),
     tick = FALSE, cex.axis = 0.85)

# Panel B: Subject-level paired means across phases
subj_means <- aggregate(dist_to_setpoint ~ subject_id + phase_clean,
                        data = s_ibd, FUN = mean, na.rm = TRUE)
wide_m <- reshape(subj_means, idvar = "subject_id",
                  timevar = "phase_clean", direction = "wide")
cols_ord <- paste0("dist_to_setpoint.", c("stable", "pre", "during", "post"))
cols_ord <- cols_ord[cols_ord %in% colnames(wide_m)]
mat_traj <- as.matrix(wide_m[, cols_ord, drop = FALSE])

plot(1, type = "n", xlim = c(0.8, 4.2),
     ylim = range(mat_traj, na.rm = TRUE),
     xaxt = "n", xlab = "",
     ylab = "Subject mean RMS log2 distance",
     main = "Within-Subject Trajectories\nAcross Flare Phases")
axis(1, at = 1:4, labels = c("stable", "pre", "during", "post"))

for (r in seq_len(nrow(mat_traj))) {
  y_r <- mat_traj[r, ]
  ok  <- !is.na(y_r)
  if (sum(ok) >= 2) {
    lines((1:4)[ok], y_r[ok], col = adjustcolor("grey40", 0.35), lwd = 1)
    points((1:4)[ok], y_r[ok], pch = 16, cex = 0.7,
           col = adjustcolor("grey20", 0.5))
  }
}
med_traj <- apply(mat_traj, 2, median, na.rm = TRUE)
lines(1:4, med_traj, col = "#E64B35", lwd = 2.5)
points(1:4, med_traj, col = "#E64B35", pch = 18, cex = 1.6)

dev.off()
cat("\nSaved outputs/Q4_partA_setpoint.png\n")

write.csv(phase_summary, "outputs/Q4_partA_summary.csv", row.names = FALSE)
cat("Saved outputs/Q4_partA_summary.csv\n")

