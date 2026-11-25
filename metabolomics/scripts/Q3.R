# ============================================================
# Q3_leadtime.R
# How far ahead can you see the metabolomic signal coming?
# ============================================================
#
# RESTART R BEFORE RUNNING.
#
# Needs:
#   data/derived/volatility.rds
#   data/derived/deviation_matrix.rds
#   data/derived/features_full_v1.parquet
#   data/derived/features_full_v1_info.parquet
#   outputs/Q1_results.csv               (Q1 bin definitions)
#
# Produces:
#   outputs/Q3_results.csv
#   outputs/Q3_main_figure.png
#   outputs/Q3_supplementary.png
#   outputs/Q3_methods_numbers.txt
#
# ============================================================
# QUESTION REFRAME (from feasibility census)
# ─────────────────────────────────────────────────────────────
# The published two-week plateau in metabolite similarity sits
# BELOW the 59.5-day sampling resolution — we cannot observe it.
# Q3 becomes: does the predictive signal (elevated deviation,
# found in Q1 at 45-75 days) persist at longer lead times, or
# is it confined to that window?
#
# APPROACH
# ─────────────────────────────────────────────────────────────
# Primary metric: median |Z-score| (same as Q1 — the metric
#   that showed signal).
# Secondary metric: ABS ridge predicted probability (from Q2 —
#   the metric that classified better).
#
# For each lead-time bin:
#   1. Within-subject permutation test (same engine as Q1)
#      — does the bin differ from stable?
#   2. Effect size (Cohen's d vs stable)
#   3. Trend test across bins — is there a monotone decline
#      from the 45-75 d peak outward?
#
# REFRAMED Q3 ANSWER:
#   - If signal is CONFINED to 45-75 d: narrow prodromal window
#   - If signal PERSISTS across 75-200 d: extended early warning
#   - If signal appears at ALL bins equally: not lead-time specific
#
# ADDITIONAL ANALYSIS: ABS signal by lead time
#   The Q2 ABS model (AUC 0.638) used all pre-onset samples.
#   Q3 asks: at which lead-time window does the ABS signal
#   peak? This tells you the optimal prediction horizon.
# ============================================================

setwd("/home/subhrab/dh607_project/metabolomics")

suppressPackageStartupMessages({
  library(arrow)
  library(glmnet)
})

dir.create("outputs", showWarnings = FALSE)
set.seed(607)

BIN_EDGES  <- c(0, 45, 75, 120, 200, 400)
BIN_LABELS <- c("0–45 d", "45–75 d", "75–120 d", "120–200 d", "200–400 d")
N_PERM     <- 5000
N_BOOT     <- 2000
TOP_SCREEN <- 50
MAD_FLOOR_PCTILE <- 5
MIN_NONZERO_MAD  <- 0.01


# ============================================================
# 0. LOAD
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

fi  <- read_parquet("data/derived/features_full_v1_info.parquet")
q1r <- read.csv("outputs/Q1_results.csv",
                stringsAsFactors = FALSE, check.names = FALSE)

cat("done\n\n")


# ============================================================
# 1. ALIGN + ASSIGN BINS
# ============================================================

common  <- Reduce(intersect, list(rownames(dev),
                                  rownames(ABS_full),
                                  samp$sample_id))
s2      <- samp[match(common, samp$sample_id), ]
DEV_all <- dev[match(common, rownames(dev)),        , drop = FALSE]
ABS_all <- ABS_full[match(common, rownames(ABS_full)), , drop = FALSE]

# assign lead-time bins to ALL samples (not just labelled)
s2$leadtime_bin <- NA_character_
bin_mask <- !is.na(s2$days_to_next_onset) &
  s2$days_to_next_onset >= BIN_EDGES[1] &
  s2$days_to_next_onset <= BIN_EDGES[length(BIN_EDGES)]
s2$leadtime_bin[bin_mask] <- as.character(
  cut(s2$days_to_next_onset[bin_mask],
      breaks = BIN_EDGES, labels = BIN_LABELS,
      include.lowest = TRUE, right = TRUE))

# stable: NA onset, NA phase
s2$is_stable <- is.na(s2$days_to_next_onset) &
  is.na(s2$episode_phase)

cat("lead-time bin counts:\n")
print(table(s2$leadtime_bin, useNA = "ifany"))
cat("stable:", sum(s2$is_stable), "\n\n")

# train mask (no outcome needed — Q3 uses permutation not classifier)
tr <- s2$split == "train" & !is.na(s2$split)


# ============================================================
# 2. ANNOTATED FEATURES + MAD CLIPPING
# ============================================================

ann_ids <- fi$feature_id[fi$has_hmdb &
                           fi$feature_id %in% colnames(DEV_all) &
                           fi$feature_id %in% colnames(ABS_all)]

DEV_ann <- DEV_all[, ann_ids, drop = FALSE]
ABS_ann <- ABS_all[, ann_ids, drop = FALSE]

# MAD clip on train
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
cat("features after MAD clip:", length(keep_feats), "\n\n")

DEV_clean <- DEV_ann[, keep_feats, drop = FALSE]
ABS_clean <- ABS_ann[, keep_feats, drop = FALSE]


# ============================================================
# 3. HELPERS
# ============================================================

auc_roc <- function(labels, scores) {
  ok <- !is.na(scores) & !is.na(labels)
  if (sum(ok) < 4) return(NA_real_)
  labels <- as.integer(labels[ok]); scores <- scores[ok]
  if (length(unique(labels)) < 2) return(NA_real_)
  ord <- order(scores, decreasing = TRUE)
  lab <- labels[ord]; n1 <- sum(lab==1); n0 <- sum(lab==0)
  tp  <- cumsum(lab==1); fp <- cumsum(lab==0)
  sum(diff(fp/n0) * (tp[-1]+tp[-length(tp)]) / 2 / n1)
}

cohens_d <- function(v1, v0) {
  v1 <- v1[!is.na(v1)]; v0 <- v0[!is.na(v0)]
  if (length(v1) < 2 || length(v0) < 2) return(NA_real_)
  sp <- sqrt((var(v1) + var(v0)) / 2)
  if (sp == 0) return(NA_real_)
  (mean(v1) - mean(v0)) / sp
}


# ============================================================
# 4. PRIMARY ANALYSIS: median |Z| BY LEAD-TIME BIN
#    Within-subject permutation test (same engine as Q1)
#    Extended: tests ALL five bins, with trend analysis
# ============================================================

cat("=== PRIMARY: median |Z| permutation by lead-time bin ===\n\n")

# stable reference values (train)
stable_z_tr <- s2$median_abs_z[tr & s2$is_stable & !is.na(s2$median_abs_z)]
cat("stable reference n:", length(stable_z_tr), "\n")
cat("stable median |Z|:", round(median(stable_z_tr), 3), "\n\n")

perm_results <- list()

for (bn in BIN_LABELS) {
  
  # eligible subjects: have >= 1 sample in this bin AND >= 1 stable
  # (train only)
  bin_subjs    <- unique(s2$subject_id[tr & !is.na(s2$leadtime_bin) &
                                         s2$leadtime_bin == bn])
  stable_subjs <- unique(s2$subject_id[tr & s2$is_stable])
  elig_subjs   <- intersect(bin_subjs, stable_subjs)
  
  n_bin  <- sum(tr & !is.na(s2$leadtime_bin) & s2$leadtime_bin == bn &
                  s2$subject_id %in% elig_subjs & !is.na(s2$median_abs_z))
  
  cat(sprintf("bin %-12s : %2d eligible subjects, %2d bin samples\n",
              bn, length(elig_subjs), n_bin))
  
  if (length(elig_subjs) < 3 || n_bin < 3) {
    cat("  SKIP\n")
    perm_results[[bn]] <- data.frame(
      bin=bn, n_bin=n_bin, n_subj=length(elig_subjs),
      obs_diff=NA, p_perm=NA, cohens_d=NA,
      stringsAsFactors=FALSE)
    next
  }
  
  # working data
  wd <- s2[tr & s2$subject_id %in% elig_subjs &
             ((!is.na(s2$leadtime_bin) & s2$leadtime_bin == bn) |
                s2$is_stable) &
             !is.na(s2$median_abs_z), ]
  wd$in_bin <- !is.na(wd$leadtime_bin) & wd$leadtime_bin == bn
  
  # observed within-subject mean diff
  obs_diffs <- sapply(elig_subjs, function(s) {
    d  <- wd[wd$subject_id == s, ]
    v1 <- d$median_abs_z[d$in_bin];  v1 <- v1[!is.na(v1)]
    v0 <- d$median_abs_z[!d$in_bin]; v0 <- v0[!is.na(v0)]
    if (length(v1)==0 || length(v0)==0) return(NA_real_)
    mean(v1) - mean(v0)
  })
  obs_stat <- mean(obs_diffs, na.rm=TRUE)
  
  # effect size (pooled Cohen's d)
  v_bin    <- wd$median_abs_z[wd$in_bin]
  v_stable <- wd$median_abs_z[!wd$in_bin]
  d_cohen  <- cohens_d(v_bin, v_stable)
  
  # permutation
  per_subj     <- split(seq_len(nrow(wd)), wd$subject_id)
  n_bin_per_s  <- sapply(elig_subjs, function(s)
    sum(wd$in_bin[wd$subject_id == s]))
  
  null_vec <- numeric(N_PERM)
  for (b in seq_len(N_PERM)) {
    perm_lbl <- wd$in_bin
    for (s in elig_subjs) {
      idx_s <- per_subj[[s]]
      perm_lbl[idx_s] <- FALSE
      n_s <- n_bin_per_s[s]
      if (n_s > 0 && n_s < length(idx_s))
        perm_lbl[idx_s[sample(length(idx_s), n_s)]] <- TRUE
    }
    diffs <- sapply(elig_subjs, function(s) {
      d  <- wd[wd$subject_id == s, ]
      pl <- perm_lbl[wd$subject_id == s]
      v1 <- d$median_abs_z[pl];  v1 <- v1[!is.na(v1)]
      v0 <- d$median_abs_z[!pl]; v0 <- v0[!is.na(v0)]
      if (length(v1)==0 || length(v0)==0) return(NA_real_)
      mean(v1) - mean(v0)
    })
    null_vec[b] <- mean(diffs, na.rm=TRUE)
  }
  
  p_val <- mean(null_vec >= obs_stat)
  
  cat(sprintf("  obs_diff=%+.3f  d=%+.3f  p=%.4f%s\n",
              obs_stat, d_cohen, p_val,
              ifelse(p_val < 0.05, " *",
                     ifelse(p_val < 0.1, " .", ""))))
  
  perm_results[[bn]] <- data.frame(
    bin=bn, n_bin=n_bin, n_subj=length(elig_subjs),
    obs_diff=round(obs_stat,3), p_perm=round(p_val,4),
    cohens_d=round(d_cohen,3),
    stringsAsFactors=FALSE)
}

res_perm <- do.call(rbind, perm_results)
res_perm$p_bh <- p.adjust(res_perm$p_perm, method="BH")

cat("\n=== PERMUTATION RESULTS (median |Z|) ===\n")
print(res_perm, row.names=FALSE)


# ============================================================
# 5. TREND TEST
#    Is there a monotone decline from the peak bin outward?
#    Spearman correlation between bin midpoint and obs_diff.
#    Tests whether signal decays with increasing lead time
#    (which would indicate a narrow prodromal window) or
#    persists (broad early warning).
# ============================================================

cat("\n=== TREND TEST ===\n")

bin_mids <- c(22.5, 60, 97.5, 160, 300)   # midpoints of 5 bins
names(bin_mids) <- BIN_LABELS

ok_bins  <- !is.na(res_perm$obs_diff)
if (sum(ok_bins) >= 3) {
  mids   <- bin_mids[res_perm$bin[ok_bins]]
  diffs  <- res_perm$obs_diff[ok_bins]
  
  # test 1: monotone decline from peak (45-75 d)
  # sort by bin midpoint (furthest → closest)
  ord     <- order(mids, decreasing=TRUE)
  sp_cor  <- cor(mids[ord], diffs[ord], method="spearman")
  sp_p    <- cor.test(mids[ord], diffs[ord],
                      method="spearman", exact=FALSE)$p.value
  
  cat(sprintf("Spearman r(lead_time, obs_diff) = %.3f  p = %.4f\n",
              sp_cor, sp_p))
  if (sp_cor > 0.5 && sp_p < 0.1) {
    cat("POSITIVE correlation: signal INCREASES with lead time\n")
    cat("→ Earlier samples show STRONGER deviation (unusual)\n")
  } else if (sp_cor < -0.3) {
    cat("NEGATIVE correlation: signal DECAYS with lead time\n")
    cat("→ Signal is STRONGEST near onset, fades at longer lead times\n")
  } else {
    cat("NO CLEAR TREND: signal does not vary systematically with lead time\n")
    cat("→ Signal is confined to specific window (consistent with Q1)\n")
  }
} else {
  cat("insufficient bins for trend test\n")
  sp_cor <- NA; sp_p <- NA
}


# ============================================================
# 6. SECONDARY: ABS SIGNAL BY LEAD-TIME BIN
#    Use ridge logistic (trained on ALL pre-onset vs stable)
#    to get predicted probabilities. Then stratify those
#    probabilities by lead-time bin.
#    Answers: at which horizon is the ABS classifier most
#    confident?
# ============================================================

cat("\n=== SECONDARY: ABS ridge predicted probability by bin ===\n")

# build labels for ABS model (same as Q2_v2)
pre_mask    <- !is.na(s2$days_to_next_onset) &
  s2$days_to_next_onset <= 400 &
  !is.na(s2$episode_phase) &
  s2$episode_phase == "pre" &
  s2$diagnosis != "nonIBD"
stable_mask <- is.na(s2$days_to_next_onset) & is.na(s2$episode_phase)

s2$outcome <- NA_integer_
s2$outcome[pre_mask]    <- 1L
s2$outcome[stable_mask] <- 0L

lab  <- !is.na(s2$outcome)
y_all <- s2$outcome[lab]
tr_lab <- tr[lab]

# screen and fit on train
auc_abs_screen <- sapply(seq_len(ncol(ABS_clean[lab,])), function(j) {
  v <- ABS_clean[lab,][tr_lab, j]
  a <- auc_roc(y_all[tr_lab], v)
  if (is.na(a)) return(NA_real_)
  max(a, 1 - a)
})
names(auc_abs_screen) <- keep_feats
top_abs <- names(sort(auc_abs_screen[!is.na(auc_abs_screen)],
                      decreasing=TRUE)[seq_len(TOP_SCREEN)])

X_tr <- ABS_clean[lab,][tr_lab,  top_abs, drop=FALSE]
X_all <- ABS_clean[lab,][, top_abs, drop=FALSE]
y_tr  <- y_all[tr_lab]

# impute + scale on train
col_mn  <- colMeans(X_tr, na.rm=TRUE)
col_sd  <- apply(X_tr, 2, sd, na.rm=TRUE); col_sd[col_sd==0] <- 1
for (j in seq_len(ncol(X_tr))) {
  X_tr[is.na(X_tr[,j]),j]   <- col_mn[j]
  X_all[is.na(X_all[,j]),j] <- col_mn[j]
}
X_tr_sc  <- scale(X_tr,  center=col_mn, scale=col_sd)
X_all_sc <- scale(X_all, center=col_mn, scale=col_sd)

cv_fit <- cv.glmnet(X_tr_sc, y_tr, family="binomial",
                    alpha=0, nfolds=5, type.measure="auc")
p_all  <- as.numeric(predict(cv_fit, newx=X_all_sc,
                             s="lambda.min", type="response"))
s2$abs_pred_prob <- NA_real_
s2$abs_pred_prob[lab] <- p_all

# mean predicted probability by bin (train)
cat("mean ABS predicted probability by lead-time bin (train):\n")
abs_by_bin <- data.frame()
for (bn in c(BIN_LABELS, "stable")) {
  if (bn == "stable") {
    rows <- s2[tr & s2$is_stable & !is.na(s2$abs_pred_prob), ]
  } else {
    rows <- s2[tr & !is.na(s2$leadtime_bin) &
                 s2$leadtime_bin == bn &
                 !is.na(s2$abs_pred_prob), ]
  }
  v <- rows$abs_pred_prob
  abs_by_bin <- rbind(abs_by_bin, data.frame(
    bin    = bn,
    n      = length(v),
    mean_p = round(mean(v, na.rm=TRUE), 3),
    sd_p   = round(sd(v, na.rm=TRUE), 3),
    stringsAsFactors=FALSE))
  cat(sprintf("  %-12s : n=%3d  mean_p=%.3f  sd=%.3f\n",
              bn, length(v),
              mean(v, na.rm=TRUE), sd(v, na.rm=TRUE)))
}


# ============================================================
# 7. COMBINED RESULTS TABLE
# ============================================================

res_q3 <- merge(res_perm,
                abs_by_bin[abs_by_bin$bin %in% BIN_LABELS,
                           c("bin","n","mean_p")],
                by="bin", all.x=TRUE)
colnames(res_q3)[colnames(res_q3)=="n"]      <- "n_abs_pred"
colnames(res_q3)[colnames(res_q3)=="mean_p"] <- "abs_mean_pred_p"

cat("\n=== Q3 COMBINED RESULTS ===\n")
print(res_q3[, c("bin","n_bin","n_subj","obs_diff","p_perm","p_bh",
                 "cohens_d","abs_mean_pred_p")],
      row.names=FALSE)


# ============================================================
# 8. FIGURES
# ============================================================

png("outputs/Q3_main_figure.png", width=1100, height=520)
par(mfrow=c(1,3), mar=c(7,5,4,2))

# --- Panel A: median |Z| obs_diff by bin --------------------
ok   <- !is.na(res_perm$obs_diff)
cols <- ifelse(!is.na(res_perm$p_perm) & res_perm$p_perm < 0.05,
               "#E64B35", "#FFAAAA")
cols[!ok] <- "grey80"

bp1 <- barplot(ifelse(ok, res_perm$obs_diff, 0),
               names.arg = paste0(BIN_LABELS, "\n(n=",
                                  res_perm$n_bin, ")"),
               col  = cols, las=2,
               ylim = c(min(c(res_perm$obs_diff, -0.1), na.rm=TRUE) - 0.05,
                        max(c(res_perm$obs_diff,  0.1), na.rm=TRUE) + 0.15),
               ylab = "mean within-subject diff (bin − stable)",
               main = "Q3: Deviation signal by lead time\n(median |Z|, train)",
               cex.names=0.75)
abline(h=0, lty=2, col="grey50")

# p-value markers
for (i in seq_len(nrow(res_perm))) {
  p <- res_perm$p_perm[i]
  if (!is.na(p)) {
    sym <- if (p<0.001) "***" else if (p<0.01) "**" else
      if (p<0.05) "*" else if (p<0.1) "." else ""
    if (sym != "") {
      y_pos <- (res_perm$obs_diff[i] %||% 0) +
        max(abs(res_perm$obs_diff), na.rm=TRUE)*0.08
      text(bp1[i], y_pos, sym, cex=1.3, col="darkred")
    }
  }
}
legend("topright",
       legend=c("p<0.05 (sig)","p≥0.05","no data"),
       fill=c("#E64B35","#FFAAAA","grey80"),
       border="black", bty="n", cex=0.8)

# --- Panel B: ABS predicted probability by bin --------------
abs_plot <- abs_by_bin
stable_p <- abs_plot$mean_p[abs_plot$bin=="stable"]
abs_plot <- abs_plot[abs_plot$bin %in% BIN_LABELS,]
abs_plot <- abs_plot[match(BIN_LABELS, abs_plot$bin),]

cols2 <- ifelse(is.na(abs_plot$mean_p), "grey80",
                ifelse(abs_plot$mean_p > stable_p + 0.05,
                       "#4DBBD5", "#CCE5FF"))

bp2 <- barplot(abs_plot$mean_p,
               names.arg = paste0(BIN_LABELS, "\n(n=", abs_plot$n, ")"),
               col  = cols2, las=2,
               ylim = c(0, max(abs_plot$mean_p, na.rm=TRUE)*1.25),
               ylab = "mean ABS predicted probability",
               main = "Q3: ABS classifier confidence\nby lead-time bin (train)",
               cex.names=0.75)
abline(h=stable_p, lty=2, col="grey50", lwd=1.5)
text(par("usr")[2]*0.98, stable_p+0.01,
     paste0("stable=", round(stable_p,3)),
     col="grey40", cex=0.75, adj=1)

# --- Panel C: p-value profile + Cohen's d -------------------
ok2    <- !is.na(res_perm$p_perm)
x_pos  <- seq_len(sum(ok2))
p_vals <- res_perm$p_perm[ok2]
d_vals <- res_perm$cohens_d[ok2]
xlabs  <- BIN_LABELS[ok2]

plot(x_pos, p_vals,
     type="b", pch=16, col="#E64B35",
     ylim=c(0,1), xaxt="n",
     xlab="lead-time bin (furthest → closest to onset)",
     ylab="permutation p-value",
     main="Q3: p-value and effect size\nprofile across bins")
axis(1, at=x_pos, labels=xlabs, las=2, cex.axis=0.75)
abline(h=0.05, lty=2, col="grey50")

# Cohen's d on secondary axis
par(new=TRUE)
plot(x_pos, d_vals,
     type="b", pch=17, col="#4DBBD5",
     ylim=c(-1,1), axes=FALSE, xlab="", ylab="")
axis(4, col="#4DBBD5", col.axis="#4DBBD5", cex.axis=0.8)
mtext("Cohen's d", side=4, line=2.5, col="#4DBBD5", cex=0.8)
abline(h=0, lty=3, col="#4DBBD5")

legend("topright",
       legend=c("p-value (left)", "Cohen's d (right)"),
       col=c("#E64B35","#4DBBD5"),
       pch=c(16,17), lty=1, bty="n", cex=0.8)

dev.off()
cat("saved outputs/Q3_main_figure.png\n")


# --- Supplementary ------------------------------------------
png("outputs/Q3_supplementary.png", width=1100, height=850)
par(mfrow=c(2,3), mar=c(7,5,4,2))

# (a) boxplots of median|Z| by bin + stable
groups_all <- c(BIN_LABELS, "stable")
plot_data <- lapply(groups_all, function(g) {
  if (g=="stable")
    s2$median_abs_z[tr & s2$is_stable & !is.na(s2$median_abs_z)]
  else
    s2$median_abs_z[tr & !is.na(s2$leadtime_bin) &
                      s2$leadtime_bin==g & !is.na(s2$median_abs_z)]
})
names(plot_data) <- groups_all
ns <- sapply(plot_data, length)
cols3 <- c(rep("#E64B35",5),"#AAAAAA")

boxplot(plot_data, col=cols3, outline=FALSE, las=2,
        ylab="median |Z-score|",
        main="Deviation by lead-time bin\n(train)",
        xaxt="n")
axis(1, at=seq_along(groups_all),
     labels=paste0(groups_all,"\nn=",ns),
     las=2, cex.axis=0.72, tick=FALSE)
abline(h=median(plot_data[["stable"]], na.rm=TRUE),
       col="grey50", lty=2)

# (b) Mahalanobis distance by bin
plot_mah <- lapply(groups_all, function(g) {
  if (g=="stable")
    s2$mahal_dist[tr & s2$is_stable & !is.na(s2$mahal_dist)]
  else
    s2$mahal_dist[tr & !is.na(s2$leadtime_bin) &
                    s2$leadtime_bin==g & !is.na(s2$mahal_dist)]
})
ns_m <- sapply(plot_mah, length)
boxplot(plot_mah, col=cols3, outline=FALSE, las=2,
        ylab="Mahalanobis distance",
        main="Mahalanobis dist by lead-time bin\n(train)",
        xaxt="n")
axis(1, at=seq_along(groups_all),
     labels=paste0(groups_all,"\nn=",ns_m),
     las=2, cex.axis=0.72, tick=FALSE)
abline(h=median(unlist(plot_mah["stable"]), na.rm=TRUE),
       col="grey50", lty=2)

# (c) Cohen's d profile
ok3 <- !is.na(res_perm$cohens_d)
barplot(res_perm$cohens_d[ok3],
        names.arg=BIN_LABELS[ok3],
        col=ifelse(res_perm$cohens_d[ok3]>0,"#E64B35","#4DBBD5"),
        las=2, ylab="Cohen's d",
        main="Effect size by bin\n(pooled Cohen's d vs stable)")
abline(h=0, lty=2)
abline(h=c(-0.2,0.2), lty=3, col="grey60")

# (d) sample counts
cnt_pre    <- sapply(BIN_LABELS, function(bn)
  sum(tr & !is.na(s2$leadtime_bin) & s2$leadtime_bin==bn))
cnt_stable <- sum(tr & s2$is_stable)
barplot(c(cnt_pre, cnt_stable),
        names.arg=c(BIN_LABELS,"stable"),
        col=c(rep("#E64B35",5),"#AAAAAA"),
        las=2, ylab="n samples (train)",
        main="Sample counts\n(train)",
        cex.names=0.75)

# (e) ABS predicted prob boxplots by bin
plot_abs <- lapply(groups_all, function(g) {
  if (g=="stable")
    s2$abs_pred_prob[tr & s2$is_stable & !is.na(s2$abs_pred_prob)]
  else
    s2$abs_pred_prob[tr & !is.na(s2$leadtime_bin) &
                       s2$leadtime_bin==g & !is.na(s2$abs_pred_prob)]
})
ns_a <- sapply(plot_abs, length)
boxplot(plot_abs, col=c(rep("#4DBBD5",5),"#AAAAAA"),
        outline=FALSE, las=2,
        ylab="ABS predicted probability",
        main="ABS classifier output by bin\n(train)",
        xaxt="n")
axis(1, at=seq_along(groups_all),
     labels=paste0(groups_all,"\nn=",ns_a),
     las=2, cex.axis=0.72, tick=FALSE)
abline(h=0.5, lty=2, col="grey50")

# (f) Q1 vs Q3 signal comparison
q1_p   <- q1r$p_medZ[match(BIN_LABELS, q1r$bin)]
q3_p   <- res_perm$p_perm
plot(seq_along(BIN_LABELS), q1_p,
     type="b", pch=16, col="#E64B35",
     ylim=c(0,1), xaxt="n",
     xlab="", ylab="permutation p-value",
     main="Q1 vs Q3 p-value profile\n(consistency check)")
lines(seq_along(BIN_LABELS), q3_p,
      type="b", pch=17, col="#00A087")
axis(1, at=seq_along(BIN_LABELS),
     labels=BIN_LABELS, las=2, cex.axis=0.75)
abline(h=0.05, lty=2, col="grey50")
legend("topright",
       legend=c("Q1 (all pre-onset)","Q3 (same, extended)"),
       col=c("#E64B35","#00A087"),
       pch=c(16,17), lty=1, bty="n", cex=0.8)

dev.off()
cat("saved outputs/Q3_supplementary.png\n")


# ============================================================
# 9. SAVE
# ============================================================

write.csv(res_q3,    "outputs/Q3_results.csv",    row.names=FALSE)
write.csv(abs_by_bin,"outputs/Q3_abs_by_bin.csv", row.names=FALSE)

sink("outputs/Q3_methods_numbers.txt")
cat("Q3 METHODS NUMBERS\n")
cat("generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("QUESTION\n")
cat("  Does predictive signal persist past the 45-75 d peak,\n")
cat("  or is it confined to that window?\n\n")

cat("PRIMARY METRIC: median |Z-score| (permutation test)\n")
for (i in seq_len(nrow(res_perm))) {
  cat(sprintf("  %-12s n=%2d  d=%s  p=%s  p_BH=%s\n",
              res_perm$bin[i], res_perm$n_bin[i],
              ifelse(is.na(res_perm$cohens_d[i]),"—",
                     sprintf("%+.3f",res_perm$cohens_d[i])),
              ifelse(is.na(res_perm$p_perm[i]),"—",
                     sprintf("%.4f",res_perm$p_perm[i])),
              ifelse(is.na(res_perm$p_bh[i]),"—",
                     sprintf("%.4f",res_perm$p_bh[i]))))
}

cat("\nTREND TEST\n")
cat("  Spearman r:", round(sp_cor,3), "\n")
cat("  p:         ", round(sp_p,4), "\n")

cat("\nSECONDARY METRIC: ABS predicted probability\n")
print(abs_by_bin, row.names=FALSE)

cat("\nINTERPRETATION\n")
sig_bins <- res_perm$bin[!is.na(res_perm$p_perm) & res_perm$p_perm < 0.05]
if (length(sig_bins)==0) {
  cat("  No bins significant in Q3 permutation test.\n")
  cat("  Signal from Q1 (45-75 d) does not replicate in extended analysis.\n")
  cat("  Possible reason: different sample composition after MAD clipping.\n")
} else if (length(sig_bins)==1 && "45–75 d" %in% sig_bins) {
  cat("  Signal CONFINED to 45-75 d window.\n")
  cat("  Narrow prodromal window: ~6-10 weeks before onset.\n")
} else {
  cat("  Signal significant in:", paste(sig_bins, collapse=", "), "\n")
  cat("  Extended early warning beyond 45-75 d.\n")
}
sink()

cat("wrote outputs/Q3_methods_numbers.txt\n")
cat("\n=== DONE — Q3 complete ===\n")
cat("Share Q3_results.csv and console output for interpretation.\n")


# null-safe infix for Panel A label positioning
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b