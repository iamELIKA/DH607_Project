## =============================================================================
## Script:  03_MTX_differential_analysis.R
## Project: IBD Microbial Transcript Activity & Candidate Proteins
## Purpose: Load the HUMAnN2 MTX gene-family table, filter features,
##          run linear mixed-effects models (lme4/lmerTest) for
##          CD vs nonIBD (primary) and UC vs nonIBD (secondary),
##          produce volcano plots and QQ plots.
##
## Pre-registered decisions (set before examining results):
##   Primary:   CD vs nonIBD
##   Secondary: UC vs nonIBD
##   Method:    lmerTest LMM, TSS + LOG, random = (1|Participant_ID)
##   Covariates: diagnosis + consent_age + sex + site_name + antibiotics_yn
##   Prevalence filter: >= 10% of analysis-set samples
##   Abundance filter:  mean RPK >= 1 (pre-TSS)
##   FDR correction: Benjamini-Hochberg
##
## Requires:
##   data/raw/MTX_genefamilies.tsv.gz     (from script 01)
##   data/processed/analysis_cohort_metadata.csv  (from script 02)
##
## Outputs:
##   data/processed/MTX_unstratified_filtered.rds
##   results/differential/CD_vs_nonIBD_results.csv
##   results/differential/UC_vs_nonIBD_results.csv
##   figures/volcano_CD_vs_nonIBD.png
##   figures/volcano_UC_vs_nonIBD.png
##   figures/qq_CD_vs_nonIBD.png
##   figures/qq_UC_vs_nonIBD.png
## =============================================================================

set.seed(42)

cat("=============================================================\n")
cat("Script 03: MTX Differential Analysis (lme4/lmerTest)\n")
cat("Session started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat("=============================================================\n\n")

## ---- 0. Packages ------------------------------------------------------------
required_pkgs <- c("data.table", "dplyr", "tidyr", "ggplot2", "stringr",
                   "lme4", "lmerTest", "ggrepel", "scales", "parallel")

missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(ggrepel)
  library(scales)
  library(parallel)
})

## ---- 1. Paths ---------------------------------------------------------------
base_dir  <- file.path(path.expand("~"), "Desktop", "607", "IBD_MTX_project")
gf_file   <- file.path(base_dir, "data", "raw",       "MTX_genefamilies.tsv.gz")
meta_file <- file.path(base_dir, "data", "processed",  "analysis_cohort_metadata.csv")
proc_dir  <- file.path(base_dir, "data", "processed")
diff_dir  <- file.path(base_dir, "transcriptomics", "microbial", "results", "differential")
figs_dir  <- file.path(base_dir, "transcriptomics", "microbial", "figures")
logs_dir  <- file.path(base_dir, "transcriptomics", "microbial", "logs")

for (d in c(proc_dir, diff_dir, figs_dir, logs_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ggplot2 theme (no cairo_pdf dependency)
theme_ibd <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#2C3E50"),
        strip.text    = element_text(colour = "white", face = "bold"),
        panel.grid.minor = element_blank())
diag_palette <- c("CD" = "#E74C3C", "UC" = "#3498DB", "nonIBD" = "#2ECC71")

## ---- 2. Load and validate cohort metadata -----------------------------------
cat("--- Loading cohort metadata ---\n")
if (!file.exists(meta_file)) stop("Run 02_metadata_and_cohort_QC.R first.")

meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("Cohort metadata:", nrow(meta), "samples,", ncol(meta), "columns\n")

# Make column names unique to avoid dplyr errors on duplicated survey columns
colnames(meta) <- make.unique(colnames(meta))

# Rename columns with spaces/special chars to be R-safe
# "Participant ID" → "Participant_ID"
# "External ID"   → "External_ID"
meta <- meta %>%
  rename(Participant_ID = `Participant ID`,
         External_ID    = `External ID`)

cat("By diagnosis:\n")
print(table(meta$diagnosis_std))

## ---- 3. Prepare model covariates --------------------------------------------
# Convert types
meta$consent_age     <- suppressWarnings(as.numeric(meta$consent_age))
meta$sex             <- factor(meta$sex)        # Female / Male
meta$site_name       <- factor(meta$site_name)
meta$antibiotics_yn  <- as.integer(meta$antibiotics_yn)  # 0/1
meta$diagnosis_std   <- factor(meta$diagnosis_std,
                                levels = c("nonIBD", "CD", "UC"))  # nonIBD = reference

# Remove samples with missing age (12 samples, 1.4%) — conservative approach
n_before <- nrow(meta)
meta <- meta %>% filter(!is.na(consent_age))
cat(sprintf("\nExcluded %d samples with missing age. Remaining: %d\n",
            n_before - nrow(meta), nrow(meta)))
cat("By diagnosis after age filter:\n")
print(table(meta$diagnosis_std))

# Record the sample IDs we will use
cohort_sample_ids <- meta$External_ID
cat("\nFinal cohort sample IDs (first 5):\n")
print(head(cohort_sample_ids, 5))

## ---- 4. Load gene-family table — structural safety check --------------------
cat("\n=============================================================\n")
cat("LOADING GENE-FAMILY TABLE\n")
cat("=============================================================\n")

if (!file.exists(gf_file)) stop("Gene-family file not found. Run script 01 first.")
cat("File size:", round(file.size(gf_file) / 1e6, 1), "MB (compressed)\n\n")

# SAFETY CHECK: Read only the header to confirm structure before full load
cat("--- Reading table header ---\n")
con <- gzfile(gf_file, "rt")
header_line <- readLines(con, n = 1)
close(con)

header_cols <- strsplit(header_line, "\t")[[1]]
feat_col_name <- header_cols[1]
raw_sample_cols <- header_cols[-1]

cat("Feature column name: '", feat_col_name, "'\n", sep = "")
cat("Number of sample columns:", length(raw_sample_cols), "\n")
cat("First 5 sample column names:\n")
print(head(raw_sample_cols, 5))

# Strip HUMAnN2 abundance suffixes from column names
strip_suffix <- function(ids) {
  ids <- sub("_Abundance-RPKs$",  "", ids)
  ids <- sub("_Abundance$",       "", ids)
  ids <- sub("_GeneFamilies$",    "", ids)
  ids <- sub("_Abundance-RELAB$", "", ids)
  ids
}
clean_sample_cols <- strip_suffix(raw_sample_cols)

# Check overlap with cohort
overlap_ids <- intersect(cohort_sample_ids, clean_sample_cols)
cat(sprintf("\nCohort samples: %d\n", length(cohort_sample_ids)))
cat(sprintf("Table columns (after suffix strip): %d\n", length(clean_sample_cols)))
cat(sprintf("Overlap: %d samples\n", length(overlap_ids)))

missing_from_table <- setdiff(cohort_sample_ids, clean_sample_cols)
if (length(missing_from_table) > 0) {
  cat(sprintf("WARNING: %d cohort samples NOT found in gene-family table:\n",
              length(missing_from_table)))
  print(head(missing_from_table, 10))
}

# Check for _TR technical replicate columns
tr_cols <- raw_sample_cols[grepl("_TR$", raw_sample_cols)]
cat(sprintf("\nTechnical replicate columns (_TR): %d — will be EXCLUDED\n",
            length(tr_cols)))

# Assess file size expectation
if (file.size(gf_file) / 1e6 < 50) {
  cat("\n*** WARNING: File is very small (", round(file.size(gf_file)/1e6, 1), " MB) ***\n")
  cat("A full HUMAnN2 merged genefamilies table for >800 samples should be\n")
  cat("much larger. This may indicate a truncated or per-sample stub file.\n")
  cat("Proceeding with load — will report actual dimensions.\n\n")
}

## ---- 5. Read the full gene-family table -------------------------------------
cat("--- Reading full gene-family table (may take several minutes) ---\n")
t0 <- proc.time()

gf_raw <- tryCatch({
  fread(file = gf_file, sep = "\t", header = TRUE,
        quote = "", fill = TRUE, data.table = TRUE,
        colClasses = "character")   # read all as character first for safety
}, error = function(e) {
  message("fread failed, trying base R: ", conditionMessage(e))
  con <- gzfile(gf_file, "rt")
  on.exit(close(con))
  as.data.table(
    read.table(con, sep = "\t", header = TRUE, quote = "",
               stringsAsFactors = FALSE, comment.char = "")
  )
})
t1 <- proc.time()
cat(sprintf("Load time: %.1f seconds\n", (t1 - t0)[["elapsed"]]))
cat("Raw table dimensions:", nrow(gf_raw), "rows ×", ncol(gf_raw), "columns\n\n")

## ---- 6. Rename columns and align to cohort ----------------------------------
# Rename feature column to a safe name
setnames(gf_raw, old = colnames(gf_raw)[1], new = "gene_family")

# Rename sample columns: strip suffixes
new_colnames <- c("gene_family", strip_suffix(colnames(gf_raw)[-1]))
setnames(gf_raw, new_colnames)

# Drop any technical replicate columns (_TR) if present
tr_to_drop <- colnames(gf_raw)[grepl("_TR$", colnames(gf_raw))]
if (length(tr_to_drop) > 0) {
  gf_raw[, (tr_to_drop) := NULL]
  cat("Dropped", length(tr_to_drop), "technical replicate columns.\n")
}

# Ensure keep_cols does not include _TR or columns that were dropped
keep_sample_ids <- setdiff(overlap_ids, tr_to_drop)
keep_sample_ids <- intersect(keep_sample_ids, colnames(gf_raw))

keep_cols <- c("gene_family", keep_sample_ids)
gf_cohort <- gf_raw[, ..keep_cols]
cat("After restricting to cohort samples:", nrow(gf_cohort), "features ×",
    ncol(gf_cohort) - 1, "samples\n")

## ---- 7. Convert to numeric and check value range ----------------------------
cat("\n--- Converting to numeric and checking value scale ---\n")
sample_cols_keep <- setdiff(colnames(gf_cohort), "gene_family")

# Convert sample columns to numeric
for (col in sample_cols_keep) {
  set(gf_cohort, j = col, value = suppressWarnings(as.numeric(gf_cohort[[col]])))
}

# Check first non-special row values
feat_vals_sample <- as.matrix(gf_cohort[!gene_family %in% c("UNMAPPED","UNINTEGRATED"),
                                         sample_cols_keep, with = FALSE][1:min(100, .N)])
cat("Value range in first 100 non-special features:\n")
cat(sprintf("  Min:    %.4f\n", min(feat_vals_sample, na.rm = TRUE)))
cat(sprintf("  Median: %.4f\n", median(feat_vals_sample, na.rm = TRUE)))
cat(sprintf("  Max:    %.2f\n", max(feat_vals_sample, na.rm = TRUE)))

if (max(feat_vals_sample, na.rm = TRUE) > 10) {
  cat("  → Values appear to be RPKs (large numbers). Correct.\n")
  cat("     TSS normalisation will be applied before log-transform.\n")
} else if (max(feat_vals_sample, na.rm = TRUE) <= 1.001) {
  cat("  → Values appear to be relative abundances (0–1).\n")
  cat("     Adjust normalization to 'NONE' and keep LOG transform.\n")
  cat("     *** STOP and review before proceeding ***\n")
}

## ---- 8. Separate feature types ----------------------------------------------
cat("\n--- Classifying HUMAnN2 feature types ---\n")

feat_ids <- gf_cohort$gene_family

# Taxon-stratified rows contain a pipe | after the gene family ID
# e.g. "UniRef90_A0A000|Bacteroides_fragilis"
# Unstratified have no pipe (or pipe only in the name/description part)
# UNMAPPED, UNINTEGRATED are special rows
is_unmapped     <- feat_ids == "UNMAPPED"
is_unintegrated <- feat_ids == "UNINTEGRATED"
# Stratified: has a "|" followed by a taxon name (not part of the ID description)
# UniRef90 IDs format: "UniRef90_XXXX: description" (colon+space before description)
# Stratified format:   "UniRef90_XXXX: description|Genus_species"
# Simple heuristic: stratified rows have "|" after the gene family base ID
is_stratified <- grepl("\\|", feat_ids) & !is_unmapped & !is_unintegrated
is_unstratified <- !is_stratified & !is_unmapped & !is_unintegrated

cat(sprintf("  Unstratified (community total):  %d\n", sum(is_unstratified)))
cat(sprintf("  Taxon-stratified:                %d\n", sum(is_stratified)))
cat(sprintf("  UNMAPPED:                        %d\n", sum(is_unmapped)))
cat(sprintf("  UNINTEGRATED:                    %d\n", sum(is_unintegrated)))

# Primary analysis uses UNSTRATIFIED features only
# (avoids double-counting; taxon-level analysis is secondary)
gf_unstrat <- gf_cohort[is_unstratified]
cat(sprintf("\nUsing %d unstratified features for primary analysis.\n", nrow(gf_unstrat)))

cat("\nSample unstratified feature IDs (first 5):\n")
print(head(gf_unstrat$gene_family, 5))

## ---- 9. Prevalence and abundance filtering ----------------------------------
cat("\n--- Applying prevalence and abundance filters ---\n")
n_samples_total <- length(sample_cols_keep)
n_features_before <- nrow(gf_unstrat)

# Prevalence: fraction of samples with value > 0
prev_vec <- rowMeans(as.matrix(gf_unstrat[, sample_cols_keep, with = FALSE]) > 0,
                     na.rm = TRUE)

# Abundance: mean across all samples
abund_vec <- rowMeans(as.matrix(gf_unstrat[, sample_cols_keep, with = FALSE]),
                      na.rm = TRUE)

# Auto-detect scale: RPK (max > 1) vs Relative Abundance (max <= 1)
is_relab <- max(feat_vals_sample, na.rm = TRUE) <= 1.001

# Filters (pre-registered: prevalence >= 10%, abundance filter adjusted for data scale)
PREV_THRESH  <- 0.10   # present in ≥ 10% of samples
ABUND_THRESH <- if (is_relab) 1e-7 else 1.0  # 1e-7 (0.00001%) for relative abundance, 1.0 for RPK

keep_feat <- prev_vec >= PREV_THRESH & abund_vec >= ABUND_THRESH
gf_filt   <- gf_unstrat[keep_feat]

cat(sprintf("  Scale detected:         %s\n", if (is_relab) "Relative Abundance (0-1)" else "RPK (raw counts)"))
cat(sprintf("  Prevalence threshold:  %.0f%% (≥ %.0f samples)\n",
            100 * PREV_THRESH, ceiling(PREV_THRESH * n_samples_total)))
cat(sprintf("  Abundance threshold:   mean ≥ %s\n", format(ABUND_THRESH, scientific = TRUE)))
cat(sprintf("  Features before filter: %d\n", n_features_before))
cat(sprintf("  Features after filter:  %d\n", nrow(gf_filt)))
cat(sprintf("  Removed:                %d (%.1f%%)\n",
            n_features_before - nrow(gf_filt),
            100 * (1 - nrow(gf_filt) / n_features_before)))

# Sparsity of filtered set
sparsity_filt <- mean(
  as.matrix(gf_filt[, sample_cols_keep, with = FALSE]) == 0, na.rm = TRUE)
cat(sprintf("  Sparsity after filter:  %.1f%% zeros\n", 100 * sparsity_filt))

## ---- 10. Build input matrices -----------------------------------------------
cat("\n--- Building analysis input matrices ---\n")

# Transpose: samples as ROWS, features as COLUMNS
feat_mat <- as.matrix(gf_filt[, sample_cols_keep, with = FALSE])
rownames(feat_mat) <- gf_filt$gene_family
feat_mat_t <- t(feat_mat)  # now: samples × features
cat("Feature matrix (transposed): ", nrow(feat_mat_t), "samples ×",
    ncol(feat_mat_t), "features\n")

# Sanitise feature names: replace ":" and space with "_"
colnames(feat_mat_t) <- gsub("[: ]", "_", colnames(feat_mat_t))
# Store mapping: sanitised name → original ID
feat_name_map <- data.frame(
  sanitised = colnames(feat_mat_t),
  original  = gf_filt$gene_family,
  stringsAsFactors = FALSE
)
write.csv(feat_name_map,
          file.path(proc_dir, "feature_name_mapping.csv"),
          row.names = FALSE)

# Align metadata to the samples in feat_mat_t
meta_aligned <- meta %>%
  filter(External_ID %in% rownames(feat_mat_t)) %>%
  arrange(match(External_ID, rownames(feat_mat_t))) %>%
  as.data.frame()
rownames(meta_aligned) <- meta_aligned$External_ID

# Verify row order matches
stopifnot(all(rownames(feat_mat_t) == rownames(meta_aligned)))
cat("Metadata aligned. Rows match: YES\n")

## ---- 11. Save filtered feature matrix for downstream scripts ----------------
saveRDS(list(
  feat_mat_t        = feat_mat_t,
  feat_name_map     = feat_name_map,
  gf_filt           = gf_filt,
  sample_ids        = rownames(feat_mat_t),
  PREV_THRESH       = PREV_THRESH,
  ABUND_THRESH      = ABUND_THRESH,
  n_features_before = n_features_before,
  n_features_after  = nrow(gf_filt)
), file.path(proc_dir, "MTX_unstratified_filtered.rds"))
cat("Saved filtered matrix to: data/processed/MTX_unstratified_filtered.rds\n")

## ---- 12. Function: run LMM analysis with lme4/lmerTest ----------------------
##
## Steps per feature:
##   1) TSS normalise (divide by sample total across all features)
##   2) log2(x + pseudocount) transform
##   3) Fit: lmer(abundance ~ diagnosis_std + consent_age + sex + site_name +
##                antibiotics_yn + (1 | Participant_ID))
##   4) Extract Satterthwaite p-values for the diagnosis term
##   5) BH FDR correction across all features
##
## Uses parallel::mclapply with mc.cores for speed.
## Falls back to lm() (no random effect) if lmer fails to converge.
## --------------------------------------------------------------------------

run_lmm_analysis <- function(feat_mat, meta_df, comparison_label,
                              ref_group = "nonIBD", cores = 32) {
  out_dir <- file.path(diff_dir, comparison_label)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  cat("\n--- Running LMM analysis:", comparison_label, "---\n")
  cat("  N samples:", nrow(feat_mat), "\n")
  cat("  N features:", ncol(feat_mat), "\n")
  cat("  Reference group:", ref_group, "\n")
  cat("  Cores:", cores, "\n\n")

  # ---- 12a. TSS normalisation (per-sample) ----------------------------------
  sample_totals <- rowSums(feat_mat, na.rm = TRUE)
  # Guard against zero-sum samples
  sample_totals[sample_totals == 0] <- 1
  feat_tss <- feat_mat / sample_totals

  # ---- 12b. Log2 transform --------------------------------------------------
  # Pseudocount = half the minimum non-zero value across the entire TSS matrix
  nonzero_vals <- feat_tss[feat_tss > 0]
  if (length(nonzero_vals) == 0) {
    warning("All values are zero in feature matrix. Aborting: ", comparison_label)
    return(NULL)
  }
  pseudocount <- min(nonzero_vals) / 2
  cat(sprintf("  Pseudocount (min non-zero / 2): %.2e\n", pseudocount))
  feat_log <- log2(feat_tss + pseudocount)

  # ---- 12c. Build per-feature data and fit LMMs in parallel -----------------
  meta_df <- as.data.frame(meta_df)
  # Ensure diagnosis reference level
  meta_df$diagnosis_std <- relevel(factor(meta_df$diagnosis_std), ref = ref_group)
  meta_df$site_name     <- relevel(factor(meta_df$site_name), ref = "Cedars-Sinai")
  meta_df$sex           <- factor(meta_df$sex)
  meta_df$antibiotics_yn <- as.numeric(meta_df$antibiotics_yn)
  meta_df$consent_age   <- as.numeric(meta_df$consent_age)

  feature_names <- colnames(feat_log)

  fit_one_feature <- function(feat_idx) {
    fn <- feature_names[feat_idx]
    df <- meta_df
    df$y <- feat_log[, feat_idx]

    # Try LMM with random intercept for participant
    result <- tryCatch({
      mod <- lmerTest::lmer(
        y ~ diagnosis_std + consent_age + sex + site_name + antibiotics_yn +
          (1 | Participant_ID),
        data = df,
        REML = TRUE,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl   = list(maxfun = 1e5)
        )
      )
      # Satterthwaite p-values
      coef_tbl <- as.data.frame(lmerTest::summary(mod)$coefficients)
      # Extract rows matching "diagnosis_std" (e.g., diagnosis_stdCD)
      diag_rows <- grep("^diagnosis_std", rownames(coef_tbl))
      if (length(diag_rows) == 0) {
        return(data.frame(feature = fn, coef = NA, stderr = NA, t_val = NA,
                          df_satt = NA, pval = NA, method = "lmer",
                          converged = TRUE, stringsAsFactors = FALSE))
      }
      # Take the diagnosis term row
      dr <- coef_tbl[diag_rows[1], , drop = FALSE]
      data.frame(
        feature   = fn,
        coef      = dr[1, "Estimate"],
        stderr    = dr[1, "Std. Error"],
        t_val     = dr[1, "t value"],
        df_satt   = dr[1, "df"],
        pval      = dr[1, "Pr(>|t|)"],
        method    = "lmer",
        converged = TRUE,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      # Fallback: simple lm() without random effect
      tryCatch({
        mod_lm <- lm(
          y ~ diagnosis_std + consent_age + sex + site_name + antibiotics_yn,
          data = df
        )
        coef_lm <- summary(mod_lm)$coefficients
        diag_rows <- grep("^diagnosis_std", rownames(coef_lm))
        if (length(diag_rows) == 0) {
          return(data.frame(feature = fn, coef = NA, stderr = NA, t_val = NA,
                            df_satt = NA, pval = NA, method = "lm_fallback",
                            converged = FALSE, stringsAsFactors = FALSE))
        }
        dr <- coef_lm[diag_rows[1], , drop = FALSE]
        data.frame(
          feature   = fn,
          coef      = dr[1, "Estimate"],
          stderr    = dr[1, "Std. Error"],
          t_val     = dr[1, "t value"],
          df_satt   = NA_real_,
          pval      = dr[1, "Pr(>|t|)"],
          method    = "lm_fallback",
          converged = FALSE,
          stringsAsFactors = FALSE
        )
      }, error = function(e2) {
        data.frame(feature = fn, coef = NA, stderr = NA, t_val = NA,
                   df_satt = NA, pval = NA, method = "failed",
                   converged = FALSE, stringsAsFactors = FALSE)
      })
    })
    result
  }

  # Run in parallel
  cat("  Fitting", length(feature_names), "features across", cores, "cores ...\n")
  t0 <- proc.time()

  results_list <- parallel::mclapply(
    seq_along(feature_names),
    fit_one_feature,
    mc.cores = cores
  )

  t1 <- proc.time()
  cat(sprintf("  LMM fitting completed in %.1f seconds\n", (t1 - t0)[["elapsed"]]))

  # Combine into data.frame
  res <- do.call(rbind, results_list)

  # Report convergence
  n_lmer    <- sum(res$method == "lmer", na.rm = TRUE)
  n_fallback <- sum(res$method == "lm_fallback", na.rm = TRUE)
  n_failed  <- sum(res$method == "failed", na.rm = TRUE)
  cat(sprintf("  Converged (lmer):   %d\n", n_lmer))
  cat(sprintf("  Fallback (lm):      %d\n", n_fallback))
  cat(sprintf("  Failed:             %d\n", n_failed))

  # BH FDR correction (only on non-NA p-values)
  res$qval <- NA_real_
  valid_p <- !is.na(res$pval)
  res$qval[valid_p] <- p.adjust(res$pval[valid_p], method = "BH")

  # Add comparison label
  res$comparison <- comparison_label

  # Sort by p-value
  res <- res[order(res$pval), ]

  cat(sprintf("  FDR < 0.05:  %d\n", sum(res$qval < 0.05, na.rm = TRUE)))
  cat(sprintf("  FDR < 0.20:  %d\n", sum(res$qval < 0.20, na.rm = TRUE)))

  # Save LMM log
  writeLines(
    c(paste("Comparison:", comparison_label),
      paste("N samples:", nrow(feat_mat)),
      paste("N features:", ncol(feat_mat)),
      paste("Pseudocount:", pseudocount),
      paste("Converged (lmer):", n_lmer),
      paste("Fallback (lm):", n_fallback),
      paste("Failed:", n_failed),
      paste("FDR<0.05:", sum(res$qval < 0.05, na.rm = TRUE)),
      paste("FDR<0.20:", sum(res$qval < 0.20, na.rm = TRUE)),
      paste("Completed:", Sys.time())),
    file.path(out_dir, "lmm_analysis_log.txt")
  )

  res
}

## ---- 13. Primary analysis: CD vs nonIBD -------------------------------------
# Subset to CD + nonIBD only
cd_ctrl_ids <- meta_aligned %>%
  filter(diagnosis_std %in% c("CD", "nonIBD")) %>%
  pull(External_ID)

feat_cd_ctrl    <- feat_mat_t[cd_ctrl_ids, ]
meta_cd_ctrl    <- meta_aligned[cd_ctrl_ids, ]

cat(sprintf("\nCD vs nonIBD subset: %d samples (%d CD, %d nonIBD)\n",
            nrow(feat_cd_ctrl),
            sum(meta_cd_ctrl$diagnosis_std == "CD"),
            sum(meta_cd_ctrl$diagnosis_std == "nonIBD")))

res_cd <- run_lmm_analysis(
  feat_mat         = feat_cd_ctrl,
  meta_df          = meta_cd_ctrl,
  comparison_label = "CD_vs_nonIBD",
  ref_group        = "nonIBD",
  cores            = 32
)

## ---- 14. Secondary analysis: UC vs nonIBD -----------------------------------
uc_ctrl_ids <- meta_aligned %>%
  filter(diagnosis_std %in% c("UC", "nonIBD")) %>%
  pull(External_ID)

feat_uc_ctrl <- feat_mat_t[uc_ctrl_ids, ]
meta_uc_ctrl <- meta_aligned[uc_ctrl_ids, ]

cat(sprintf("\nUC vs nonIBD subset: %d samples (%d UC, %d nonIBD)\n",
            nrow(feat_uc_ctrl),
            sum(meta_uc_ctrl$diagnosis_std == "UC"),
            sum(meta_uc_ctrl$diagnosis_std == "nonIBD")))

res_uc <- run_lmm_analysis(
  feat_mat         = feat_uc_ctrl,
  meta_df          = meta_uc_ctrl,
  comparison_label = "UC_vs_nonIBD",
  ref_group        = "nonIBD",
  cores            = 32
)

## ---- 15. Merge results and add original feature names -----------------------
add_original_names <- function(res_df, name_map) {
  if (is.null(res_df)) return(NULL)
  res_df %>%
    left_join(name_map, by = c("feature" = "sanitised")) %>%
    rename(gene_family_id = original) %>%
    # Parse taxon stratum (NA for unstratified)
    mutate(
      taxonomic_stratum = if_else(
        grepl("\\|", gene_family_id),
        str_extract(gene_family_id, "(?<=\\|).*$"),
        "community_total"
      ),
      base_gene_family = str_remove(gene_family_id, "\\|.*$")
    ) %>%
    arrange(qval, pval)
}

res_cd_clean <- add_original_names(res_cd, feat_name_map)
res_uc_clean <- add_original_names(res_uc, feat_name_map)

## Save results tables
if (!is.null(res_cd_clean)) {
  write.csv(res_cd_clean,
            file.path(diff_dir, "CD_vs_nonIBD_results.csv"),
            row.names = FALSE)
  cat("\nTop CD vs nonIBD results (FDR < 0.25):\n")
  top_cd <- res_cd_clean %>% filter(qval < 0.25) %>%
    select(gene_family_id, coef, stderr, pval, qval, taxonomic_stratum) %>%
    head(20)
  print(as.data.frame(top_cd))
}

if (!is.null(res_uc_clean)) {
  write.csv(res_uc_clean,
            file.path(diff_dir, "UC_vs_nonIBD_results.csv"),
            row.names = FALSE)
  cat("\nTop UC vs nonIBD results (FDR < 0.25):\n")
  top_uc <- res_uc_clean %>% filter(qval < 0.25) %>%
    select(gene_family_id, coef, stderr, pval, qval, taxonomic_stratum) %>%
    head(20)
  print(as.data.frame(top_uc))
}

## ---- 16. Sensitivity analysis: prevalence threshold 20% --------------------
cat("\n=== Sensitivity analysis: prevalence threshold 20% ===\n")
keep_feat_20 <- prev_vec >= 0.20 & abund_vec >= ABUND_THRESH
gf_filt_20   <- gf_unstrat[keep_feat_20]
cat(sprintf("Features at ≥20%% prevalence: %d (vs %d at ≥10%%)\n",
            nrow(gf_filt_20), nrow(gf_filt)))

# Only run if substantially different from primary
if (nrow(gf_filt_20) < nrow(gf_filt) * 0.95) {
  feat_20      <- t(as.matrix(gf_filt_20[, sample_cols_keep, with = FALSE]))
  rownames(feat_20) <- rownames(feat_mat_t)
  colnames(feat_20) <- gsub("[: ]", "_", gf_filt_20$gene_family)

  res_cd_20 <- run_lmm_analysis(
    feat_mat         = feat_20[cd_ctrl_ids, ],
    meta_df          = meta_cd_ctrl,
    comparison_label = "CD_vs_nonIBD_prev20pct",
    ref_group        = "nonIBD",
    cores            = 32
  )

  # Compare: how many FDR<0.05 features overlap?
  if (!is.null(res_cd) && !is.null(res_cd_20)) {
    sig_10 <- res_cd$feature[res_cd$qval < 0.05]
    sig_20 <- res_cd_20$feature[res_cd_20$qval < 0.05]
    cat(sprintf("  FDR<0.05 in ≥10%% filter:  %d\n", length(sig_10)))
    cat(sprintf("  FDR<0.05 in ≥20%% filter:  %d\n", length(sig_20)))
    cat(sprintf("  Overlap:                   %d\n", length(intersect(sig_10, sig_20))))
    write.csv(data.frame(
      feature = sig_10,
      in_20pct_sig = sig_10 %in% sig_20
    ), file.path(diff_dir, "sensitivity_prev_filter_comparison.csv"),
    row.names = FALSE)
  }
}

## ---- 17. Visualisations -----------------------------------------------------
make_volcano <- function(res_df, title_str, fdr_line = 0.05,
                         label_n_top = 15) {
  if (is.null(res_df) || nrow(res_df) == 0) {
    message("No results for: ", title_str); return(NULL)
  }

  df <- res_df %>%
    mutate(
      sig_label = case_when(
        qval < 0.05  ~ "FDR < 0.05",
        qval < 0.20  ~ "FDR < 0.20",
        TRUE         ~ "NS"
      ),
      sig_label = factor(sig_label, levels = c("FDR < 0.05","FDR < 0.20","NS")),
      neg_log10_p = -log10(pmax(pval, 1e-30)),
      label_this  = rank(pval) <= label_n_top
    )

  p <- ggplot(df, aes(x = coef, y = neg_log10_p, colour = sig_label)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(fdr_line), linetype = "dotted",
               colour = "grey30") +
    scale_colour_manual(
      values = c("FDR < 0.05" = "#C0392B", "FDR < 0.20" = "#E67E22", "NS" = "grey70"),
      name = NULL
    ) +
    geom_text_repel(
      data = df %>% filter(label_this),
      aes(label = str_trunc(base_gene_family, 30)),
      size = 2.5, max.overlaps = 20, colour = "grey10"
    ) +
    labs(title = title_str,
         x = "LMM coefficient (log2-scale effect size)",
         y = expression(-log[10](p-value))) +
    theme_ibd
  p
}

make_qq <- function(res_df, title_str) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  n  <- nrow(res_df)
  df <- res_df %>%
    arrange(pval) %>%
    mutate(
      expected = -log10((seq_len(n)) / (n + 1)),
      observed = -log10(pmax(pval, 1e-30))
    )
  lambda <- median(qchisq(1 - res_df$pval, df = 1), na.rm = TRUE) /
            qchisq(0.5, df = 1)

  ggplot(df, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    geom_point(alpha = 0.4, size = 1) +
    labs(title = paste(title_str, sprintf("(λ = %.2f)", lambda)),
         x = expression("Expected " * -log[10](p)),
         y = expression("Observed " * -log[10](p))) +
    theme_ibd
}

# Generate plots — show in console/viewer, then save as PNG images
if (!is.null(res_cd_clean)) {
  p_volc_cd <- make_volcano(res_cd_clean, "CD vs nonIBD — MTX gene families")
  p_qq_cd   <- make_qq(res_cd_clean, "QQ: CD vs nonIBD")

  if (!is.null(p_volc_cd)) {
    # 1. Show in console / plot window
    tryCatch(print(p_volc_cd), error = function(e) invisible(NULL))
    # 2. Save as PNG image
    volc_cd_path <- file.path(figs_dir, "volcano_CD_vs_nonIBD.png")
    ggsave(volc_cd_path, plot = p_volc_cd, width = 9, height = 6, dpi = 300)
    cat("Saved image:", volc_cd_path, "\n")
  }
  if (!is.null(p_qq_cd)) {
    # 1. Show in console / plot window
    tryCatch(print(p_qq_cd), error = function(e) invisible(NULL))
    # 2. Save as PNG image
    qq_cd_path <- file.path(figs_dir, "qq_CD_vs_nonIBD.png")
    ggsave(qq_cd_path, plot = p_qq_cd, width = 5, height = 5, dpi = 300)
    cat("Saved image:", qq_cd_path, "\n")
  }
}

if (!is.null(res_uc_clean)) {
  p_volc_uc <- make_volcano(res_uc_clean, "UC vs nonIBD — MTX gene families")
  p_qq_uc   <- make_qq(res_uc_clean, "QQ: UC vs nonIBD")

  if (!is.null(p_volc_uc)) {
    # 1. Show in console / plot window
    tryCatch(print(p_volc_uc), error = function(e) invisible(NULL))
    # 2. Save as PNG image
    volc_uc_path <- file.path(figs_dir, "volcano_UC_vs_nonIBD.png")
    ggsave(volc_uc_path, plot = p_volc_uc, width = 9, height = 6, dpi = 300)
    cat("Saved image:", volc_uc_path, "\n")
  }
  if (!is.null(p_qq_uc)) {
    # 1. Show in console / plot window
    tryCatch(print(p_qq_uc), error = function(e) invisible(NULL))
    # 2. Save as PNG image
    qq_uc_path <- file.path(figs_dir, "qq_UC_vs_nonIBD.png")
    ggsave(qq_uc_path, plot = p_qq_uc, width = 5, height = 5, dpi = 300)
    cat("Saved image:", qq_uc_path, "\n")
  }
}

# Clean up default Rplots.pdf if created by headless R batch mode during print()
if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}

## ---- 18. Summary table of significant features ------------------------------
cat("\n=== SUMMARY OF SIGNIFICANT ASSOCIATIONS ===\n")
if (!is.null(res_cd_clean)) {
  cat("CD vs nonIBD:\n")
  cat(sprintf("  Features tested:          %d\n", nrow(res_cd_clean)))
  cat(sprintf("  FDR < 0.05 (primary):     %d\n",
              sum(res_cd_clean$qval < 0.05, na.rm = TRUE)))
  cat(sprintf("  FDR < 0.20 (exploratory): %d\n",
              sum(res_cd_clean$qval < 0.20, na.rm = TRUE)))
  cat(sprintf("  Up in CD (coef > 0):      %d\n",
              sum(res_cd_clean$qval < 0.05 & res_cd_clean$coef > 0, na.rm = TRUE)))
  cat(sprintf("  Down in CD (coef < 0):    %d\n",
              sum(res_cd_clean$qval < 0.05 & res_cd_clean$coef < 0, na.rm = TRUE)))
}
if (!is.null(res_uc_clean)) {
  cat("UC vs nonIBD:\n")
  cat(sprintf("  Features tested:          %d\n", nrow(res_uc_clean)))
  cat(sprintf("  FDR < 0.05 (primary):     %d\n",
              sum(res_uc_clean$qval < 0.05, na.rm = TRUE)))
  cat(sprintf("  FDR < 0.20 (exploratory): %d\n",
              sum(res_uc_clean$qval < 0.20, na.rm = TRUE)))
}

## ---- 19. Save session info --------------------------------------------------
session_log <- file.path(logs_dir, paste0("03_session_", format(Sys.Date(), "%Y%m%d"), ".txt"))
sink(session_log); sessionInfo(); sink()

cat("\n=============================================================\n")
cat("SCRIPT 03 COMPLETE\n")
cat("=============================================================\n")
cat("Key outputs:\n")
cat("  data/processed/MTX_unstratified_filtered.rds  ← used by scripts 04 & 05\n")
cat("  results/differential/CD_vs_nonIBD_results.csv\n")
cat("  results/differential/UC_vs_nonIBD_results.csv\n")
cat("  figures/volcano_*.png, figures/qq_*.png\n")
cat("\nNEXT STEPS:\n")
cat("  1. Review the summary table printed above.\n")
cat("  2. Check QQ plots: lambda should be close to 1 for well-calibrated tests.\n")
cat("     Lambda >> 1 suggests confounding or model mis-specification.\n")
cat("  3. Inspect top CD vs nonIBD results for biological plausibility\n")
cat("     (butyrate pathways should be DOWN, oxidative stress pathways UP).\n")
cat("  4. Share the SUMMARY section output before running script 04 (MGX adjustment).\n")
cat("=============================================================\n")
cat("Session ended:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

sink("03.txt")

