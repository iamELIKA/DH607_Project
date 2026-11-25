## =============================================================================
## Script:  04_MGX_adjustment_and_transcriptional_activity.R
## Project: IBD Microbial Transcript Activity & Candidate Proteins
## Purpose: Load Enzyme Commission (EC) tables for both MTX (RNA) and MGX (DNA),
##          match same-visit MTX-MGX paired samples, and compute DNA-adjusted
##          transcriptional activity (RNA ~ DNA + diagnosis + covariates).
##          Distinguishes true transcriptional regulation from DNA abundance shifts.
##
## Inputs:
##   data/raw/MTX_ecs.tsv.gz
##   data/raw/MGX_ecs.tsv.gz
##   data/processed/MTX_MGX_same_visit_pairs.csv
##   data/processed/analysis_cohort_metadata.csv
##
## Outputs:
##   results/differential/CD_vs_nonIBD_MGX_adjusted_results.csv
##   results/differential/UC_vs_nonIBD_MGX_adjusted_results.csv
##   results/differential/transcriptional_activity_summary.csv
##   figures/rna_vs_dna_scatter_CD.png
##   figures/volcano_MGX_adjusted_CD.png
##   figures/volcano_MGX_adjusted_UC.png
##   logs/04_session_<date>.txt
## =============================================================================

set.seed(42)

cat("=============================================================\n")
cat("Script 04: MGX Adjustment & Microbial Transcriptional Activity (EC-level)\n")
cat("Session started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat("=============================================================\n\n")

## ---- 0. Packages ------------------------------------------------------------
required_pkgs <- c("data.table", "dplyr", "tidyr", "ggplot2", "scales", "ggrepel", "parallel")
missing_pkgs  <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(parallel)
})

## ---- 1. Paths & Theme -------------------------------------------------------
base_dir   <- file.path(path.expand("~"), "Desktop", "607", "IBD_MTX_project")
data_raw   <- file.path(base_dir, "data", "raw")
data_proc  <- file.path(base_dir, "data", "processed")
diff_dir   <- file.path(base_dir, "results", "differential")
figs_dir   <- file.path(base_dir, "figures")
logs_dir   <- file.path(base_dir, "logs")

for (d in c(data_raw, data_proc, diff_dir, figs_dir, logs_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

theme_ibd <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#2C3E50"),
        strip.text       = element_text(colour = "white", face = "bold"),
        panel.grid.minor = element_blank())

## ---- 2. Verify Input Files --------------------------------------------------
mtx_ec_file <- file.path(data_raw, "MTX_ecs.tsv.gz")
mgx_ec_file <- file.path(data_raw, "MGX_ecs.tsv.gz")
pairs_file  <- file.path(data_proc, "MTX_MGX_same_visit_pairs.csv")
meta_file   <- file.path(data_proc, "analysis_cohort_metadata.csv")

if (!file.exists(mtx_ec_file) || file.size(mtx_ec_file) < 100000) {
  stop("Missing or incomplete data/raw/MTX_ecs.tsv.gz. Download from ibdmdb.org/results and place in data/raw/.")
}
if (!file.exists(mgx_ec_file) || file.size(mgx_ec_file) < 100000) {
  stop("Missing or incomplete data/raw/MGX_ecs.tsv.gz. Download from ibdmdb.org/results and place in data/raw/.")
}
if (!file.exists(pairs_file)) stop("Pairs file missing. Run Script 02 first.")
if (!file.exists(meta_file))  stop("Metadata file missing. Run Script 02 first.")

cat(sprintf("Found MTX_ecs.tsv.gz: %.1f MB\n", file.size(mtx_ec_file) / 1e6))
cat(sprintf("Found MGX_ecs.tsv.gz: %.1f MB\n", file.size(mgx_ec_file) / 1e6))

## ---- 3. Load Metadata & Paired Samples --------------------------------------
pairs <- read.csv(pairs_file, stringsAsFactors = FALSE)
pairs_clean <- pairs %>% distinct(mtx_id, .keep_all = TRUE)
cat(sprintf("Clean unique same-visit MTX-MGX pairs: %d\n", nrow(pairs_clean)))

meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
colnames(meta) <- make.unique(colnames(meta))

meta_paired <- meta %>%
  filter(`External ID` %in% pairs_clean$mtx_id) %>%
  left_join(pairs_clean, by = c("External ID" = "mtx_id")) %>%
  filter(diagnosis_std %in% c("CD", "UC", "nonIBD"))

cat(sprintf("Paired cohort: %d samples (CD=%d, UC=%d, nonIBD=%d)\n",
            nrow(meta_paired),
            sum(meta_paired$diagnosis_std == "CD"),
            sum(meta_paired$diagnosis_std == "UC"),
            sum(meta_paired$diagnosis_std == "nonIBD")))

if (nrow(meta_paired) == 0) {
  stop("No paired samples found. Verify analysis_cohort_metadata.csv has valid diagnosis_std.")
}

## ---- 4. Helper to Read & Normalise EC Tables --------------------------------
strip_suffix <- function(ids) {
  ids <- sub("_Abundance-RPKs$",  "", ids)
  ids <- sub("_Abundance-RELAB$", "", ids)
  ids <- sub("_Abundance$",       "", ids)
  ids <- sub("_GeneFamilies$",    "", ids)
  ids <- sub("_ecs$",             "", ids)
  ids <- sub("_ECs$",             "", ids)
  ids
}

load_ec_matrix <- function(filepath, target_sample_ids, label = "EC") {
  cat(sprintf("\n--- Loading %s EC table (%s) ---\n", label, basename(filepath)))
  con <- gzfile(filepath, "rt")
  h_line <- readLines(con, n = 1)
  close(con)

  raw_cols   <- strsplit(h_line, "\t")[[1]]
  clean_cols <- c(raw_cols[1], strip_suffix(raw_cols[-1]))
  matched    <- intersect(target_sample_ids, clean_cols[-1])

  cat(sprintf("%s: %d / %d requested sample columns matched\n",
              label, length(matched), length(target_sample_ids)))

  if (length(matched) == 0) {
    stop(sprintf("No sample IDs matched in %s. Check column naming.", basename(filepath)))
  }

  select_idx <- c(1L, which(clean_cols %in% matched))
  dt <- fread(filepath, header = TRUE, sep = "\t", select = select_idx)
  setnames(dt, c("ec_id", strip_suffix(colnames(dt)[-1])))

  # Keep only unstratified ECs (no "|"), drop UNMAPPED / UNGROUPED
  dt <- dt[!grepl("|", ec_id, fixed = TRUE) &
             !grepl("^UNMAPPED|^UNGROUPED", ec_id)]

  # Deduplicate any duplicate sample columns if present
  samp_cols <- unique(intersect(matched, colnames(dt)))
  mat <- as.matrix(dt[, samp_cols, with = FALSE])
  rownames(mat) <- dt$ec_id
  storage.mode(mat) <- "numeric"

  cat(sprintf("%s unstratified matrix: %d ECs x %d samples\n",
              label, nrow(mat), ncol(mat)))
  mat
}

normalise_ec_matrix <- function(mat, label = "EC") {
  max_val <- max(mat[1:min(50, nrow(mat)), ], na.rm = TRUE)
  cat(sprintf("%s max value (first 50 ECs): %.6f\n", label, max_val))

  if (max_val <= 1.001) {
    cat(sprintf("%s appears to be relative abundance (max <= 1). Applying log2 directly.\n", label))
    pseudo <- if (any(mat > 0, na.rm = TRUE)) min(mat[mat > 0], na.rm = TRUE) / 2 else 1e-9
    log_mat <- log2(mat + pseudo)
  } else {
    cat(sprintf("%s appears to be RPK/CPM scale (max > 1). Applying TSS + log2.\n", label))
    col_tot <- colSums(mat, na.rm = TRUE)
    col_tot[col_tot == 0] <- 1
    tss_mat <- sweep(mat, 2, col_tot, "/")
    pseudo  <- if (any(tss_mat > 0, na.rm = TRUE)) min(tss_mat[tss_mat > 0], na.rm = TRUE) / 2 else 1e-9
    log_mat <- log2(tss_mat + pseudo)
  }
  log_mat
}

mtx_raw <- load_ec_matrix(mtx_ec_file, meta_paired$`External ID`, label = "MTX")
mgx_raw <- load_ec_matrix(mgx_ec_file, meta_paired$mgx_id,        label = "MGX")

## ---- 5. Filter & Align Paired Samples and ECs -------------------------------
# Filter MTX ECs: >= 10% prevalence and mean abundance threshold
mtx_max <- max(mtx_raw[1:min(50, nrow(mtx_raw)), ], na.rm = TRUE)
abund_thresh <- if (mtx_max <= 1.001) 1e-7 else 1.0

prev_pass  <- rowMeans(mtx_raw > 0, na.rm = TRUE) >= 0.10
abund_pass <- rowMeans(mtx_raw, na.rm = TRUE) >= abund_thresh
keep_ecs   <- intersect(rownames(mtx_raw)[prev_pass & abund_pass], rownames(mgx_raw))

cat(sprintf("\nECs passing MTX prevalence (>=10%%) & abundance (>=%.1e) & present in MGX: %d\n",
            abund_thresh, length(keep_ecs)))

if (length(keep_ecs) == 0) {
  stop("Zero ECs passed filtering and overlap. Check EC IDs and abundance scale.")
}

mtx_log <- normalise_ec_matrix(mtx_raw[keep_ecs, , drop = FALSE], label = "MTX")
mgx_log <- normalise_ec_matrix(mgx_raw[keep_ecs, , drop = FALSE], label = "MGX")

sample_map <- meta_paired %>%
  filter(`External ID` %in% colnames(mtx_log),
         mgx_id %in% colnames(mgx_log))

cat(sprintf("Final aligned paired samples: %d (CD=%d, UC=%d, nonIBD=%d)\n",
            nrow(sample_map),
            sum(sample_map$diagnosis_std == "CD"),
            sum(sample_map$diagnosis_std == "UC"),
            sum(sample_map$diagnosis_std == "nonIBD")))

## ---- 6. Fit Unadjusted & DNA-Adjusted Models --------------------------------
run_ec_comparison <- function(case_group, sample_df, mtx_mat, mgx_mat, n_cores) {
  cat(sprintf("\n--- Running %s vs nonIBD (unadjusted + DNA-adjusted) across %d cores ---\n",
              case_group, n_cores))

  sub_df <- sample_df %>% filter(diagnosis_std %in% c(case_group, "nonIBD"))
  ref_site <- if ("Cedars-Sinai" %in% sub_df$site_name) "Cedars-Sinai" else sort(unique(sub_df$site_name))[1]

  sub_df$diagnosis_std  <- relevel(factor(sub_df$diagnosis_std), ref = "nonIBD")
  sub_df$site_name      <- relevel(factor(sub_df$site_name), ref = ref_site)
  sub_df$sex            <- factor(sub_df$sex)
  sub_df$antibiotics_yn <- as.numeric(sub_df$antibiotics_yn)
  sub_df$consent_age    <- as.numeric(sub_df$consent_age)

  ecs <- rownames(mtx_mat)

  fit_one_ec <- function(idx) {
    ec <- ecs[idx]
    df <- sub_df
    df$rna <- mtx_mat[ec, df$`External ID`]
    df$dna <- mgx_mat[ec, df$mgx_id]
    df <- df[!is.na(df$rna) & !is.na(df$dna), ]

    if (nrow(df) < 10) {
      return(data.frame(
        ec_id = ec, n_samples = nrow(df),
        coef_unadj = NA, pval_unadj = NA,
        coef_dna = NA, pval_dna = NA,
        coef_diag_adj = NA, stderr_adj = NA, pval_diag_adj = NA, r_squared = NA,
        stringsAsFactors = FALSE
      ))
    }

    tryCatch({
      mod_u <- lm(rna ~ diagnosis_std + consent_age + sex + site_name + antibiotics_yn, data = df)
      co_u  <- summary(mod_u)$coefficients
      du    <- grep("^diagnosis_std", rownames(co_u))

      mod_a <- lm(rna ~ dna + diagnosis_std + consent_age + sex + site_name + antibiotics_yn, data = df)
      co_a  <- summary(mod_a)$coefficients
      da    <- grep("^diagnosis_std", rownames(co_a))
      dna_r <- which(rownames(co_a) == "dna")

      data.frame(
        ec_id         = ec,
        n_samples     = nrow(df),
        coef_unadj    = if (length(du) > 0) co_u[du[1], "Estimate"] else NA,
        pval_unadj    = if (length(du) > 0) co_u[du[1], "Pr(>|t|)"] else NA,
        coef_dna      = if (length(dna_r) > 0) co_a[dna_r[1], "Estimate"] else NA,
        pval_dna      = if (length(dna_r) > 0) co_a[dna_r[1], "Pr(>|t|)"] else NA,
        coef_diag_adj = if (length(da) > 0) co_a[da[1], "Estimate"] else NA,
        stderr_adj    = if (length(da) > 0) co_a[da[1], "Std. Error"] else NA,
        pval_diag_adj = if (length(da) > 0) co_a[da[1], "Pr(>|t|)"] else NA,
        r_squared     = summary(mod_a)$r.squared,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        ec_id = ec, n_samples = nrow(df),
        coef_unadj = NA, pval_unadj = NA,
        coef_dna = NA, pval_dna = NA,
        coef_diag_adj = NA, stderr_adj = NA, pval_diag_adj = NA, r_squared = NA,
        stringsAsFactors = FALSE
      )
    })
  }

  res <- do.call(rbind, parallel::mclapply(seq_along(ecs), fit_one_ec, mc.cores = n_cores))
  res$qval_unadj    <- p.adjust(res$pval_unadj, method = "BH")
  res$qval_diag_adj <- p.adjust(res$pval_diag_adj, method = "BH")

  res <- res %>%
    mutate(
      regulation_status = case_when(
        !is.na(qval_diag_adj) & qval_diag_adj < 0.05 ~ "Active Transcriptional Regulation",
        !is.na(qval_unadj) & qval_unadj < 0.05 & (is.na(qval_diag_adj) | qval_diag_adj >= 0.05) ~ "Abundance Driven (Passive)",
        TRUE ~ "NS"
      )
    ) %>%
    arrange(pval_diag_adj)

  res
}

n_cores <- min(32L, parallel::detectCores(logical = FALSE))
cd_adj <- run_ec_comparison("CD", sample_map, mtx_log, mgx_log, n_cores)
uc_adj <- run_ec_comparison("UC", sample_map, mtx_log, mgx_log, n_cores)

write.csv(cd_adj, file.path(diff_dir, "CD_vs_nonIBD_MGX_adjusted_results.csv"), row.names = FALSE, quote = TRUE)
write.csv(uc_adj, file.path(diff_dir, "UC_vs_nonIBD_MGX_adjusted_results.csv"), row.names = FALSE, quote = TRUE)
cat("Saved CD and UC MGX-adjusted EC results to:", diff_dir, "\n")

summary_df <- bind_rows(
  cd_adj %>% count(regulation_status) %>% mutate(comparison = "CD_vs_nonIBD", pct = round(100 * n / sum(n), 1)),
  uc_adj %>% count(regulation_status) %>% mutate(comparison = "UC_vs_nonIBD", pct = round(100 * n / sum(n), 1))
) %>% select(comparison, regulation_status, n, pct)

write.csv(summary_df, file.path(diff_dir, "transcriptional_activity_summary.csv"), row.names = FALSE, quote = TRUE)

## ---- 7. Visualisations ------------------------------------------------------
cat("\n--- Generating DNA adjustment visualisations ---\n")

reg_colours <- c(
  "Active Transcriptional Regulation" = "#C0392B",
  "Abundance Driven (Passive)"        = "#3498DB",
  "NS"                                = "grey70"
)

# 7a. Scatter: Unadjusted vs DNA-adjusted effect sizes (CD)
p_scatter <- ggplot(cd_adj %>% filter(!is.na(coef_diag_adj), !is.na(coef_unadj)),
                    aes(x = coef_unadj, y = coef_diag_adj, colour = regulation_status)) +
  geom_point(alpha = 0.55, size = 1.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  scale_colour_manual(values = reg_colours, name = "Regulation Category") +
  labs(
    title    = "Unadjusted vs DNA-Adjusted EC Effect Sizes (CD vs nonIBD)",
    subtitle = "Red = active transcriptional regulation independent of DNA abundance",
    x        = "Unadjusted MTX coefficient (log2)",
    y        = "DNA-Adjusted MTX coefficient (log2)"
  ) +
  theme_ibd +
  theme(legend.position = "bottom")

tryCatch(print(p_scatter), error = function(e) invisible(NULL))
ggsave(file.path(figs_dir, "rna_vs_dna_scatter_CD.png"), plot = p_scatter, width = 8, height = 7, dpi = 300)

# 7b. Volcano plots for CD and UC
make_volcano <- function(df, comp_label, out_name) {
  p <- ggplot(df %>% filter(!is.na(pval_diag_adj)),
              aes(x = coef_diag_adj,
                  y = -log10(pmax(pval_diag_adj, 1e-30)),
                  colour = regulation_status)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey30") +
    scale_colour_manual(values = reg_colours, name = "Regulation Category") +
    labs(
      title    = paste0(comp_label, ": Metagenomically Adjusted EC Transcription"),
      subtitle = "Distinguishing True Transcriptional Regulation from DNA Abundance Shifts",
      x        = "DNA-Adjusted Coefficient (log2 fold-change in RNA per DNA copy)",
      y        = expression(-log[10](p-value))
    ) +
    theme_ibd +
    theme(legend.position = "bottom")

  tryCatch(print(p), error = function(e) invisible(NULL))
  ggsave(file.path(figs_dir, out_name), plot = p, width = 9, height = 6, dpi = 300)
  cat("Saved:", file.path("figures", out_name), "\n")
}

make_volcano(cd_adj, "CD vs nonIBD", "volcano_MGX_adjusted_CD.png")
make_volcano(uc_adj, "UC vs nonIBD", "volcano_MGX_adjusted_UC.png")

if (file.exists("Rplots.pdf")) invisible(file.remove("Rplots.pdf"))

## ---- 8. Summary & Session Log -----------------------------------------------
cat("\n=== SUMMARY OF TRANSCRIPTIONAL ACTIVITY (EC-LEVEL) ===\n")
print(summary_df)

log_file <- file.path(logs_dir, sprintf("04_session_%s.txt", format(Sys.Date(), "%Y%m%d")))
sink(log_file)
cat("Script 04: MGX Adjustment & Transcriptional Activity (EC-level)\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(summary_df)
cat("\n")
print(sessionInfo())
sink()
cat("\nSession log saved to:", log_file, "\n")

cat("\n=============================================================\n")
cat("SCRIPT 04 COMPLETE\n")
cat("=============================================================\n")
cat("Session ended:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
