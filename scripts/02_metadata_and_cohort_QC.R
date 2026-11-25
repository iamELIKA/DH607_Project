## =============================================================================
## Script:  02_metadata_and_cohort_QC.R
## Project: IBD Microbial Transcript Activity & Candidate Proteins
## Purpose: Parse and characterise the study metadata; define the analysis
##          cohort; report sample counts, repeated-measures structure, and
##          multi-omics data-type overlap.
##          Column names are hard-coded from confirmed output of 01_download_and_audit.R
## Requires: data/metadata/hmp2_metadata_2018-08-20.csv (downloaded by script 01)
## Outputs:
##   results/qc/cohort_summary.csv
##   results/qc/sample_overlap.csv
##   results/qc/metadata_QC_report.txt
##   figures/qc_diagnosis_distribution.pdf
##   figures/qc_samples_per_participant.pdf
##   figures/qc_multiomics_overlap.pdf
## =============================================================================

set.seed(42)

cat("=============================================================\n")
cat("Script 02: Metadata & Cohort QC\n")
cat("Session started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat("=============================================================\n\n")

## ---- 0. Packages ------------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "stringr", "scales", "readr",
                   "forcats", "patchwork")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  message("Installing: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr)
  library(scales); library(readr); library(forcats); library(patchwork)
})

## ---- 1. Paths ---------------------------------------------------------------
base_dir   <- file.path(path.expand("~"), "Desktop", "607", "IBD_MTX_project")
meta_file  <- file.path(base_dir, "data", "metadata", "hmp2_metadata_2018-08-20.csv")
results_qc <- file.path(base_dir, "results", "qc")
figs_dir   <- file.path(base_dir, "figures")
logs_dir   <- file.path(base_dir, "logs")
proc_dir   <- file.path(base_dir, "data", "processed")

for (d in c(results_qc, figs_dir, logs_dir, proc_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## ---- 2. Hard-coded column names (confirmed from script 01 output) -----------
## If these need updating, change ONLY here; nothing else in this script uses
## raw column-name strings.
COL_SAMPLE_ID   <- "External ID"
COL_PART_ID     <- "Participant ID"
COL_DIAGNOSIS   <- "diagnosis"
COL_DATA_TYPE   <- "data_type"
COL_AGE         <- "consent_age"
COL_SEX         <- "sex"
COL_SITE        <- "site_name"
COL_VISIT_NUM   <- "visit_num"
COL_WEEK_NUM    <- "week_num"
COL_ANTIBIOTICS <- "Antibiotics"   # 2-week recall binary
COL_HBI         <- "hbi"           # Harvey-Bradshaw Index (CD activity)
COL_SCCAI       <- "sccai"         # Simple Clinical Colitis Activity Index (UC)
COL_IS_INFLAMED <- "is_inflamed"
COL_FECALCAL    <- "fecalcal_ng_ml"
COL_IMMUNOSUPP  <- "Immunosuppressants (e.g. oral corticosteroids)"
COL_REMICADE    <- "Remicade (Infliximab)"
COL_HUMIRA      <- "Humira (Adalimumab)"

# Data-type labels expected in COL_DATA_TYPE
DTYPE_MTX <- "metatranscriptomics"   # may differ — audit will check
DTYPE_MGX <- "metagenomics"
DTYPE_MPX <- "proteomics"
DTYPE_MBX <- "metabolomics"
DTYPE_HTX <- "host_transcriptomics"

## ---- 3. Load metadata -------------------------------------------------------
cat("--- Loading metadata ---\n")
if (!file.exists(meta_file)) stop("Metadata file not found. Run script 01 first.")

meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("Dimensions:", nrow(meta), "rows ×", ncol(meta), "columns\n\n")

## Verify all expected columns are present
expected_cols <- c(COL_SAMPLE_ID, COL_PART_ID, COL_DIAGNOSIS, COL_DATA_TYPE,
                   COL_AGE, COL_SEX, COL_SITE, COL_VISIT_NUM, COL_WEEK_NUM)
missing_expected <- setdiff(expected_cols, colnames(meta))
if (length(missing_expected) > 0) {
  warning("These expected columns are MISSING from metadata: ",
          paste(missing_expected, collapse = ", "))
} else {
  cat("All primary columns confirmed present.\n\n")
}

## ---- 4. Inspect data_type values -------------------------------------------
cat("=== Raw data_type values ===\n")
dtype_raw <- table(meta[[COL_DATA_TYPE]], useNA = "always")
print(dtype_raw)
cat("\n")

# Map actual data_type strings to short labels
# Adjust these pattern strings if the table above shows different names
dtype_patterns <- list(
  MTX = c("metatranscriptomics", "MTX"),
  MGX = c("metagenomics", "MGX"),
  MPX = c("proteomics", "MPX"),
  MBX = c("metabolomics", "MBX"),
  HTX = c("host_transcriptomics", "HTX", "host transcriptomics"),
  biopsy = c("biopsy", "RNA"),
  viromics = c("viromics"),
  serology = c("serology"),
  `16S` = c("16S", "amplicon")
)

classify_dtype <- function(x, patterns) {
  for (nm in names(patterns)) {
    if (any(grepl(paste(patterns[[nm]], collapse = "|"), x, ignore.case = TRUE))) {
      return(nm)
    }
  }
  return("other")
}
meta$dtype_label <- sapply(meta[[COL_DATA_TYPE]], classify_dtype, patterns = dtype_patterns)
cat("Classified data_type labels:\n")
print(table(meta$dtype_label, useNA = "always"))
cat("\n")

## ---- 5. Inspect diagnosis values -------------------------------------------
cat("=== diagnosis column values ===\n")
diag_raw <- table(meta[[COL_DIAGNOSIS]], useNA = "always")
print(diag_raw)
cat("\n")

# Standardise diagnosis labels
meta$diagnosis_std <- case_when(
  grepl("^CD$|Crohn|crohn", meta[[COL_DIAGNOSIS]], ignore.case = TRUE) ~ "CD",
  grepl("^UC$|ulcerative|Ulcerative", meta[[COL_DIAGNOSIS]], ignore.case = TRUE) ~ "UC",
  grepl("nonIBD|non.IBD|non_IBD|control|Control|healthy", meta[[COL_DIAGNOSIS]],
        ignore.case = TRUE) ~ "nonIBD",
  is.na(meta[[COL_DIAGNOSIS]]) ~ NA_character_,
  TRUE ~ "other"
)
cat("Standardised diagnosis labels:\n")
print(table(meta$diagnosis_std, useNA = "always"))
cat("\n")

## ---- 6. MTX-specific subset -------------------------------------------------
cat("=== MTX sample characterisation ===\n")
mtx_meta <- meta %>% filter(dtype_label == "MTX")
cat("Total MTX rows:", nrow(mtx_meta), "\n")

# Unique MTX participants
mtx_participants <- unique(mtx_meta[[COL_PART_ID]])
cat("Unique MTX participants:", length(mtx_participants), "\n\n")

# Diagnosis breakdown for MTX samples
cat("MTX samples by diagnosis:\n")
mtx_diag_tbl <- table(mtx_meta$diagnosis_std, useNA = "always")
print(mtx_diag_tbl)
cat("\n")

cat("MTX participants by diagnosis:\n")
mtx_pid_diag <- mtx_meta %>%
  select(all_of(c(COL_PART_ID, "diagnosis_std"))) %>%
  distinct()
print(table(mtx_pid_diag$diagnosis_std, useNA = "always"))
cat("\n")

## ---- 7. Repeated-measures structure -----------------------------------------
cat("=== Repeated-measures structure (MTX only) ===\n")
spp <- mtx_meta %>%
  group_by(.data[[COL_PART_ID]]) %>%
  summarise(n_samples = n(), .groups = "drop")

cat("Samples per participant (MTX):\n")
print(summary(spp$n_samples))
cat("\n")
cat("Distribution:\n")
print(table(spp$n_samples))
cat("\n")

# Participants with only 1 sample (cannot contribute to within-person random effect well)
single_sample_pids <- spp %>% filter(n_samples == 1)
cat("Participants with only 1 MTX sample:", nrow(single_sample_pids),
    "(will still be included; random effect will be estimated from between-person variance)\n\n")

## ---- 8. Sample-size check for primary comparison ----------------------------
cat("=== Primary comparison: CD vs nonIBD ===\n")
cd_meta   <- mtx_meta %>% filter(diagnosis_std == "CD")
ctrl_meta <- mtx_meta %>% filter(diagnosis_std == "nonIBD")
uc_meta   <- mtx_meta %>% filter(diagnosis_std == "UC")

cat(sprintf("CD:     %d samples from %d participants\n",
            nrow(cd_meta), length(unique(cd_meta[[COL_PART_ID]]))))
cat(sprintf("UC:     %d samples from %d participants\n",
            nrow(uc_meta), length(unique(uc_meta[[COL_PART_ID]]))))
cat(sprintf("nonIBD: %d samples from %d participants\n",
            nrow(ctrl_meta), length(unique(ctrl_meta[[COL_PART_ID]]))))
cat("\nMinimum recommended for longitudinal MaAsLin2: 30+ samples per group.\n")

## ---- 9. Covariate completeness in MTX samples --------------------------------
cat("\n=== Covariate completeness (MTX samples only) ===\n")
covariate_cols <- c(COL_AGE, COL_SEX, COL_SITE, COL_ANTIBIOTICS,
                    COL_HBI, COL_SCCAI, COL_IS_INFLAMED, COL_FECALCAL)
cov_available  <- intersect(covariate_cols, colnames(mtx_meta))

cov_summary <- data.frame(
  column      = cov_available,
  n_non_na    = sapply(cov_available, function(cn) sum(!is.na(mtx_meta[[cn]]) & mtx_meta[[cn]] != "")),
  pct_complete = NA_real_,
  stringsAsFactors = FALSE
)
cov_summary$pct_complete <- round(100 * cov_summary$n_non_na / nrow(mtx_meta), 1)
cat("Covariate completeness:\n")
print(cov_summary)
cat("\n")

## ---- 10. Sex and age distributions ------------------------------------------
if (COL_SEX %in% colnames(mtx_meta)) {
  cat("Sex distribution (MTX participants):\n")
  sex_tbl <- mtx_meta %>%
    select(all_of(c(COL_PART_ID, COL_SEX, "diagnosis_std"))) %>%
    distinct(across(all_of(c(COL_PART_ID, "diagnosis_std"))), .keep_all = TRUE)
  print(table(sex_tbl[[COL_SEX]], sex_tbl$diagnosis_std, useNA = "always"))
  cat("\n")
}

if (COL_AGE %in% colnames(mtx_meta)) {
  mtx_meta[[COL_AGE]] <- suppressWarnings(as.numeric(mtx_meta[[COL_AGE]]))
  age_per_pid <- mtx_meta %>%
    select(all_of(c(COL_PART_ID, COL_AGE, "diagnosis_std"))) %>%
    distinct(.data[[COL_PART_ID]], .keep_all = TRUE)
  cat("Age summary (MTX participants, per diagnosis):\n")
  age_smry <- age_per_pid %>%
    group_by(diagnosis_std) %>%
    summarise(
      n       = n(),
      median  = median(.data[[COL_AGE]], na.rm = TRUE),
      IQR_lo  = quantile(.data[[COL_AGE]], 0.25, na.rm = TRUE),
      IQR_hi  = quantile(.data[[COL_AGE]], 0.75, na.rm = TRUE),
      missing = sum(is.na(.data[[COL_AGE]])),
      .groups = "drop"
    )
  print(age_smry)
  cat("\n")
}

## ---- 11. Site / cohort distribution -----------------------------------------
if (COL_SITE %in% colnames(mtx_meta)) {
  cat("Site distribution (MTX samples):\n")
  print(table(mtx_meta[[COL_SITE]], useNA = "always"))
  cat("\n")
}

## ---- 12. Disease-activity score distribution --------------------------------
cat("=== Disease activity scores ===\n")
if (COL_HBI %in% colnames(mtx_meta)) {
  cd_hbi <- mtx_meta %>% filter(diagnosis_std == "CD") %>%
    pull(.data[[COL_HBI]]) %>% suppressWarnings(as.numeric())
  cat("HBI in CD samples — summary:\n")
  print(summary(cd_hbi))
  cat(sprintf("  HBI ≥ 5 (active CD): %d / %d samples (%.1f%%)\n",
              sum(cd_hbi >= 5, na.rm = TRUE), sum(!is.na(cd_hbi)),
              100 * mean(cd_hbi >= 5, na.rm = TRUE)))
  cat("\n")
}

if (COL_SCCAI %in% colnames(mtx_meta)) {
  uc_sccai <- mtx_meta %>% filter(diagnosis_std == "UC") %>%
    pull(.data[[COL_SCCAI]]) %>% suppressWarnings(as.numeric())
  cat("SCCAI in UC samples — summary:\n")
  print(summary(uc_sccai))
  cat(sprintf("  SCCAI ≥ 5 (active UC): %d / %d samples (%.1f%%)\n",
              sum(uc_sccai >= 5, na.rm = TRUE), sum(!is.na(uc_sccai)),
              100 * mean(uc_sccai >= 5, na.rm = TRUE)))
  cat("\n")
}

## ---- 13. Antibiotic use in MTX samples --------------------------------------
if (COL_ANTIBIOTICS %in% colnames(mtx_meta)) {
  abx_vals <- mtx_meta[[COL_ANTIBIOTICS]]
  cat("Antibiotic use distribution (MTX samples):\n")
  print(table(abx_vals, mtx_meta$diagnosis_std, useNA = "always"))
  cat("\n")
  # Create a clean binary flag
  mtx_meta$antibiotics_yn <- case_when(
    grepl("^yes$", abx_vals, ignore.case = TRUE) ~ 1L,
    grepl("^no$",  abx_vals, ignore.case = TRUE) ~ 0L,
    TRUE ~ NA_integer_
  )
  cat(sprintf("Antibiotic flag: %d Yes, %d No, %d missing\n",
              sum(mtx_meta$antibiotics_yn == 1, na.rm = TRUE),
              sum(mtx_meta$antibiotics_yn == 0, na.rm = TRUE),
              sum(is.na(mtx_meta$antibiotics_yn))))
}

## ---- 14. Multi-omics overlap ------------------------------------------------
cat("\n=== Multi-omics data-layer overlap ===\n")

get_layer_info <- function(df, pid_col, dtype_lbl, target_label) {
  sub <- df %>% filter(dtype_label == target_label)
  list(
    label        = target_label,
    n_samples    = nrow(sub),
    participants = unique(sub[[pid_col]])
  )
}

layers <- c("MTX", "MGX", "MPX", "MBX", "HTX")
layer_info <- lapply(layers, function(lbl) get_layer_info(meta, COL_PART_ID, dtype_label, lbl))
names(layer_info) <- layers

cat(sprintf("  %-12s  %-12s  %-20s\n", "Layer", "N_samples", "N_participants"))
cat(rep("-", 50), "\n", sep = "")
for (lbl in layers) {
  li <- layer_info[[lbl]]
  cat(sprintf("  %-12s  %-12d  %-20d\n",
              lbl, li$n_samples, length(li$participants)))
}

cat("\nParticipants with matched layers (intersection with MTX):\n")
mtx_pids <- layer_info$MTX$participants
for (lbl in setdiff(layers, "MTX")) {
  overlap_pids  <- intersect(mtx_pids, layer_info[[lbl]]$participants)
  cat(sprintf("  MTX ∩ %-6s: %d participants\n", lbl, length(overlap_pids)))
}

# Sample-level overlap (same External ID in both MTX and MGX metadata)
# Note: same participant does NOT guarantee same-visit samples
cat("\nNOTE: Participant overlap above does not guarantee same-visit pairing.\n")
cat("Same-visit matching requires merging on Participant ID + week_num/visit_num.\n")

# Attempt same-visit overlap for MTX vs MGX
if ("MGX" %in% names(layer_info) && length(layer_info$MGX$participants) > 0) {
  mtx_visits <- mtx_meta %>%
    select(all_of(c(COL_PART_ID, COL_WEEK_NUM, COL_SAMPLE_ID))) %>%
    rename(pid = 1, week = 2, mtx_id = 3)

  mgx_meta <- meta %>%
    filter(dtype_label == "MGX") %>%
    select(all_of(c(COL_PART_ID, COL_WEEK_NUM, COL_SAMPLE_ID))) %>%
    rename(pid = 1, week = 2, mgx_id = 3)

  same_visit_pairs <- inner_join(mtx_visits, mgx_meta, by = c("pid", "week"))
  cat(sprintf("\nSame-visit MTX–MGX pairs: %d (from %d participants)\n",
              nrow(same_visit_pairs),
              length(unique(same_visit_pairs$pid))))

  # Save same-visit pairs for script 04
  write.csv(same_visit_pairs,
            file.path(proc_dir, "MTX_MGX_same_visit_pairs.csv"),
            row.names = FALSE, quote = FALSE)
  cat("Saved to: data/processed/MTX_MGX_same_visit_pairs.csv\n")
}

## ---- 15. Define analysis cohort and save ------------------------------------
cat("\n=== Defining analysis cohort ===\n")

# Eligibility criteria:
# 1. dtype_label is MTX
# 2. diagnosis_std is CD, UC, or nonIBD (no NA, no 'other')
# 3. Has a valid Participant ID
# 4. Has diagnosis confirmed (not missing)

analysis_cohort <- mtx_meta %>%
  filter(
    dtype_label    == "MTX",
    diagnosis_std  %in% c("CD", "UC", "nonIBD"),
    !is.na(.data[[COL_PART_ID]]),
    !is.na(diagnosis_std)
  )

cat(sprintf("Analysis cohort: %d MTX samples from %d participants\n",
            nrow(analysis_cohort),
            length(unique(analysis_cohort[[COL_PART_ID]]))))
cat("\nBy diagnosis:\n")
print(table(analysis_cohort$diagnosis_std))

# Save cohort metadata
cohort_save_path <- file.path(proc_dir, "analysis_cohort_metadata.csv")
write.csv(analysis_cohort, cohort_save_path, row.names = FALSE, quote = TRUE)
cat("Saved cohort metadata to:", cohort_save_path, "\n")

## ---- 16. Summary report table -----------------------------------------------
summary_tbl <- analysis_cohort %>%
  group_by(diagnosis_std) %>%
  summarise(
    n_samples      = n(),
    n_participants = n_distinct(.data[[COL_PART_ID]]),
    median_samples_per_participant = median(table(.data[[COL_PART_ID]])),
    pct_female     = if (COL_SEX %in% colnames(.)) {
      round(100 * mean(grepl("^F|female", .data[[COL_SEX]], ignore.case = TRUE),
                       na.rm = TRUE), 1)
    } else NA_real_,
    median_age     = if (COL_AGE %in% colnames(.)) {
      median(suppressWarnings(as.numeric(.data[[COL_AGE]])), na.rm = TRUE)
    } else NA_real_,
    pct_antibiotics = if ("antibiotics_yn" %in% colnames(.)) {
      round(100 * mean(.data[["antibiotics_yn"]], na.rm = TRUE), 1)
    } else NA_real_,
    .groups = "drop"
  )

cat("\n=== COHORT SUMMARY TABLE ===\n")
print(as.data.frame(summary_tbl))

write.csv(summary_tbl,
          file.path(results_qc, "cohort_summary.csv"),
          row.names = FALSE, quote = FALSE)
cat("Saved to: results/qc/cohort_summary.csv\n")

## ---- 17. Figures ------------------------------------------------------------
cat("\n--- Generating QC figures ---\n")

theme_ibd <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#2C3E50"),
        strip.text = element_text(colour = "white", face = "bold"),
        panel.grid.minor = element_blank())

diag_palette <- c("CD" = "#E74C3C", "UC" = "#3498DB", "nonIBD" = "#2ECC71")

# A) Sample counts by diagnosis
p1 <- analysis_cohort %>%
  count(diagnosis_std) %>%
  ggplot(aes(x = fct_reorder(diagnosis_std, n), y = n, fill = diagnosis_std)) +
  geom_col(width = 0.6, colour = "grey20") +
  geom_text(aes(label = n), hjust = -0.2, size = 4) +
  scale_fill_manual(values = diag_palette, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  coord_flip() +
  labs(title = "MTX samples by diagnosis", x = NULL, y = "Number of samples") +
  theme_ibd

# B) Samples per participant
p2 <- analysis_cohort %>%
  group_by(.data[[COL_PART_ID]], diagnosis_std) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = n, fill = diagnosis_std)) +
  geom_histogram(binwidth = 1, colour = "white", position = "stack") +
  scale_fill_manual(values = diag_palette, name = "Diagnosis") +
  labs(title = "MTX samples per participant",
       x = "Number of samples", y = "Number of participants") +
  theme_ibd

# C) Multi-omics participant overlap bar chart
overlap_df <- data.frame(
  layer   = layers,
  n_pids  = sapply(layers, function(l) length(layer_info[[l]]$participants)),
  mtx_overlap = sapply(layers, function(l) {
    length(intersect(mtx_pids, layer_info[[l]]$participants))
  })
) %>%
  pivot_longer(cols = c(n_pids, mtx_overlap),
               names_to = "type", values_to = "count") %>%
  mutate(type = recode(type,
                       n_pids = "Total participants",
                       mtx_overlap = "Overlap with MTX"))

p3 <- ggplot(overlap_df, aes(x = layer, y = count, fill = type)) +
  geom_col(position = "dodge", colour = "grey20", width = 0.7) +
  scale_fill_manual(values = c("Total participants" = "#95A5A6",
                               "Overlap with MTX"   = "#2C3E50"),
                    name = NULL) +
  labs(title = "Multi-omics participant coverage",
       x = "Data layer", y = "Number of participants") +
  theme_ibd

# D) Disease activity over time for MTX samples (CD: HBI)
if (COL_HBI %in% colnames(analysis_cohort) && COL_WEEK_NUM %in% colnames(analysis_cohort)) {
  cd_act <- analysis_cohort %>%
    filter(diagnosis_std == "CD") %>%
    mutate(hbi_num   = suppressWarnings(as.numeric(.data[[COL_HBI]])),
           week_num  = suppressWarnings(as.numeric(.data[[COL_WEEK_NUM]]))) %>%
    filter(!is.na(hbi_num), !is.na(week_num))

  if (nrow(cd_act) > 10) {
    p4_hbi <- ggplot(cd_act, aes(x = week_num, y = hbi_num)) +
      geom_point(alpha = 0.3, colour = "#E74C3C", size = 1) +
      geom_smooth(method = "loess", se = TRUE, colour = "#C0392B") +
      geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40") +
      annotate("text", x = max(cd_act$week_num, na.rm = TRUE) * 0.9,
               y = 5.5, label = "Active (HBI ≥ 5)", size = 3, colour = "grey30") +
      labs(title = "CD disease activity (HBI) over time",
           x = "Study week", y = "Harvey-Bradshaw Index") +
      theme_ibd
    ggsave(file.path(figs_dir, "qc_CD_HBI_over_time.pdf"),
           p4_hbi, width = 7, height = 4, device = cairo_pdf)
    cat("Saved: figures/qc_CD_HBI_over_time.pdf\n")
  }
}

# Save figures
ggsave(file.path(figs_dir, "qc_diagnosis_distribution.pdf"),
       p1, width = 6, height = 4, device = cairo_pdf)
ggsave(file.path(figs_dir, "qc_samples_per_participant.pdf"),
       p2, width = 7, height = 4, device = cairo_pdf)
ggsave(file.path(figs_dir, "qc_multiomics_overlap.pdf"),
       p3, width = 7, height = 5, device = cairo_pdf)

cat("Saved: figures/qc_diagnosis_distribution.pdf\n")
cat("Saved: figures/qc_samples_per_participant.pdf\n")
cat("Saved: figures/qc_multiomics_overlap.pdf\n")

## ---- 18. Write QC text report -----------------------------------------------
qc_report_path <- file.path(results_qc, "metadata_QC_report.txt")
sink(qc_report_path)
cat("IBD MTX Project — Metadata QC Report\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("=== COLUMN MAPPING ===\n")
cat("Sample ID:     ", COL_SAMPLE_ID, "\n")
cat("Participant:   ", COL_PART_ID, "\n")
cat("Diagnosis:     ", COL_DIAGNOSIS, "\n")
cat("Data type:     ", COL_DATA_TYPE, "\n")
cat("Age:           ", COL_AGE, "\n")
cat("Sex:           ", COL_SEX, "\n")
cat("Site:          ", COL_SITE, "\n\n")

cat("=== COHORT ===\n")
print(as.data.frame(summary_tbl))

cat("\n=== MULTIOMICS PARTICIPANT COVERAGE ===\n")
for (lbl in layers) {
  li <- layer_info[[lbl]]
  cat(sprintf("  %-6s: %4d samples, %3d participants (%d overlap MTX)\n",
              lbl, li$n_samples, length(li$participants),
              length(intersect(mtx_pids, li$participants))))
}

cat("\n=== STATISTICAL SUITABILITY ===\n")
cat("Primary comparison: CD vs nonIBD\n")
cd_n  <- nrow(cd_meta); cd_p  <- length(unique(cd_meta[[COL_PART_ID]]))
ctr_n <- nrow(ctrl_meta); ctr_p <- length(unique(ctrl_meta[[COL_PART_ID]]))
cat(sprintf("  CD:     %d samples, %d participants\n", cd_n, cd_p))
cat(sprintf("  nonIBD: %d samples, %d participants\n", ctr_n, ctr_p))
cat(sprintf("  Total analysis set: %d samples\n\n", cd_n + ctr_n))
cat("Secondary comparison: UC vs nonIBD\n")
uc_n <- nrow(uc_meta); uc_p <- length(unique(uc_meta[[COL_PART_ID]]))
cat(sprintf("  UC:     %d samples, %d participants\n", uc_n, uc_p))
cat(sprintf("  nonIBD: %d samples, %d participants\n", ctr_n, ctr_p))
sink()
cat("QC report saved to:", qc_report_path, "\n")

## ---- 19. Save session info --------------------------------------------------
session_log <- file.path(logs_dir, paste0("02_session_", format(Sys.Date(), "%Y%m%d"), ".txt"))
sink(session_log)
sessionInfo()
sink()

cat("\n=============================================================\n")
cat("SCRIPT 02 COMPLETE\n")
cat("=============================================================\n")
cat("Key outputs:\n")
cat("  data/processed/analysis_cohort_metadata.csv  ← used by script 03\n")
cat("  data/processed/MTX_MGX_same_visit_pairs.csv  ← used by script 04 (if MGX available)\n")
cat("  results/qc/cohort_summary.csv\n")
cat("  results/qc/metadata_QC_report.txt\n")
cat("  figures/qc_*.pdf\n")
cat("\n")
cat("NEXT STEPS:\n")
cat("  1. Review cohort_summary.csv to confirm sample sizes are adequate.\n")
cat("  2. Note which covariates have missing data; decide on imputation vs exclusion.\n")
cat("  3. Share the printed cohort summary before running script 03.\n")
cat("  4. ALSO share the Section B output from script 01 (gene-family table structure)\n")
cat("     if you have not done so yet — it determines how script 03 loads the data.\n")
cat("=============================================================\n")
cat("Session ended:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
