## =============================================================================
## Script:  01_download_and_audit.R
## Project: IBD Microbial Transcript Activity & Candidate Proteins
## Purpose: Download the MTX gene-family table, pathway table, and metadata.
##          Inspect dimensions, identifiers, value distributions, and sample
##          overlap WITHOUT modifying any data.
## Usage:   Source interactively in RStudio. Read the printed output carefully
##          and share it before running any downstream analysis script.
## Author:  [Your name]
## Date:    2026-10-29
## =============================================================================

## ---- 0. Session setup -------------------------------------------------------

# Record start time and R version for reproducibility
cat("=============================================================\n")
cat("Session started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat("=============================================================\n\n")

# Required packages — install only if missing
required_pkgs <- c("data.table", "dplyr", "tidyr", "stringr", "ggplot2", "tools")

missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(data.table)   # fast file reading for large tables
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(tools)        # for md5sum
})

## ---- 1. Project paths -------------------------------------------------------
# Adjust base_dir if your project lives elsewhere
base_dir    <- file.path(path.expand("~"), "Desktop", "607", "IBD_MTX_project")
data_raw    <- file.path(base_dir, "data", "raw")
data_meta   <- file.path(base_dir, "data", "metadata")
results_qc  <- file.path(base_dir, "results", "qc")
logs_dir    <- file.path(base_dir, "logs")

# Create directories if needed
for (d in c(data_raw, data_meta, results_qc, logs_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## ---- 2. File URLs -----------------------------------------------------------
urls <- list(
  genefamilies = paste0(
    "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/",
    "HMP2/MTX/2017-12-14/genefamilies.tsv.gz"
  ),
  pathabundance = paste0(
    "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/",
    "HMP2/MTX/2017-12-14/pathabundance_relab.tsv.gz"
  ),
  metadata = paste0(
    "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/metadata/",
    "hmp2_metadata_2018-08-20.csv"
  )
)

dest_files <- list(
  genefamilies = file.path(data_raw,  "MTX_genefamilies.tsv.gz"),
  pathabundance = file.path(data_raw,  "MTX_pathabundance_relab.tsv.gz"),
  metadata      = file.path(data_meta, "hmp2_metadata_2018-08-20.csv")
)

## ---- 3. Download with integrity check ---------------------------------------
download_file_safe <- function(url, dest, overwrite = FALSE) {
  if (file.exists(dest) && !overwrite) {
    message("  [SKIP] Already exists: ", basename(dest))
    return(invisible(dest))
  }
  message("  [DOWNLOAD] ", basename(dest), " ...")
  tryCatch({
    download.file(url, dest, mode = "wb", quiet = FALSE)
    if (!file.exists(dest) || file.size(dest) < 1000) {
      stop("Download failed or file is suspiciously small.")
    }
    message("  [OK] Downloaded: ", basename(dest),
            " (", round(file.size(dest) / 1e6, 1), " MB)")
  }, error = function(e) {
    stop("Download error for ", basename(dest), ": ", conditionMessage(e))
  })
  invisible(dest)
}

cat("--- Downloading files (skip if already present) ---\n")
for (nm in names(urls)) {
  download_file_safe(urls[[nm]], dest_files[[nm]])
}

## Record MD5 checksums for reproducibility
cat("\n--- MD5 checksums ---\n")
md5s <- sapply(dest_files, function(f) {
  if (file.exists(f)) md5sum(f) else "FILE_MISSING"
})
print(md5s)

## ---- 4. Inspect file sizes --------------------------------------------------
cat("\n--- File sizes ---\n")
for (nm in names(dest_files)) {
  f <- dest_files[[nm]]
  cat(sprintf("  %-40s  %7.1f MB\n", basename(f),
              ifelse(file.exists(f), file.size(f) / 1e6, NA)))
}

## ---- 5. Read metadata (small file, read fully) ------------------------------
cat("\n=============================================================\n")
cat("SECTION A: METADATA INSPECTION\n")
cat("=============================================================\n")

meta <- read.csv(dest_files$metadata, stringsAsFactors = FALSE, check.names = FALSE)
cat("\nMetadata dimensions (rows × columns):", nrow(meta), "×", ncol(meta), "\n")
cat("\nAll column names:\n")
print(colnames(meta))

cat("\nFirst 3 rows × first 20 columns:\n")
print(head(meta[, seq_len(min(20, ncol(meta)))], 3))

cat("\nClass of each column (first 30):\n")
col_classes <- sapply(meta[, seq_len(min(30, ncol(meta)))], class)
print(col_classes)

# --- Identify key columns (case-insensitive search) --------------------------
find_col <- function(df, patterns) {
  cols <- colnames(df)
  for (pat in patterns) {
    hits <- cols[grepl(pat, cols, ignore.case = TRUE)]
    if (length(hits) > 0) return(hits)
  }
  return(character(0))
}

cat("\n--- Candidate key-column matches ---\n")
key_searches <- list(
  diagnosis       = c("diagnosis", "disease", "condition", "group"),
  participant_id  = c("participant", "subject", "external_id", "patient", "person"),
  sample_id       = c("sample_id", "sample.id", "sampleID", "^id$"),
  visit_num       = c("visit", "week", "collection", "timepoint"),
  collection_date = c("date", "collect"),
  disease_activity = c("activity", "cdai", "hbi", "harvey", "mayo", "pcdai"),
  age             = c("age"),
  sex             = c("sex", "gender"),
  site            = c("site", "cohort", "center", "location"),
  antibiotics     = c("antibiotic", "antibiot"),
  data_type       = c("data_type", "datatype", "type", "assay", "omic")
)

found_cols <- lapply(key_searches, find_col, df = meta)
for (nm in names(found_cols)) {
  fc <- found_cols[[nm]]
  cat(sprintf("  %-20s → %s\n", nm,
              if (length(fc) > 0) paste(fc, collapse = ", ") else "[NOT FOUND]"))
}

# --- Tabulate diagnosis / data_type if found ----------------------------------
if (length(found_cols$diagnosis) > 0) {
  diag_col <- found_cols$diagnosis[1]
  cat("\n--- Diagnosis distribution (column:", diag_col, ") ---\n")
  print(table(meta[[diag_col]], useNA = "always"))
}

if (length(found_cols$data_type) > 0) {
  dtype_col <- found_cols$data_type[1]
  cat("\n--- Data type distribution (column:", dtype_col, ") ---\n")
  print(table(meta[[dtype_col]], useNA = "always"))
}

if (length(found_cols$participant_id) > 0) {
  pid_col <- found_cols$participant_id[1]
  cat("\n--- Unique participants:", length(unique(meta[[pid_col]])), "\n")
  cat("--- Total rows (samples):", nrow(meta), "\n")
  spp <- table(meta[[pid_col]])
  cat("--- Samples per participant (summary):\n")
  print(summary(as.integer(spp)))
}

## ---- 6. Inspect gene-family table header (avoid reading full file) ----------
cat("\n=============================================================\n")
cat("SECTION B: MTX GENE-FAMILY TABLE INSPECTION\n")
cat("=============================================================\n")

cat("\nReading first 20 lines of compressed gene-family table ...\n")
gf_header_raw <- tryCatch(
  readLines(gzcon(url(dest_files$genefamilies)), n = 5),
  error = function(e) {
    # Fallback: read from local file
    con <- gzfile(dest_files$genefamilies, "rt")
    on.exit(close(con))
    readLines(con, n = 5)
  }
)
cat("First 5 lines (raw):\n")
cat(paste(substr(gf_header_raw, 1, 200), collapse = "\n"), "\n")

# Parse the header line to get column count
header_line <- gf_header_raw[1]
header_cols <- strsplit(header_line, "\t")[[1]]
cat("\nNumber of columns (= 1 feature column + N sample columns):", length(header_cols), "\n")
cat("First column name (expected: '# Gene Family' or similar):\n  '", header_cols[1], "'\n", sep = "")
cat("\nSample ID examples (first 5 column names after feature column):\n")
cat(paste(" ", header_cols[seq(2, min(6, length(header_cols)))]), sep = "\n")

# Check if column names have abundance-unit suffixes
has_suffix <- any(grepl("_Abundance|_RPK|_relab", header_cols, ignore.case = TRUE))
cat("\nColumn names have abundance-unit suffix (_Abundance, _RPKs, etc.):", has_suffix, "\n")
if (has_suffix) {
  cat("NOTE: You will need to strip these suffixes to match sample IDs to metadata.\n")
}

## ---- 7. Read first 500 rows of gene-family table for structural checks ------
cat("\nReading first 500 feature rows for structural inspection ...\n")
gf_head <- tryCatch({
  con <- gzfile(dest_files$genefamilies, "rt")
  on.exit(close(con))
  fread(text = paste(readLines(con, n = 501), collapse = "\n"),
        sep = "\t", header = TRUE, fill = TRUE)
}, error = function(e) {
  message("fread failed, trying read.table: ", conditionMessage(e))
  con <- gzfile(dest_files$genefamilies, "rt")
  on.exit(close(con))
  read.table(con, sep = "\t", header = TRUE, nrows = 500,
             quote = "", fill = TRUE, stringsAsFactors = FALSE, comment.char = "")
})

cat("Preview table dimensions (should be 500 rows × N+1 columns):", dim(gf_head), "\n")
cat("\nFeature column values — first 10 entries:\n")
feat_col_name <- colnames(gf_head)[1]
cat("  Feature column name:", feat_col_name, "\n")
print(head(gf_head[[feat_col_name]], 10))

# --- Identify feature types --------------------------------------------------
feat_ids <- gf_head[[feat_col_name]]

# Patterns in HUMAnN2 output:
# Unstratified: "UniRef90_A0A000|name" (no pipe after first pipe) — actually "UniRef90_XXXX: name"
# Stratified:   "UniRef90_A0A000|Bacteroides_fragilis"
# Special rows:  "UNMAPPED", "UNINTEGRATED"
is_stratified  <- grepl("\\|[^|]+$", feat_ids) &
                  !grepl("^UNMAPPED$|^UNINTEGRATED$", feat_ids, ignore.case = TRUE)
is_unmapped    <- grepl("^UNMAPPED$", feat_ids, ignore.case = TRUE)
is_unintegrated <- grepl("^UNINTEGRATED$", feat_ids, ignore.case = TRUE)
is_unstratified <- !is_stratified & !is_unmapped & !is_unintegrated

cat("\n--- Feature type breakdown (first 500 rows) ---\n")
cat(sprintf("  Unstratified (community total):  %d\n", sum(is_unstratified)))
cat(sprintf("  Taxon-stratified (|organism):    %d\n", sum(is_stratified)))
cat(sprintf("  UNMAPPED:                        %d\n", sum(is_unmapped)))
cat(sprintf("  UNINTEGRATED:                    %d\n", sum(is_unintegrated)))

cat("\nSample unstratified feature IDs:\n")
print(head(feat_ids[is_unstratified], 5))
cat("\nSample stratified feature IDs:\n")
print(head(feat_ids[is_stratified], 5))

# --- Check value range -------------------------------------------------------
cat("\n--- Value distribution in first 500 rows (all sample columns) ---\n")
num_cols <- colnames(gf_head)[-1]
vals_sample <- as.matrix(gf_head[, ..num_cols])  # data.table subsetting
cat("Summary of all numeric values:\n")
print(summary(as.vector(vals_sample)))

n_zero   <- sum(vals_sample == 0, na.rm = TRUE)
n_total  <- length(vals_sample)
cat(sprintf("  Fraction zero:        %.1f%%\n", 100 * n_zero / n_total))
cat(sprintf("  Fraction NA/missing:  %.1f%%\n",
            100 * sum(is.na(vals_sample)) / n_total))
cat(sprintf("  Maximum value:        %.2f\n", max(vals_sample, na.rm = TRUE)))
cat(sprintf("  Median non-zero:      %.4f\n",
            median(vals_sample[vals_sample > 0], na.rm = TRUE)))

cat("\n*** INTERPRETATION NOTE ***\n")
cat("If values are >> 1 (e.g., hundreds to millions), they are likely RPKs.\n")
cat("If values are between 0 and 1, they may be relative abundances.\n")
cat("The expected unit for HUMAnN2 genefamilies output is RPKs.\n")

## ---- 8. Count total rows (scan file, do NOT load everything) ----------------
cat("\n--- Counting total rows in gene-family table (may take a minute) ---\n")
n_lines <- tryCatch({
  con <- gzfile(dest_files$genefamilies, "rt")
  on.exit(close(con))
  nl <- 0L
  repeat {
    chunk <- readLines(con, n = 50000L)
    if (length(chunk) == 0L) break
    nl <- nl + length(chunk)
  }
  nl
}, error = function(e) {
  message("Could not count lines: ", conditionMessage(e))
  NA_integer_
})
n_features <- n_lines - 1L  # subtract header
cat(sprintf("  Total lines (including header): %d\n", n_lines))
cat(sprintf("  Total feature rows:             %d\n", n_features))

## ---- 9. Inspect pathway table header ----------------------------------------
cat("\n=============================================================\n")
cat("SECTION C: MTX PATHWAY TABLE INSPECTION\n")
cat("=============================================================\n")

pa_header_raw <- tryCatch({
  con <- gzfile(dest_files$pathabundance, "rt")
  on.exit(close(con))
  readLines(con, n = 5)
}, error = function(e) NA_character_)

if (!is.na(pa_header_raw[1])) {
  cat("First 5 lines of pathway table (first 200 chars each):\n")
  cat(paste(substr(pa_header_raw, 1, 200), collapse = "\n"), "\n")

  pa_header_cols <- strsplit(pa_header_raw[1], "\t")[[1]]
  cat("\nNumber of pathway table columns:", length(pa_header_cols), "\n")
  cat("Pathway table sample columns (first 5):\n")
  cat(paste(" ", pa_header_cols[seq(2, min(6, length(pa_header_cols)))]), sep = "\n")
}

## ---- 10. Sample overlap: MTX table vs metadata ------------------------------
cat("\n=============================================================\n")
cat("SECTION D: SAMPLE ID OVERLAP (MTX ↔ METADATA)\n")
cat("=============================================================\n")

# Get sample IDs from the gene-family table (strip any suffix)
gf_sample_ids_raw <- header_cols[-1]  # exclude feature column

# Strip common HUMAnN2 suffixes
strip_suffix <- function(ids) {
  ids <- sub("_Abundance-RPKs$",  "", ids)
  ids <- sub("_Abundance$",       "", ids)
  ids <- sub("_GeneFamilies$",    "", ids)
  ids <- sub("_Abundance-RELAB$", "", ids)
  ids
}
gf_sample_ids <- strip_suffix(gf_sample_ids_raw)
cat("Gene-family table sample count (columns):", length(gf_sample_ids), "\n")
cat("Sample ID examples (after suffix removal):\n")
print(head(gf_sample_ids, 8))

# Get sample IDs from metadata
# Try common ID column names
meta_id_candidates <- c("External ID", "External.ID", "sample_id", "sample.id",
                        "SampleID", "Sample_ID", found_cols$sample_id)
meta_id_col <- NA
for (cn in meta_id_candidates) {
  if (cn %in% colnames(meta)) {
    meta_id_col <- cn; break
  }
}
if (is.na(meta_id_col)) {
  cat("WARNING: Could not automatically identify metadata sample ID column.\n")
  cat("Please inspect colnames(meta) and set meta_id_col manually.\n")
  meta_sample_ids <- character(0)
} else {
  cat("\nMetadata sample ID column used:", meta_id_col, "\n")
  meta_sample_ids <- meta[[meta_id_col]]
  cat("Metadata sample count:", length(meta_sample_ids), "\n")
}

# Filter metadata for MTX rows only (if data_type column exists)
if (!is.na(meta_id_col) && length(found_cols$data_type) > 0) {
  dtype_col <- found_cols$data_type[1]
  mtx_meta <- meta[grepl("MTX|metatranscript", meta[[dtype_col]], ignore.case = TRUE), ]
  cat("\nMetadata rows labelled as MTX:", nrow(mtx_meta), "\n")
  mtx_meta_ids <- mtx_meta[[meta_id_col]]
} else {
  mtx_meta_ids <- meta_sample_ids
  cat("\nCould not filter by data type; using all metadata sample IDs.\n")
}

# Overlap
in_both <- intersect(gf_sample_ids, mtx_meta_ids)
only_gf  <- setdiff(gf_sample_ids, mtx_meta_ids)
only_meta <- setdiff(mtx_meta_ids, gf_sample_ids)

cat(sprintf("\n  Samples in gene-family table:          %d\n", length(gf_sample_ids)))
cat(sprintf("  MTX samples in metadata:               %d\n", length(mtx_meta_ids)))
cat(sprintf("  Overlap (in both):                     %d\n", length(in_both)))
cat(sprintf("  In gene-family table only:             %d\n", length(only_gf)))
cat(sprintf("  In metadata only:                      %d\n", length(only_meta)))

if (length(only_gf) > 0) {
  cat("\n  Sample IDs in table but not metadata (first 10):\n")
  print(head(only_gf, 10))
}
if (length(only_meta) > 0) {
  cat("\n  Sample IDs in metadata but not table (first 10):\n")
  print(head(only_meta, 10))
}

# Check for _TR (technical replicate) samples in table
tr_samples <- gf_sample_ids_raw[grepl("_TR$", gf_sample_ids_raw)]
cat(sprintf("\n  Technical replicate samples (_TR suffix) in table: %d\n",
            length(tr_samples)))
if (length(tr_samples) > 0) {
  cat("  Examples:", paste(head(tr_samples, 5), collapse = ", "), "\n")
  cat("  ACTION REQUIRED: Decide whether to average, discard, or keep _TR samples.\n")
}

## ---- 11. Check for other data types in metadata (MGX, MPX, MBX, HTX) -------
cat("\n=============================================================\n")
cat("SECTION E: MULTI-OMICS DATA TYPE OVERVIEW\n")
cat("=============================================================\n")

if (length(found_cols$data_type) > 0) {
  dtype_col  <- found_cols$data_type[1]
  pid_col    <- if (length(found_cols$participant_id) > 0)
                  found_cols$participant_id[1] else NULL
  dtype_tbl  <- table(meta[[dtype_col]], useNA = "always")
  cat("All data types in metadata:\n")
  print(dtype_tbl)

  # For each data layer, count samples and unique participants
  layers <- c("MTX", "MGX", "MPX", "MBX", "HTX", "16S", "viromics")
  cat("\nData layer summary:\n")
  cat(sprintf("  %-15s  %-12s  %-18s\n", "Layer", "N_samples",
              if (!is.null(pid_col)) "N_participants" else ""))
  for (lyr in layers) {
    sub_meta <- meta[grepl(lyr, meta[[dtype_col]], ignore.case = TRUE), ]
    if (nrow(sub_meta) > 0) {
      np <- if (!is.null(pid_col)) length(unique(sub_meta[[pid_col]])) else NA
      cat(sprintf("  %-15s  %-12d  %s\n", lyr, nrow(sub_meta),
                  if (!is.na(np)) as.character(np) else ""))
    }
  }

  # Cross-tabulate participants with matched MTX and MGX
  if (!is.null(pid_col)) {
    mtx_pids <- unique(meta[grepl("MTX", meta[[dtype_col]], ignore.case = TRUE), pid_col])
    mgx_pids <- unique(meta[grepl("MGX", meta[[dtype_col]], ignore.case = TRUE), pid_col])
    mpx_pids <- unique(meta[grepl("MPX", meta[[dtype_col]], ignore.case = TRUE), pid_col])
    mbx_pids <- unique(meta[grepl("MBX", meta[[dtype_col]], ignore.case = TRUE), pid_col])
    htx_pids <- unique(meta[grepl("HTX", meta[[dtype_col]], ignore.case = TRUE), pid_col])

    cat(sprintf("\n  Participants with MTX:               %d\n", length(mtx_pids)))
    cat(sprintf("  Participants with MTX + MGX:         %d\n", length(intersect(mtx_pids, mgx_pids))))
    cat(sprintf("  Participants with MTX + MPX:         %d\n", length(intersect(mtx_pids, mpx_pids))))
    cat(sprintf("  Participants with MTX + MBX:         %d\n", length(intersect(mtx_pids, mbx_pids))))
    cat(sprintf("  Participants with MTX + HTX:         %d\n", length(intersect(mtx_pids, htx_pids))))
  }
}

## ---- 12. Save audit log -----------------------------------------------------
audit_log <- file.path(logs_dir, paste0("01_audit_", format(Sys.Date(), "%Y%m%d"), ".txt"))
# Capture relevant summary to a log file
log_content <- c(
  paste0("Audit date: ", Sys.time()),
  paste0("R version: ", R.version.string),
  "",
  paste0("MTX genefamilies size (MB): ", round(file.size(dest_files$genefamilies)/1e6, 1)),
  paste0("MTX pathway table size (MB): ", round(file.size(dest_files$pathabundance)/1e6, 1)),
  paste0("Metadata size (MB): ", round(file.size(dest_files$metadata)/1e6, 3)),
  "",
  paste0("Gene-family table columns: ", length(header_cols)),
  paste0("Gene-family table total feature rows: ", n_features),
  "",
  paste0("Metadata rows: ", nrow(meta)),
  paste0("Metadata columns: ", ncol(meta)),
  "",
  paste0("Overlap (MTX table ↔ MTX metadata): ", length(in_both)),
  paste0("Technical replicate samples (_TR): ", length(tr_samples))
)
writeLines(log_content, audit_log)
cat("\n--- Audit log saved to:", audit_log, "---\n")

## ---- 13. Final message -------------------------------------------------------
cat("\n=============================================================\n")
cat("AUDIT COMPLETE. NEXT STEPS:\n")
cat("=============================================================\n")
cat("1. Review ALL printed output above carefully.\n")
cat("2. Confirm:\n")
cat("   a. The first column name of the gene-family table\n")
cat("      (expected: '# Gene Family')\n")
cat("   b. That values are RPKs (large positive numbers, not fractions)\n")
cat("   c. The exact column name in metadata for diagnosis\n")
cat("   d. The exact column name in metadata for participant/subject ID\n")
cat("   e. The number of samples overlapping between table and metadata\n")
cat("   f. Whether _TR technical replicate samples are present\n")
cat("3. Share this output before running 02_metadata_and_cohort_QC.R\n")
cat("   so that downstream scripts use the CORRECT column names.\n")
cat("=============================================================\n")
cat("Session ended:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
