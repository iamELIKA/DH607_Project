## =============================================================================
## Script:  h1_download_and_visualise_HTX.R (Script h1)
## Project: IBD Microbial Transcript Activity & Candidate Proteins
## Purpose: Download host transcriptomics (HTX) count data from IBDMDB / GEO,
##          perform data cleanliness & QC checks, align with patient metadata,
##          and generate data visualisations:
##            1. Library size / sequencing depth distribution
##            2. Principal Component Analysis (PCA) by diagnosis and biopsy site
##            3. Hierarchical clustering heatmap of top variable host genes
##            4. Bar plots of relative gene expression across disease states
##
## Display: All plots are displayed to the console/active plot device
##          and saved in high-resolution image format (.png, 300 DPI), not PDF.
##
## Outputs:
##   data/raw/host_tx_counts.tsv.gz
##   data/processed/host_tx_log2cpm_filtered.rds
##   results/qc/host_transcriptomics_QC_report.txt
##   figures/htx_qc_library_size.png
##   figures/htx_pca_diagnosis_tissue.png
##   figures/htx_heatmap_top50_genes.png
##   figures/htx_barplot_top_variable_genes.png
## =============================================================================

set.seed(42)

cat("=============================================================\n")
cat("Script h1: Download, QC & Visualise Host Transcriptomics (HTX)\n")
cat("Session started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n")
cat("=============================================================\n\n")

## ---- 0. Packages ------------------------------------------------------------
required_pkgs <- c("data.table", "dplyr", "tidyr", "ggplot2", "stringr",
                   "scales", "pheatmap", "ggrepel")

missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(pheatmap)
  library(ggrepel)
})

## ---- 1. Paths & Directories -------------------------------------------------
base_dir   <- file.path(path.expand("~"), "Desktop", "607", "IBD_MTX_project")
data_raw   <- file.path(base_dir, "data", "raw")
data_meta  <- file.path(base_dir, "data", "metadata")
data_proc  <- file.path(base_dir, "data", "processed")
results_qc <- file.path(base_dir, "results", "qc")
figs_dir   <- file.path(base_dir, "figures")
logs_dir   <- file.path(base_dir, "logs")

for (d in c(data_raw, data_meta, data_proc, results_qc, figs_dir, logs_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

theme_htx <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "#2C3E50"),
        strip.text       = element_text(colour = "white", face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", size = 13))

diag_colors <- c("CD" = "#E74C3C", "UC" = "#3498DB", "nonIBD" = "#2ECC71")

## ---- 2. Download Host Transcriptomics Data ----------------------------------
htx_dest <- file.path(data_raw, "host_tx_counts.tsv.gz")

htx_urls <- c(
  "https://www.ibdmdb.org/downloads/products/HMP2/HTX/host_tx_counts.tsv.gz",
  "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/HTX/host_tx_counts.tsv.gz",
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111889/suppl/GSE111889_host_tx_counts.tsv.gz"
)

if (!file.exists(htx_dest) || file.size(htx_dest) < 1000) {
  cat("--- Downloading host transcriptomics count table ---\n")
  download_success <- FALSE
  for (u in htx_urls) {
    cat("Attempting download from:", u, "\n")
    tryCatch({
      download.file(u, destfile = htx_dest, mode = "wb", quiet = FALSE)
      if (file.exists(htx_dest) && file.size(htx_dest) > 1000) {
        cat("Successfully downloaded:", basename(htx_dest),
            "(", round(file.size(htx_dest) / 1e6, 2), "MB)\n\n")
        download_success <- TRUE
        break
      }
    }, error = function(e) {
      cat("  Mirror failed:", conditionMessage(e), "\n")
    })
  }
  if (!download_success) {
    stop("Could not download host transcriptomics counts from any mirror. Check network connection.")
  }
} else {
  cat("Host transcriptomics file already exists locally:", basename(htx_dest),
      "(", round(file.size(htx_dest) / 1e6, 2), "MB)\n\n")
}

## ---- 3. Load Metadata & Check Alignment -------------------------------------
cat("--- Loading study metadata ---\n")
meta_file <- file.path(data_meta, "hmp2_metadata_2018-08-20.csv")
if (!file.exists(meta_file)) stop("Metadata file missing at: ", meta_file)

meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
colnames(meta) <- make.unique(colnames(meta))

# Subset to host transcriptomics entries
htx_meta <- meta %>%
  filter(data_type == "host_transcriptomics") %>%
  rename(
    Sample_ID      = `External ID`,
    Participant_ID = `Participant ID`
  ) %>%
  mutate(
    diagnosis_std = case_when(
      grepl("^CD$|Crohn", diagnosis, ignore.case = TRUE) ~ "CD",
      grepl("^UC$|Ulcerative", diagnosis, ignore.case = TRUE) ~ "UC",
      grepl("nonIBD|control|healthy", diagnosis, ignore.case = TRUE) ~ "nonIBD",
      TRUE ~ "other"
    ),
    tissue_location = if ("biopsy_location" %in% colnames(.)) {
      case_when(
        grepl("Ileum", biopsy_location, ignore.case = TRUE) ~ "Ileum",
        grepl("Rectum", biopsy_location, ignore.case = TRUE) ~ "Rectum",
        grepl("Colon|Cecum", biopsy_location, ignore.case = TRUE) ~ "Colon",
        TRUE ~ "Other"
      )
    } else { "Biopsy" },
    consent_age = suppressWarnings(as.numeric(consent_age))
  )

cat("HTX metadata records:", nrow(htx_meta), "samples across",
    length(unique(htx_meta$Participant_ID)), "participants\n")
cat("Diagnoses breakdown:\n")
print(table(htx_meta$diagnosis_std))
cat("Tissue biopsy location breakdown:\n")
print(table(htx_meta$tissue_location))
cat("\n")

## ---- 4. Load & Inspect Host Count Table -------------------------------------
cat("--- Reading host transcriptomics count table (fread) ---\n")
t0 <- proc.time()
htx_raw <- fread(htx_dest, header = TRUE, sep = "\t", data.table = TRUE)
t1 <- proc.time()
cat(sprintf("Loaded in %.1f seconds. Raw table: %d genes × %d columns\n",
            (t1 - t0)[["elapsed"]], nrow(htx_raw), ncol(htx_raw)))

gene_col_name <- colnames(htx_raw)[1]
raw_sample_cols <- colnames(htx_raw)[-1]
cat("Gene identifier column name:", gene_col_name, "\n")
cat("Sample columns count:", length(raw_sample_cols), "\n")
cat("First 5 sample columns:\n")
print(head(raw_sample_cols, 5))

# Identify overlap with metadata Sample_IDs
overlap_samples <- intersect(htx_meta$Sample_ID, raw_sample_cols)
cat(sprintf("\nOverlap with HTX metadata Sample_IDs: %d / %d samples\n",
            length(overlap_samples), nrow(htx_meta)))

if (length(overlap_samples) == 0) {
  # Some files have GSM or modified IDs; check matching
  cat("Checking alternative column match...\n")
  clean_cols <- sub("\\.counts$|\\.bam$", "", raw_sample_cols)
  overlap_samples <- intersect(htx_meta$Sample_ID, clean_cols)
}

if (length(overlap_samples) < 50) {
  warning("Sample matching with metadata is low. Proceeding with all available count columns.")
  use_samples <- raw_sample_cols
} else {
  use_samples <- overlap_samples
}

## ---- 5. Data Cleanliness & QC Verification ----------------------------------
cat("\n=============================================================\n")
cat("DATA CLEANLINESS & INTEGRITY CHECK\n")
cat("=============================================================\n")

# Extract count matrix
count_dt <- htx_raw[, c(gene_col_name, use_samples), with = FALSE]
gene_ids <- count_dt[[1]]
count_mat <- as.matrix(count_dt[, -1, with = FALSE])
rownames(count_mat) <- gene_ids

# Ensure numeric matrix
mode(count_mat) <- "numeric"

# 1. Missing values check
na_count <- sum(is.na(count_mat))
nan_count <- sum(is.nan(count_mat))
inf_count <- sum(is.infinite(count_mat))
neg_count <- sum(count_mat < 0, na.rm = TRUE)

cat("1. Missing / NA values:     ", na_count, "\n")
cat("2. NaN values:              ", nan_count, "\n")
cat("3. Infinite values:         ", inf_count, "\n")
cat("4. Negative count values:   ", neg_count, "\n")

# Replace any stray NA with 0 if present
if (na_count > 0) {
  count_mat[is.na(count_mat)] <- 0
  cat("  -> Imputed NAs as 0.\n")
}

# 2. Library sizes (sequencing depth per sample)
lib_sizes <- colSums(count_mat)
cat(sprintf("5. Library size per sample:\n"))
cat(sprintf("     Min depth:    %s reads\n", comma(min(lib_sizes))))
cat(sprintf("     Median depth: %s reads\n", comma(median(lib_sizes))))
cat(sprintf("     Max depth:    %s reads\n", comma(max(lib_sizes))))

zero_depth_samples <- sum(lib_sizes == 0)
cat("6. Zero-depth samples:      ", zero_depth_samples, "\n")
if (zero_depth_samples > 0) {
  keep_cols <- names(lib_sizes)[lib_sizes > 0]
  count_mat <- count_mat[, keep_cols, drop = FALSE]
  lib_sizes <- lib_sizes[keep_cols]
}

# 3. Gene expression filtering (remove unexpressed / extremely low genes)
# Keep genes with at least 10 reads in at least 10% of samples
min_samples <- max(3, floor(0.10 * ncol(count_mat)))
gene_expressed <- rowSums(count_mat >= 10) >= min_samples
cat(sprintf("7. Total host genes evaluated: %d\n", nrow(count_mat)))
cat(sprintf("8. Expressed genes kept (>=10 reads in >=%d samples): %d (%.1f%%)\n",
            min_samples, sum(gene_expressed), 100 * mean(gene_expressed)))

count_filtered <- count_mat[gene_expressed, , drop = FALSE]

## ---- 6. Normalization: log2-CPM Transformation ------------------------------
cat("\n--- Normalising: Computing log2-CPM ---\n")
# log2-CPM = log2((count / library_size) * 1e6 + 1)
cpm_mat <- sweep(count_filtered, 2, lib_sizes[colnames(count_filtered)], "/") * 1e6
log2_cpm <- log2(cpm_mat + 1)

# Align metadata with sample columns in log2_cpm
aligned_meta <- htx_meta %>%
  filter(Sample_ID %in% colnames(log2_cpm)) %>%
  arrange(match(Sample_ID, colnames(log2_cpm)))

# If metadata is incomplete for some samples, construct dummy labels
if (nrow(aligned_meta) < ncol(log2_cpm)) {
  missing_in_meta <- setdiff(colnames(log2_cpm), aligned_meta$Sample_ID)
  dummy_rows <- data.frame(
    Sample_ID       = missing_in_meta,
    Participant_ID  = paste0("Unknown_", seq_along(missing_in_meta)),
    diagnosis_std   = "Unknown",
    tissue_location = "Biopsy",
    consent_age     = NA_real_,
    stringsAsFactors = FALSE
  )
  aligned_meta <- bind_rows(aligned_meta, dummy_rows) %>%
    arrange(match(Sample_ID, colnames(log2_cpm)))
}

# Save processed log2-CPM matrix
saveRDS(list(
  log2_cpm     = log2_cpm,
  count_mat    = count_filtered,
  metadata     = aligned_meta,
  lib_sizes    = lib_sizes[colnames(log2_cpm)]
), file.path(data_proc, "host_tx_log2cpm_filtered.rds"))
cat("Saved processed data: data/processed/host_tx_log2cpm_filtered.rds\n\n")

## ---- 7. Visualisation 1: Library Size Distribution --------------------------
cat("--- Visualisation 1: Sequencing depth distribution ---\n")
qc_df <- data.frame(
  Sample_ID       = names(lib_sizes[colnames(log2_cpm)]),
  Library_Size    = lib_sizes[colnames(log2_cpm)],
  stringsAsFactors = FALSE
) %>%
  left_join(aligned_meta, by = "Sample_ID")

p1_qc <- ggplot(qc_df, aes(x = diagnosis_std, y = Library_Size, fill = diagnosis_std)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "white") +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.8) +
  scale_y_continuous(labels = comma, trans = "log10") +
  scale_fill_manual(values = diag_colors, na.value = "grey60", name = "Diagnosis") +
  facet_wrap(~ tissue_location, scales = "free_y") +
  labs(
    title = "Host Transcriptomics: Sequencing Depth by Diagnosis & Tissue",
    subtitle = paste("N =", nrow(qc_df), "mucosal biopsy samples"),
    x = "Diagnosis",
    y = expression("Total Mapped Reads (" * log[10] * " scale)")
  ) +
  theme_htx +
  theme(legend.position = "none")

# Show in console / plot window
tryCatch(print(p1_qc), error = function(e) invisible(NULL))

# Save image (PNG, 300 DPI)
img1_path <- file.path(figs_dir, "htx_qc_library_size.png")
ggsave(img1_path, plot = p1_qc, width = 8, height = 5, dpi = 300)
cat("Saved image:", img1_path, "\n\n")

## ---- 8. Visualisation 2: Principal Component Analysis (PCA) -----------------
cat("--- Visualisation 2: Principal Component Analysis (PCA) ---\n")
# Select top 1,000 most variable genes by variance across samples
gene_vars <- apply(log2_cpm, 1, var)
top1000_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(1000, length(gene_vars))]

pca_fit <- prcomp(t(log2_cpm[top1000_genes, ]), center = TRUE, scale. = TRUE)
pca_var <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2) * 100

pca_df <- as.data.frame(pca_fit$x[, 1:4])
pca_df$Sample_ID <- rownames(pca_df)
pca_df <- left_join(pca_df, aligned_meta, by = "Sample_ID")

p2_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = diagnosis_std, shape = tissue_location)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = diag_colors, na.value = "grey60", name = "Diagnosis") +
  scale_shape_manual(values = c("Ileum" = 17, "Rectum" = 16, "Colon" = 15, "Other" = 3, "Biopsy" = 18),
                     name = "Tissue Location") +
  labs(
    title = "Host Mucosal Transcriptome: PCA (Top 1,000 Variable Genes)",
    subtitle = "Separation across intestinal tissue segments and IBD diagnosis",
    x = sprintf("PC1 (%.1f%% variance)", pca_var[1]),
    y = sprintf("PC2 (%.1f%% variance)", pca_var[2])
  ) +
  theme_htx +
  theme(legend.box = "horizontal")

# Show in console / plot window
tryCatch(print(p2_pca), error = function(e) invisible(NULL))

# Save image (PNG, 300 DPI)
img2_path <- file.path(figs_dir, "htx_pca_diagnosis_tissue.png")
ggsave(img2_path, plot = p2_pca, width = 8, height = 6, dpi = 300)
cat("Saved image:", img2_path, "\n\n")

## ---- 9. Visualisation 3: Hierarchical Clustering Heatmap --------------------
cat("--- Visualisation 3: Heatmap of Top 50 Variable Host Genes ---\n")
top50_genes <- names(sort(gene_vars, decreasing = TRUE))[1:50]
heatmap_mat <- log2_cpm[top50_genes, ]

# Row z-score transformation
heatmap_z <- t(scale(t(heatmap_mat)))
heatmap_z[heatmap_z > 3]  <- 3
heatmap_z[heatmap_z < -3] <- -3

# Annotation columns
anno_col <- data.frame(
  Diagnosis = aligned_meta$diagnosis_std,
  Tissue    = aligned_meta$tissue_location,
  row.names = aligned_meta$Sample_ID
)

anno_colors <- list(
  Diagnosis = diag_colors,
  Tissue    = c("Ileum" = "#E67E22", "Rectum" = "#9B59B6", "Colon" = "#1ABC9C", "Other" = "#95A5A6", "Biopsy" = "#BDC3C7")
)

# Show heatmap in console / plot window
tryCatch({
  pheatmap(
    heatmap_z,
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = "Hierarchical Clustering: Top 50 Variable Host Genes (Z-score)",
    color = colorRampPalette(c("#2980B9", "#ECF0F1", "#C0392B"))(100)
  )
}, error = function(e) invisible(NULL))

# Save image (PNG, 300 DPI)
img3_path <- file.path(figs_dir, "htx_heatmap_top50_genes.png")
png(img3_path, width = 2400, height = 2400, res = 300)
pheatmap(
  heatmap_z,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  main = "Hierarchical Clustering: Top 50 Variable Host Genes (Z-score)",
  color = colorRampPalette(c("#2980B9", "#ECF0F1", "#C0392B"))(100)
)
dev.off()
cat("Saved image:", img3_path, "\n\n")

## ---- 10. Visualisation 4: Bar Plots of Relative Host Gene Expression --------
cat("--- Visualisation 4: Bar plots of top variable genes across diagnosis ---\n")
# Select top 6 most variable host genes
top6_genes <- names(sort(gene_vars, decreasing = TRUE))[1:6]

bar_df <- as.data.frame(t(log2_cpm[top6_genes, ]))
bar_df$Sample_ID <- rownames(bar_df)
bar_df_long <- bar_df %>%
  pivot_longer(cols = all_of(top6_genes), names_to = "Gene", values_to = "log2_CPM") %>%
  left_join(aligned_meta, by = "Sample_ID") %>%
  filter(diagnosis_std %in% c("CD", "UC", "nonIBD"))

# Compute mean and standard error per diagnosis
summary_df <- bar_df_long %>%
  group_by(Gene, diagnosis_std) %>%
  summarise(
    mean_expr = mean(log2_CPM, na.rm = TRUE),
    se_expr   = sd(log2_CPM, na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    .groups   = "drop"
  )

p4_bars <- ggplot(summary_df, aes(x = diagnosis_std, y = mean_expr, fill = diagnosis_std)) +
  geom_col(width = 0.6, alpha = 0.85, colour = "black") +
  geom_errorbar(aes(ymin = pmax(0, mean_expr - se_expr), ymax = mean_expr + se_expr),
                width = 0.2, linewidth = 0.7) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = diag_colors, name = "Diagnosis") +
  labs(
    title = "Relative Host Gene Expression Across IBD Subtypes",
    subtitle = "Mean log2-CPM (+/- Standard Error) for Top Variable Host Mucosal Genes",
    x = "Diagnosis",
    y = expression("Relative Expression (Mean " * log[2] * "-CPM)")
  ) +
  theme_htx +
  theme(legend.position = "bottom")

# Show in console / plot window
tryCatch(print(p4_bars), error = function(e) invisible(NULL))

# Save image (PNG, 300 DPI)
img4_path <- file.path(figs_dir, "htx_barplot_top_variable_genes.png")
ggsave(img4_path, plot = p4_bars, width = 9, height = 6, dpi = 300)
cat("Saved image:", img4_path, "\n\n")

## Clean up default Rplots.pdf if created by headless R batch mode during print()
if (file.exists("Rplots.pdf")) {
  invisible(file.remove("Rplots.pdf"))
}

## ---- 11. Save QC Report -----------------------------------------------------
report_file <- file.path(results_qc, "host_transcriptomics_QC_report.txt")
sink(report_file)
cat("=============================================================\n")
cat("Host Transcriptomics (HTX) QC & Exploration Report\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=============================================================\n\n")
cat("1. Raw Table Dimensions: ", nrow(htx_raw), "genes ×", ncol(htx_raw) - 1, "samples\n")
cat("2. Analyzed Samples:     ", ncol(count_mat), "\n")
cat("3. Expressed Genes Kept: ", nrow(count_filtered), "(out of", nrow(count_mat), "total)\n")
cat("4. Missing / NAs:        ", na_count, "\n")
cat("5. Library Size Depth:   Median =", comma(median(lib_sizes)), "reads (Range:", comma(min(lib_sizes)), "-", comma(max(lib_sizes)), ")\n\n")
cat("6. Sample Breakdown by Diagnosis:\n")
print(table(aligned_meta$diagnosis_std))
cat("\n7. Sample Breakdown by Tissue Location:\n")
print(table(aligned_meta$tissue_location))
cat("\n8. Top 10 Most Variable Host Genes:\n")
print(head(sort(gene_vars, decreasing = TRUE), 10))
cat("\n9. Visualisation Outputs (PNG Images):\n")
cat("   - figures/htx_qc_library_size.png\n")
cat("   - figures/htx_pca_diagnosis_tissue.png\n")
cat("   - figures/htx_heatmap_top50_genes.png\n")
cat("   - figures/htx_barplot_top_variable_genes.png\n")
cat("=============================================================\n")
sink()

cat("QC report saved to:", report_file, "\n")
cat("=============================================================\n")
cat("SCRIPT h1 COMPLETE\n")
cat("=============================================================\n")
cat("Session ended:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
