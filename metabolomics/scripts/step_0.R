#Finding the table of the metadata

f <- "data/raw/ST000923_AN001513_Results.txt"

cat("size:", round(file.size(f)/1e6, 1), "MB\n")

first <- readLines(f, n = 40)
for (i in seq_along(first)) {
  cat(sprintf("%2d| %s\n", i, substr(first[i], 1, 150)))
}


#Fetching metabolomics files from ibdmdb
setwd("/home/ritwikan/metabolomics")
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

options(timeout = 3600)

base <- "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MBX/"

files <- c(
  "HMP2_metabolomics_w_metadata.biom.gz",
  "HMP2_metabolomics.biom",
  "C18n_Metabolites_ID_AfterPublication_2021-05-17.xlsx"
)

for (fn in files) {
  dest <- file.path("data/raw", fn)
  if (file.exists(dest)) { cat(fn, "- already have\n"); next }
  cat("\ndownloading", fn, "\n")
  try(download.file(paste0(base, fn), dest, mode = "wb"))
}

# what landed
data.frame(
  file = basename(list.files("data/raw", full.names = TRUE)),
  MB   = round(file.size(list.files("data/raw", full.names = TRUE))/1e6, 2)
)

#installign requisite packages
install.packages("BiocManager")
BiocManager::install("biomformat")
install.packages(c("readxl", "R.utils"))


library(biomformat)

R.utils::gunzip("data/raw/HMP2_metabolomics_w_metadata.biom.gz",
                destname = "data/raw/HMP2_metabolomics_w_metadata.biom",
                remove = FALSE)

bi <- read_biom("data/raw/HMP2_metabolomics_w_metadata.biom")

smeta <- sample_metadata(bi)
colnames(smeta)
head(smeta, 3)
head(rownames(smeta), 5)

# what did we actually get?
dim(biom_data(bi))
head(colnames(biom_data(bi)), 5)      # sample IDs
head(rownames(biom_data(bi)), 5)      # feature IDs

# is there feature-level metadata at least?
fm <- observation_metadata(bi)
if (is.null(fm)) cat("no feature metadata either\n") else {
  cat("feature metadata columns:\n"); print(colnames(fm)); print(head(fm, 3))
}

# look at the raw structure
str(bi, max.level = 2)

# 1. Is there feature metadata? (HMDB IDs, QC CV, m/z, RT)
fm <- observation_metadata(bi)
if (is.null(fm)) {
  cat("none\n")
} else {
  print(colnames(fm))
  print(head(fm, 3))
}

# 2. Method breakdown from the feature prefixes
feat <- rownames(biom_data(bi))
method <- sub("_.*$", "", feat)
print(table(method))

# 3. Pull the matrix out and save it
mat <- as.matrix(biom_data(bi))
cat("dim:", dim(mat), "\n")
cat("fraction zero/missing:", round(mean(mat == 0 | is.na(mat)), 3), "\n")

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
saveRDS(mat, "data/derived/mbx_matrix.rds")

# ============================================================
# 02_feature_info.R
# Flatten the feature metadata and answer CENSUS 6
# ============================================================
#
# observation_metadata() gives a LIST (one named vector per
# feature), not a data frame. We turn it into a table.
# ============================================================

setwd("/home/ritwikan/metabolomics")

library(biomformat)

bi  <- read_biom("data/raw/HMP2_metabolomics_w_metadata.biom")
fm  <- observation_metadata(bi)
mat <- readRDS("data/derived/mbx_matrix.rds")


# ============================================================
# 1. FLATTEN THE LIST INTO A DATA FRAME
# ============================================================
#
# Pull each field by NAME, so a feature missing a field gives
# NA rather than silently shifting the columns across.

grab <- function(x, field) {
  v <- x[[field]]
  if (is.null(v)) NA else as.character(v)
}

feature_info <- data.frame(
  feature_id = names(fm),
  RT         = sapply(fm, grab, "RT"),
  mz         = sapply(fm, grab, "m/z"),
  metabolite = sapply(fm, grab, "Metabolite"),
  hmdb       = sapply(fm, grab, "HMDB (*Representative ID)"),
  qc_cv      = sapply(fm, grab, "Pooled QC sample CV"),
  method     = sapply(fm, grab, "Method"),
  stringsAsFactors = FALSE
)

# numeric columns came through as text
feature_info$RT    <- as.numeric(feature_info$RT)
feature_info$mz    <- as.numeric(feature_info$mz)
feature_info$qc_cv <- as.numeric(feature_info$qc_cv)

rownames(feature_info) <- NULL

cat("=== FEATURE TABLE ===\n")
cat("rows:", nrow(feature_info), "\n")
print(head(feature_info))

# sanity: same order as the matrix?
cat("\nfeature order matches matrix:",
    identical(feature_info$feature_id, rownames(mat)), "\n")


# ============================================================
# 2. CENSUS 6 - WHAT FRACTION IS ANNOTATED?
# ============================================================
#
# A feature is usable for the single-cell bridge only if it has
# a real HMDB ID. Blanks, NA, and "redundant ion" do not count.

has_hmdb <- !is.na(feature_info$hmdb) &
  feature_info$hmdb != "" &
  !grepl("redundant", feature_info$hmdb, ignore.case = TRUE)

has_name <- !is.na(feature_info$metabolite) &
  feature_info$metabolite != ""

cat("\n\n=== CENSUS 6: ANNOTATION ===\n")
cat("total features:          ", nrow(feature_info), "\n")
cat("with a metabolite name:  ", sum(has_name),
    sprintf("(%.2f%%)", 100 * mean(has_name)), "\n")
cat("with an HMDB ID:         ", sum(has_hmdb),
    sprintf("(%.2f%%)", 100 * mean(has_hmdb)), "\n")

cat("\n--- by method ---\n")
by_method <- data.frame(
  method     = feature_info$method,
  annotated  = has_hmdb
)
tab <- aggregate(annotated ~ method, data = by_method,
                 FUN = function(x) c(n = length(x),
                                     annotated = sum(x),
                                     pct = round(100 * mean(x), 2)))
print(do.call(data.frame, tab))

feature_info$has_hmdb <- has_hmdb


# ============================================================
# 3. THE QC CV PICTURE
# ============================================================
#
# This drives Person A's filter threshold. Look before choosing.

cat("\n\n=== POOLED QC CV ===\n")
cat("features with a CV value:", sum(!is.na(feature_info$qc_cv)), "\n")
print(summary(feature_info$qc_cv))

for (thr in c(0.15, 0.20, 0.25, 0.30)) {
  keep <- !is.na(feature_info$qc_cv) & feature_info$qc_cv <= thr
  cat(sprintf("CV <= %.2f : %6d features kept (%.1f%%)\n",
              thr, sum(keep), 100 * mean(keep)))
}

dir.create("output", showWarnings = FALSE)
png("output/qc_cv_distribution.png", width = 700, height = 500)
hist(feature_info$qc_cv[feature_info$qc_cv < 1],
     breaks = 60,
     main = "Pooled QC sample CV",
     xlab = "CV")
abline(v = 0.25, col = "red", lwd = 2)
dev.off()


# ============================================================
# 4. MISSINGNESS PER FEATURE
# ============================================================
#
# 38.8% of the matrix is zero/missing overall. That is
# left-censored (below detection), not random.

detect_frac <- rowMeans(mat > 0 & !is.na(mat))
feature_info$detect_frac <- detect_frac

cat("\n\n=== DETECTION FREQUENCY ===\n")
print(summary(detect_frac))

for (thr in c(0.25, 0.50, 0.75)) {
  cat(sprintf("detected in >= %.0f%% of samples : %6d features\n",
              100 * thr, sum(detect_frac >= thr)))
}

cat("\nAnnotated AND well-detected (HMDB + >=50%% of samples):",
    sum(has_hmdb & detect_frac >= 0.5), "\n")
cat("^ this is roughly what the single-cell bridge will work with\n")


# ============================================================
# 5. SAVE
# ============================================================

saveRDS(feature_info, "data/derived/feature_info.rds")
write.csv(feature_info, "output/feature_info.csv", row.names = FALSE)

cat("\n=== DONE ===\n")
cat("Put the CENSUS 6 numbers in the week 1 page.\n")

fi <- readRDS("data/derived/feature_info.rds")

cat("CV among ANNOTATED features:\n")
print(summary(fi$qc_cv[fi$has_hmdb]))

cat("\nCV among UNANNOTATED:\n")
print(summary(fi$qc_cv[!fi$has_hmdb]))

cat("\nannotated features surviving each CV threshold:\n")
for (thr in c(0.15, 0.20, 0.25, 0.30)) {
  k <- fi$has_hmdb & !is.na(fi$qc_cv) & fi$qc_cv <= thr
  cat(sprintf("  CV <= %.2f : %d of 522\n", thr, sum(k)))
}

fi <- readRDS("data/derived/feature_info.rds")

cat("CV among ANNOTATED features:\n")
print(summary(fi$qc_cv[fi$has_hmdb]))

cat("\nCV among UNANNOTATED:\n")
print(summary(fi$qc_cv[!fi$has_hmdb]))

cat("\nannotated features surviving each CV threshold:\n")
for (thr in c(0.15, 0.20, 0.25, 0.30)) {
  k <- fi$has_hmdb & !is.na(fi$qc_cv) & fi$qc_cv <= thr
  cat(sprintf("  CV <= %.2f : %d of 522\n", thr, sum(k)))
}


fi <- readRDS("data/derived/feature_info.rds")

cat("CV among ANNOTATED features:\n")
print(summary(fi$qc_cv[fi$has_hmdb]))

cat("\nCV among UNANNOTATED:\n")
print(summary(fi$qc_cv[!fi$has_hmdb]))

cat("\nannotated features surviving each CV threshold:\n")
for (thr in c(0.15, 0.20, 0.25, 0.30)) {
  k <- fi$has_hmdb & !is.na(fi$qc_cv) & fi$qc_cv <= thr
  cat(sprintf("  CV <= %.2f : %d of 522\n", thr, sum(k)))
}

file.edit("episode_definition.md")


setwd("/home/ritwikan/metabolomics")
options(timeout = 1800)

download.file("https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MGX/2018-05-04/taxonomic_profiles.tsv.gz", "data/raw/taxonomic_profiles.tsv.gz", mode = "wb")

cat("size:", round(file.size("data/raw/taxonomic_profiles.tsv.gz")/1e6, 1), "MB\n")

tax <- read.delim(gzfile("data/raw/taxonomic_profiles.tsv.gz"),
                  nrows = 3, check.names = FALSE)
cat("columns:", ncol(tax), "\n")
cat("first taxon:", substr(tax[[1]][1], 1, 70), "\n")
cat("first 3 samples:", paste(colnames(tax)[2:4], collapse = ", "), "\n")


#downloading hmp2 metadata
setwd("/home/ritwikan/metabolomics")
options(timeout = 1800)

download.file(
  "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/metadata/hmp2_metadata_2018-08-20.csv",
  "data/raw/hmp2_metadata.csv", mode = "wb")

cat("size:", round(file.size("data/raw/hmp2_metadata.csv")/1e6, 2), "MB\n")

meta <- read.csv("data/raw/hmp2_metadata.csv", check.names = FALSE)

cat("rows:", nrow(meta), " cols:", ncol(meta), "\n\n")
print(colnames(meta))

cat("\n--- data types ---\n")
print(table(meta$data_type))

mbx_ids <- colnames(readRDS("data/derived/mbx_matrix.rds"))
tax_ids <- colnames(read.delim(gzfile("data/raw/taxonomic_profiles.tsv.gz"),
                               nrows = 1, check.names = FALSE))[-1]

cat("\nmetabolomics matched:", sum(mbx_ids %in% meta[["External ID"]]), "of 546\n")
cat("taxonomy matched:", sum(tax_ids %in% meta[["External ID"]]), "of 1638\n")
cat("subjects with metabolomics:",
    length(unique(meta[["Participant ID"]][meta[["External ID"]] %in% mbx_ids])), "\n")

##########################################################################################################
# ============================================================
# 03_dysbiosis_score.R
# Dysbiosis score, Lloyd-Price et al. 2019
# ============================================================
#
# RESTART R BEFORE RUNNING THIS. A half-run session can leave a
# corrupted `abund` in memory and the results will be nonsense.
#
# VALIDATION TARGETS (from the paper):
#   17.1% dysbiotic overall,  24.3% in CD,  11.6% in UC
# ============================================================

setwd("/home/ritwikan/metabolomics")

library(vegan)

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)


# ============================================================
# 1. LOAD
# ============================================================

meta <- read.csv("data/raw/hmp2_metadata.csv", check.names = FALSE,
                 stringsAsFactors = FALSE)

tax_raw <- read.delim(gzfile("data/raw/taxonomic_profiles.tsv.gz"),
                      check.names = FALSE, stringsAsFactors = FALSE)

cat("taxonomy table:", nrow(tax_raw), "rows x", ncol(tax_raw), "cols\n")


# ============================================================
# 2. SPECIES ROWS ONLY
# ============================================================

taxon_names <- tax_raw[[1]]

is_species <- grepl("s__", taxon_names, fixed = TRUE) &
  !grepl("t__", taxon_names, fixed = TRUE)

cat("all rows:", nrow(tax_raw), "-> species rows:", sum(is_species), "\n")

abund <- as.matrix(tax_raw[is_species, -1])
rownames(abund) <- taxon_names[is_species]

# vegan wants rows = samples
abund <- t(abund)
storage.mode(abund) <- "numeric"

cat("matrix:", nrow(abund), "samples x", ncol(abund), "species\n")


# ============================================================
# 3. CLEAN  -  ORDER MATTERS HERE
# ============================================================
#
# Do these BEFORE dividing. A single zero-sum row produces NaN,
# which then turns colSums() into NA, and subsetting a matrix
# with an NA logical fills whole columns with NA. One bad row
# destroys the entire matrix if you normalise first.

# 3a. missing cells are genuine absences
n_na <- sum(is.na(abund))
cat("NA cells in source:", n_na, "\n")
abund[is.na(abund)] <- 0

# 3b. drop samples with nothing detected
tot <- rowSums(abund)
cat("samples with zero total:", sum(tot == 0), "\n")
abund <- abund[tot > 0, , drop = FALSE]

# 3c. NOW it is safe to normalise
cat("row sums before (first 5):", round(rowSums(abund)[1:5], 3), "\n")
abund <- abund / rowSums(abund)

# 3d. drop species absent everywhere
keep <- colSums(abund) > 0
keep[is.na(keep)] <- FALSE          # belt and braces
abund <- abund[, keep, drop = FALSE]

# 3e. fail loudly if anything is still wrong
stopifnot(!any(is.na(abund)))
stopifnot(nrow(abund) > 1000, ncol(abund) > 100)

cat("after cleaning:", nrow(abund), "samples x", ncol(abund), "species\n")


# ============================================================
# 4. ATTACH DIAGNOSIS AND SUBJECT
# ============================================================
#
# External ID is unique per sample PER ASSAY, so look up only
# the metagenomics rows.

meta_mgx <- meta[meta$data_type == "metagenomics", ]

idx <- match(rownames(abund), meta_mgx[["External ID"]])

si <- data.frame(
  sample_id  = rownames(abund),
  subject_id = meta_mgx[["Participant ID"]][idx],
  diagnosis  = meta_mgx[["diagnosis"]][idx],
  week       = meta_mgx[["week_num"]][idx],
  stringsAsFactors = FALSE
)
si$day <- si$week * 7

cat("\nunmatched samples:", sum(is.na(si$subject_id)), "\n")

ok <- !is.na(si$subject_id)
abund <- abund[ok, , drop = FALSE]
si     <- si[ok, ]

cat("samples kept:", nrow(abund), "\n")
cat("subjects:", length(unique(si$subject_id)), "\n\n")
print(table(si$diagnosis, useNA = "ifany"))

stopifnot(nrow(si) == nrow(abund))


# ============================================================
# 5. BRAY-CURTIS
# ============================================================

cat("\ncomputing Bray-Curtis (takes a minute)...\n")

bc <- as.matrix(vegdist(abund, method = "bray"))

cat("distance matrix:", nrow(bc), "x", ncol(bc), "\n")


# ============================================================
# 6. SCORE EVERY SAMPLE
# ============================================================
#
# median distance to non-IBD samples FROM OTHER SUBJECTS

is_nonibd <- si$diagnosis == "nonIBD"
cat("non-IBD reference samples:", sum(is_nonibd), "\n")

score <- rep(NA_real_, nrow(si))

for (i in seq_len(nrow(si))) {
  use <- is_nonibd & (si$subject_id != si$subject_id[i])
  score[i] <- median(bc[i, use])
}

si$dysbiosis_score <- score

cat("\nscore summary:\n")
print(summary(score))


# ============================================================
# 7. THRESHOLD
# ============================================================

threshold <- quantile(score[is_nonibd], 0.90, na.rm = TRUE)

cat("\n*** threshold (90th pct of non-IBD) =",
    round(threshold, 4), "***\n")

si$dysbiotic <- si$dysbiosis_score > threshold


# ============================================================
# 8. VALIDATION
# ============================================================

cat("\n=== VALIDATION vs Lloyd-Price et al. 2019 ===\n\n")

tab <- data.frame(
  diagnosis   = names(table(si$diagnosis)),
  n_samples   = as.integer(table(si$diagnosis)),
  n_dysbiotic = as.integer(tapply(si$dysbiotic, si$diagnosis, sum)),
  fraction    = round(as.numeric(tapply(si$dysbiotic, si$diagnosis, mean)), 3)
)
print(tab)

cat("\noverall dysbiotic:", sum(si$dysbiotic), "of", nrow(si),
    sprintf("(%.1f%%)", 100 * mean(si$dysbiotic)), "\n")

cat("\nPAPER REPORTS:\n")
cat("  overall  17.1%  (272 samples)\n")
cat("  CD       24.3%  (n=178)\n")
cat("  UC       11.6%  (n=51)\n")
cat("\nCHECK: is CD roughly double UC? If not, STOP and debug.\n")


# ============================================================
# 9. SAVE
# ============================================================

saveRDS(si,        "data/derived/dysbiosis.rds")
saveRDS(threshold, "data/derived/dysbiosis_threshold.rds")
write.csv(si, "output/dysbiosis_scores.csv", row.names = FALSE)

png("output/dysbiosis_by_diagnosis.png", width = 700, height = 500)
boxplot(dysbiosis_score ~ diagnosis, data = si,
        main = "Dysbiosis score by diagnosis",
        ylab = "median Bray-Curtis to non-IBD")
abline(h = threshold, col = "red", lwd = 2, lty = 2)
dev.off()

cat("\n=== DONE ===\n")



# ============================================================
# 04_episodes_and_census.R
# Census numbers 1-5, and the episodes they need
# ============================================================
#
# RESTART R BEFORE RUNNING. Needs 03_dysbiosis_score.R done.
#
# Rules below MUST match episode_definition.md.
# ============================================================

setwd("/home/ritwikan/metabolomics")

dir.create("output", showWarnings = FALSE)

# ---- frozen constants (match episode_definition.md) --------
MIN_SEPARATION_DAYS  <- 30
MIN_PAIRS_PER_BIN    <- 30
MIN_SUBJECTS_PER_BIN <- 15


# ============================================================
# 0. LOAD
# ============================================================

meta <- read.csv("data/raw/hmp2_metadata.csv", check.names = FALSE,
                 stringsAsFactors = FALSE)

dys <- readRDS("data/derived/dysbiosis.rds")   # metagenomic samples
mat <- readRDS("data/derived/mbx_matrix.rds")  # metabolomics

mbx_ids <- colnames(mat)


# ---- metabolomic samples, with subject and day -------------

meta_mbx <- meta[meta$data_type == "metabolomics", ]
idx <- match(mbx_ids, meta_mbx[["External ID"]])

mbx <- data.frame(
  sample_id  = mbx_ids,
  subject_id = meta_mbx[["Participant ID"]][idx],
  diagnosis  = meta_mbx[["diagnosis"]][idx],
  week       = meta_mbx[["week_num"]][idx],
  stringsAsFactors = FALSE
)
mbx$day <- mbx$week * 7

cat("metabolomic samples:", nrow(mbx), "\n")
cat("  missing subject:", sum(is.na(mbx$subject_id)), "\n")
cat("  missing week:   ", sum(is.na(mbx$day)), "\n")

mbx <- mbx[!is.na(mbx$subject_id) & !is.na(mbx$day), ]
mbx <- mbx[order(mbx$subject_id, mbx$day), ]

cat("usable:", nrow(mbx), "samples from",
    length(unique(mbx$subject_id)), "subjects\n")


# ============================================================
# CENSUS 1 - TIMEPOINTS PER SUBJECT
# ============================================================

n_per <- table(mbx$subject_id)

cat("\n\n=== CENSUS 1: metabolomic timepoints per subject ===\n")
cat("subjects total:", length(n_per), "\n")
cat("  >=3 timepoints:", sum(n_per >= 3), "\n")
cat("  >=4 timepoints:", sum(n_per >= 4), "\n")
cat("  >=5 timepoints:", sum(n_per >= 5), "\n")

cat("\ndistribution:\n")
print(table(as.integer(n_per)))

# by diagnosis
dx <- mbx$diagnosis[match(names(n_per), mbx$subject_id)]
cat("\nsubjects with >=3, by diagnosis:\n")
print(table(dx[n_per >= 3]))


# ============================================================
# CENSUS 2 - GAP BETWEEN CONSECUTIVE SAMPLES
# ============================================================

gaps <- c()
for (s in unique(mbx$subject_id)) {
  d <- sort(mbx$day[mbx$subject_id == s])
  if (length(d) >= 2) gaps <- c(gaps, diff(d))
}

cat("\n\n=== CENSUS 2: days between consecutive samples ===\n")
cat("number of gaps:", length(gaps), "\n")
cat("median:", median(gaps), " IQR:",
    quantile(gaps, .25), "-", quantile(gaps, .75), "\n\n")
print(summary(gaps))

png("output/gap_distribution.png", width = 700, height = 500)
hist(gaps, breaks = 40,
     main = "Days between consecutive metabolomic samples",
     xlab = "days")
abline(v = median(gaps), col = "red", lwd = 2)
dev.off()
cat("\nsaved output/gap_distribution.png\n")
cat("LOOK AT IT before setting the bin edges below.\n")


# ============================================================
# BUILD EPISODES
# ============================================================
#
# Episodes are defined on the METAGENOMIC samples (where the
# dysbiosis score lives), not the metabolomic ones.

d_all <- dys[!is.na(dys$day), ]
d_all <- d_all[order(d_all$subject_id, d_all$day), ]

cat("\n\n=== building episodes ===\n")
cat("samples with a dysbiosis score and a day:", nrow(d_all), "\n")

episodes <- data.frame()

for (s in unique(d_all$subject_id)) {
  
  d <- d_all[d_all$subject_id == s, ]
  d <- d[order(d$day), ]
  if (nrow(d) < 2) next
  
  already_dys <- d$dysbiotic[1]
  
  i <- 2
  while (i <= nrow(d)) {
    
    # ONSET: below threshold before, at/above now
    if (d$dysbiotic[i] && !d$dysbiotic[i - 1]) {
      
      onset_day <- d$day[i]
      
      # walk forward while still dysbiotic
      j <- i
      while (j < nrow(d) && d$dysbiotic[j + 1]) j <- j + 1
      last_dys_day <- d$day[j]
      
      # RESOLUTION: first sub-threshold sample after
      if (j < nrow(d)) {
        resolution_day <- d$day[j + 1]
        confirmed <- (j + 2 <= nrow(d)) && !d$dysbiotic[j + 2]
      } else {
        resolution_day <- NA
        confirmed <- FALSE
      }
      
      episodes <- rbind(episodes, data.frame(
        subject_id     = s,
        diagnosis      = d$diagnosis[1],
        onset_day      = onset_day,
        last_dys_day   = last_dys_day,
        resolution_day = resolution_day,
        confirmed_resolution = confirmed,
        already_dysbiotic_at_start = already_dys,
        stringsAsFactors = FALSE
      ))
      
      i <- j + 1
    } else {
      i <- i + 1
    }
  }
}

cat("raw episodes:", nrow(episodes), "\n")


# ---- merge episodes closer than MIN_SEPARATION_DAYS --------

episodes <- episodes[order(episodes$subject_id, episodes$onset_day), ]
keep <- rep(TRUE, nrow(episodes))

if (nrow(episodes) > 1) {
  for (i in 2:nrow(episodes)) {
    same <- episodes$subject_id[i] == episodes$subject_id[i - 1]
    prev <- episodes$resolution_day[i - 1]
    if (same && !is.na(prev) &&
        (episodes$onset_day[i] - prev) < MIN_SEPARATION_DAYS) {
      keep[i] <- FALSE
      episodes$resolution_day[i - 1] <- episodes$resolution_day[i]
      episodes$last_dys_day[i - 1]   <- episodes$last_dys_day[i]
    }
  }
}
episodes <- episodes[keep, ]

cat("after merging (<", MIN_SEPARATION_DAYS, "days apart):",
    nrow(episodes), "\n")
cat("unresolved (censored):", sum(is.na(episodes$resolution_day)), "\n")
cat("subjects with an episode:",
    length(unique(episodes$subject_id)), "\n\n")
print(table(episodes$diagnosis))

cat("\nPAPER REPORTS: 78 full periods + 9 censored\n")


# ============================================================
# CENSUS 3 - ONSETS USABLE FOR Q1
# ============================================================
#
# Need >=2 metabolomic samples BEFORE the onset.

usable3 <- c()
for (i in seq_len(nrow(episodes))) {
  n_before <- sum(mbx$subject_id == episodes$subject_id[i] &
                    mbx$day < episodes$onset_day[i])
  if (n_before >= 2) usable3 <- c(usable3, i)
}

cat("\n\n=== CENSUS 3: onsets usable for Q1 ===\n")
cat("total episodes:", nrow(episodes), "\n")
cat("with >=2 metabolomic samples before onset:", length(usable3), "\n")
if (length(usable3) > 0) {
  cat("from", length(unique(episodes$subject_id[usable3])), "subjects\n")
  print(table(episodes$diagnosis[usable3]))
}
cat("\nGATE: single digits -> Q1 is a case series, not a stat claim.\n")


# ============================================================
# CENSUS 4 - COMPLETE EPISODES FOR Q4
# ============================================================
#
# Need >=1 metabolomic sample in EACH of pre / during / post.

usable4 <- c()
for (i in seq_len(nrow(episodes))) {
  
  if (is.na(episodes$resolution_day[i])) next
  
  m <- mbx[mbx$subject_id == episodes$subject_id[i], ]
  
  n_pre    <- sum(m$day <  episodes$onset_day[i])
  n_during <- sum(m$day >= episodes$onset_day[i] &
                    m$day <= episodes$last_dys_day[i])
  n_post   <- sum(m$day >= episodes$resolution_day[i])
  
  if (n_pre >= 1 && n_during >= 1 && n_post >= 1) usable4 <- c(usable4, i)
}

cat("\n\n=== CENSUS 4: complete episodes for Q4 ===\n")
cat("with >=1 metabolomic sample in pre AND during AND post:",
    length(usable4), "\n")
if (length(usable4) > 0) {
  cat("from", length(unique(episodes$subject_id[usable4])), "subjects\n")
  print(table(episodes$diagnosis[usable4]))
}
cat("\nGATE: under ~15 -> demote Q4 to descriptive IN THE PROPOSAL NOW.\n")


# ============================================================
# CENSUS 5 - LEAD-TIME BINS
# ============================================================

pairs <- data.frame()
for (s in unique(mbx$subject_id)) {
  d <- sort(mbx$day[mbx$subject_id == s])
  if (length(d) < 2) next
  for (a in 1:(length(d) - 1)) {
    for (b in (a + 1):length(d)) {
      pairs <- rbind(pairs, data.frame(
        subject_id = s, from_day = d[a], to_day = d[b],
        lead_days = d[b] - d[a], stringsAsFactors = FALSE))
    }
  }
}

cat("\n\n=== CENSUS 5: lead-time bins ===\n")
cat("within-subject visit pairs:", nrow(pairs), "\n")

# >>> SET FROM output/gap_distribution.png <<<
# first bin narrow: the paper found metabolite dissimilarity
# plateaus at about two weeks, so that region matters most

bin_edges <- c(0, 45, 75, 120, 200, 400)

pairs$bin <- cut(pairs$lead_days, breaks = bin_edges, include.lowest = TRUE)

bs <- aggregate(cbind(n_pairs = lead_days) ~ bin, data = pairs, FUN = length)
bs$n_subjects <- sapply(bs$bin, function(b)
  length(unique(pairs$subject_id[pairs$bin == b])))
bs$usable <- bs$n_pairs >= 30 & bs$n_subjects >= 15
print(bs)

if (any(bs$usable))
  cat("\nlongest usable bin:", as.character(bs$bin[max(which(bs$usable))]), "\n")




# ============================================================
# SAVE
# ============================================================

saveRDS(episodes, "data/derived/episodes.rds")
saveRDS(mbx,      "data/derived/mbx_samples.rds")
saveRDS(pairs,    "data/derived/visit_pairs.rds")

write.csv(episodes, "output/episodes.csv", row.names = FALSE)
write.csv(bs,       "output/leadtime_bins.csv", row.names = FALSE)

cat("\n=== DONE ===\n")
cat("Copy census 1-5 into the week 1 page, then read the gates.\n")

