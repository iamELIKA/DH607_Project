# ============================================================
# 05_feature_matrix.R   -   Person A
# Build the analysis feature matrix
# ============================================================
#
# RESTART R BEFORE RUNNING.
#
# Needs: data/derived/feature_info.rds
#        data/derived/mbx_matrix.rds
#
# Produces: features_full_v1.parquet       (33,415 x 546)
#           features_annotated_v1.parquet  (405 x 546)
#           output/funnel.png, funnel.csv, funnel_by_method.csv
#
# Install once, then leave commented:
#   install.packages(c("arrow", "jsonlite"))
#   BiocManager::install("pcaMethods"); install.packages("imputeLCMD")
# ============================================================

setwd("/home/ritwikan/metabolomics")

library(arrow)

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)


# ============================================================
# THRESHOLDS  -  every one of these goes in the write-up
# ============================================================

CV_THRESHOLD        <- 0.25   # pooled-QC CV, primary
CV_SENSITIVITY      <- c(0.20, 0.30)
MIN_DETECT_FRAC     <- 0.50   # detected in >= this fraction of samples
DETECT_SENSITIVITY  <- c(0.25, 0.75)

IMPUTE_METHOD       <- "half_min"   # or "qrilc"
PSEUDOCOUNT         <- 1
DROP_MISSING_CV     <- TRUE   # features with no CV value: drop
DO_KEGG             <- TRUE

# values that appear in ID columns but are NOT identifiers
BAD_ID <- c("", "NA", "na", "-", "--", "N/A")


# ============================================================
# 1. LOAD
# ============================================================

fi  <- readRDS("data/derived/feature_info.rds")
mat <- readRDS("data/derived/mbx_matrix.rds")   # features x samples

stopifnot(identical(fi$feature_id, rownames(mat)))

cat("start:", nrow(mat), "features x", ncol(mat), "samples\n\n")
cat("features per method:\n")
print(table(fi$method, useNA = "ifany"))

# a few features have no Method - fall back to the name prefix
if (any(is.na(fi$method))) {
  pre <- sub("_.*$", "", fi$feature_id)
  map <- c(C18n = "C18-neg", HILn = "HILIC-neg",
           HILp = "HILIC-pos", C8p = "C8-pos")
  n_before <- sum(is.na(fi$method))
  fi$method[is.na(fi$method)] <- map[pre[is.na(fi$method)]]
  cat("\nfilled", n_before, "missing Method values from the name prefix\n")
}

methods <- sort(unique(fi$method))
cat("\nmethods:", paste(methods, collapse = ", "), "\n")


# ============================================================
# 2. PER-METHOD PROCESSING
# ============================================================
#
# Each method is filtered, imputed and transformed on its own.
# They are brought together only at step 3.

funnel    <- data.frame()
processed <- list()

for (m in methods) {
  
  cat("\n==================================================\n")
  cat("METHOD:", m, "\n")
  cat("==================================================\n")
  
  sel <- fi$method == m
  f_m <- fi[sel, ]
  x   <- mat[sel, , drop = FALSE]
  
  n0 <- nrow(x)
  cat("features in:", n0, "\n")
  
  # ---- 2a. pooled-QC CV filter -------------------------
  cv_ok <- !is.na(f_m$qc_cv) & f_m$qc_cv <= CV_THRESHOLD
  if (!DROP_MISSING_CV) cv_ok <- cv_ok | is.na(f_m$qc_cv)
  
  cat("  no CV value:", sum(is.na(f_m$qc_cv)),
      ifelse(DROP_MISSING_CV, "(dropped)", "(kept)"), "\n")
  
  x   <- x[cv_ok, , drop = FALSE]
  f_m <- f_m[cv_ok, ]
  n1 <- nrow(x)
  cat("  after CV <=", CV_THRESHOLD, ":", n1, "\n")
  
  for (thr in CV_SENSITIVITY)
    cat("    [CV <=", thr, "would keep",
        sum(!is.na(fi$qc_cv[sel]) & fi$qc_cv[sel] <= thr), "]\n")
  
  # ---- 2b. detection frequency -------------------------
  detected <- !is.na(x) & x > 0
  det_frac <- rowMeans(detected)
  
  det_ok <- det_frac >= MIN_DETECT_FRAC
  x   <- x[det_ok, , drop = FALSE]
  f_m <- f_m[det_ok, ]
  n2 <- nrow(x)
  cat("  after detected in >=", MIN_DETECT_FRAC, "of samples:", n2, "\n")
  
  for (thr in DETECT_SENSITIVITY)
    cat("    [>=", thr, "would keep", sum(det_frac >= thr), "]\n")
  
  # ---- 2c. imputation ----------------------------------
  # Missingness is LEFT-CENSORED: a value is absent because it
  # fell below the detection limit, not at random. Mean or
  # median imputation would bias those values upward.
  
  x[x == 0] <- NA
  n_missing <- sum(is.na(x))
  cat("  cells to impute:", n_missing,
      sprintf("(%.1f%%)", 100 * n_missing / length(x)), "\n")
  
  if (IMPUTE_METHOD == "half_min") {
    for (i in seq_len(nrow(x))) {
      v <- x[i, ]
      if (any(is.na(v))) v[is.na(v)] <- min(v, na.rm = TRUE) / 2
      x[i, ] <- v
    }
  } else if (IMPUTE_METHOD == "qrilc") {
    if (!requireNamespace("imputeLCMD", quietly = TRUE))
      stop("install imputeLCMD, or set IMPUTE_METHOD <- 'half_min'")
    lg <- log2(t(x) + PSEUDOCOUNT)
    lg <- imputeLCMD::impute.QRILC(lg)[[1]]
    x  <- t(2^lg - PSEUDOCOUNT)
  } else stop("unknown IMPUTE_METHOD")
  
  stopifnot(!any(is.na(x)))
  
  # ---- 2d. log2 transform ------------------------------
  # Total median scaling was applied upstream - do NOT repeat.
  x <- log2(x + PSEUDOCOUNT)
  cat("  value range after log2:", round(range(x), 2), "\n")
  
  processed[[m]] <- list(x = x, info = f_m)
  
  funnel <- rbind(
    funnel,
    data.frame(method = m, step = "0. input",     n = n0),
    data.frame(method = m, step = "1. QC CV",     n = n1),
    data.frame(method = m, step = "2. detection", n = n2))
}


# ============================================================
# 3. MERGE THE FOUR METHODS
# ============================================================

X  <- do.call(rbind, lapply(processed, function(p) p$x))
FI <- do.call(rbind, lapply(processed, function(p) p$info))
stopifnot(identical(rownames(X), FI$feature_id))

cat("\n\n==================================================\n")
cat("MERGED:", nrow(X), "features x", ncol(X), "samples\n")


# ---- 3a. what counts as a real HMDB identifier ------------
# The source column contains blanks, the literal string "NA",
# and "redundant ion" notes. Only accept things that look like
# accessions: HMDB + digits, optional trailing *.

real_hmdb <- !is.na(FI$hmdb) &
  !(trimws(FI$hmdb) %in% BAD_ID) &
  grepl("^HMDB[0-9]+\\*?$", trimws(FI$hmdb))

cat("\n=== HMDB ID CHECK ===\n")
cat("non-blank hmdb values:", sum(!is.na(FI$hmdb) & FI$hmdb != ""), "\n")
cat("that look like real accessions:", sum(real_hmdb), "\n")

rejected <- unique(FI$hmdb[!is.na(FI$hmdb) & FI$hmdb != "" & !real_hmdb])
if (length(rejected) > 0) {
  cat("rejected values:\n")
  print(head(rejected, 20))
}

FI$hmdb_clean <- ifelse(real_hmdb, sub("\\*$", "", trimws(FI$hmdb)), NA)
FI$has_hmdb   <- real_hmdb


# ---- 3b. duplicate compounds across methods ---------------
# Rule, applied mechanically: keep the instance with the LOWER
# pooled-QC CV. No case-by-case judgement.

ids  <- FI$hmdb_clean
dups <- unique(ids[!is.na(ids)][duplicated(ids[!is.na(ids)])])

cat("\n=== DEDUP ===\n")
cat("compounds detected on more than one platform:", length(dups), "\n")

drop <- character(0)
for (h in dups) {
  rows <- which(!is.na(ids) & ids == h)
  best <- rows[which.min(FI$qc_cv[rows])]
  drop <- c(drop, FI$feature_id[setdiff(rows, best)])
  cat("  ", h, ": kept", FI$feature_id[best],
      "(CV", round(FI$qc_cv[best], 3), ")\n")
}

if (length(drop) > 0) {
  k  <- !(FI$feature_id %in% drop)
  X  <- X[k, , drop = FALSE]
  FI <- FI[k, ]
}

cat("features dropped as duplicates:", length(drop), "\n")
cat("after dedup:", nrow(X), "features\n")


# ============================================================
# 4. KEGG  -  one bulk download, then join on HMDB
# ============================================================
#
# Two traps handled here:
#
# (a) refmet/<name>/all has NO kegg_id field - it returns name,
#     regno, pubchem_cid, inchi_key, exactmass, formula,
#     sys_name, main_class, sub_class, synonyms. The endpoint
#     that carries external IDs is refmet/all_ids. It returns a
#     named JSON list whose records have DIFFERING fields, so it
#     must be parsed field-by-field, not rbind-ed.
#
# (b) Our HMDB accessions use the legacy 9-character format
#     (HMDB00535); RefMet uses the current 11-character one
#     (HMDB0000535). Join on the numeric portion.

kegg_file <- "data/raw/refmet_all_ids.txt"

FI$kegg <- NA_character_

norm_hmdb <- function(x) {
  n <- sub("^HMDB0*", "", toupper(trimws(x)))
  n[n == "" | is.na(x)] <- NA
  n
}

if (DO_KEGG) {
  
  if (!file.exists(kegg_file)) {
    cat("\ndownloading RefMet external ID table...\n")
    options(timeout = 600)
    download.file("https://www.metabolomicsworkbench.org/rest/refmet/all_ids",
                  kegg_file, mode = "wb")
  }
  cat("\nRefMet ID table:", round(file.size(kegg_file) / 1e6, 2), "MB\n")
  
  if (!requireNamespace("jsonlite", quietly = TRUE))
    install.packages("jsonlite")
  
  j <- jsonlite::fromJSON(kegg_file, simplifyVector = FALSE)
  cat("RefMet records:", length(j), "\n")
  
  grab <- function(x, f) {
    v <- x[[f]]
    if (is.null(v) || length(v) == 0) NA_character_ else as.character(v)[1]
  }
  
  ref <- data.frame(
    hmdb = vapply(j, grab, character(1), "hmdb_id"),
    kegg = vapply(j, grab, character(1), "kegg_id"),
    stringsAsFactors = FALSE
  )
  ref$hmdb[ref$hmdb %in% BAD_ID] <- NA
  ref$kegg[ref$kegg %in% BAD_ID] <- NA
  
  cat("  records with an HMDB ID:", sum(!is.na(ref$hmdb)), "\n")
  cat("  records with a KEGG ID: ", sum(!is.na(ref$kegg)), "\n")
  
  FI$kegg <- ref$kegg[match(norm_hmdb(FI$hmdb_clean), norm_hmdb(ref$hmdb))]
  
  cat("\nspot-check - verify these by eye:\n")
  print(head(FI[FI$has_hmdb & !is.na(FI$kegg),
                c("feature_id", "metabolite", "hmdb_clean", "kegg")], 8))
}


# ============================================================
# 5. IDENTIFIER SUMMARY
# ============================================================

n_named <- sum(!is.na(FI$metabolite) & FI$metabolite != "")

cat("\n=== IDENTIFIERS ===\n")
cat("features:                   ", nrow(FI), "\n")
cat("with a metabolite name:     ", n_named, "\n")
cat("with an HMDB ID:            ", sum(FI$has_hmdb), "\n")
cat("named but no HMDB ID:       ",
    sum(!is.na(FI$metabolite) & FI$metabolite != "" & !FI$has_hmdb), "\n")
cat("with a KEGG ID:             ", sum(!is.na(FI$kegg)), "\n")
cat("HMDB present, KEGG missing: ",
    sum(FI$has_hmdb & is.na(FI$kegg)), "of", sum(FI$has_hmdb), "\n")


# ============================================================
# 6. EMIT THE TWO MATRICES
# ============================================================

to_df <- function(m) cbind(feature_id = rownames(m),
                           as.data.frame(m), stringsAsFactors = FALSE)

write_parquet(to_df(X), "data/derived/features_full_v1.parquet")
write_parquet(FI,       "data/derived/features_full_v1_info.parquet")

ann <- FI$has_hmdb
write_parquet(to_df(X[ann, , drop = FALSE]),
              "data/derived/features_annotated_v1.parquet")
write_parquet(FI[ann, ],
              "data/derived/features_annotated_v1_info.parquet")

saveRDS(FI, "data/derived/feature_info_final.rds")

cat("\n=== WRITTEN ===\n")
cat("features_full_v1.parquet      :", nrow(X), "x", ncol(X), "\n")
cat("features_annotated_v1.parquet :", sum(ann), "x", ncol(X), "\n")


# ============================================================
# 7. FUNNEL
# ============================================================

overall <- data.frame(
  step = c("0. input",
           paste0("1. QC CV <= ", CV_THRESHOLD),
           paste0("2. detected >= ", MIN_DETECT_FRAC),
           "3. dedup across methods",
           "4. HMDB annotated",
           "5. KEGG mapped"),
  n = c(nrow(mat),
        sum(funnel$n[funnel$step == "1. QC CV"]),
        sum(funnel$n[funnel$step == "2. detection"]),
        nrow(X),
        sum(FI$has_hmdb),
        sum(!is.na(FI$kegg))),
  stringsAsFactors = FALSE
)

cat("\n=== FUNNEL ===\n")
print(overall)

write.csv(overall, "output/funnel.csv", row.names = FALSE)
write.csv(funnel,  "output/funnel_by_method.csv", row.names = FALSE)

png("output/funnel.png", width = 900, height = 600)
par(mar = c(5, 16, 4, 6))
bp <- barplot(rev(overall$n), horiz = TRUE, log = "x",
              names.arg = rev(overall$step), las = 1,
              xlab = "features (log scale)",
              main = "Feature filtering funnel", col = "salmon")
text(rev(overall$n), bp, labels = rev(overall$n), pos = 4, xpd = TRUE)
dev.off()

cat("\nsaved output/funnel.png\n")
cat("\n=== DONE ===\n")

