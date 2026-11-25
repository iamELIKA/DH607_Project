# ============================================================
# 07_sample_table.R   -   Person B
# Sample table, timeline, episode alignment, splits
# ============================================================
#
# RESTART R BEFORE RUNNING.
#
# Needs: data/raw/hmp2_metadata.csv
#        data/derived/mbx_matrix.rds
#        data/derived/dysbiosis.rds      (script 03)
#        data/derived/episodes.rds       (script 04)
#
# Produces: data/derived/samples.rds
#           output/samples.csv           <- the project backbone
#           output/split.csv
#           output/data_dictionary.csv
#           output/cohort_flow.png, cohort_flow.csv
# ============================================================

setwd("/home/ritwikan/metabolomics")

SEED <- 607          # fixed, never changed
TEST_FRACTION <- 0.30

dir.create("output", showWarnings = FALSE)


# ============================================================
# 0. LOAD
# ============================================================

meta <- read.csv("data/raw/hmp2_metadata.csv", check.names = FALSE,
                 stringsAsFactors = FALSE)
mat  <- readRDS("data/derived/mbx_matrix.rds")
dys  <- readRDS("data/derived/dysbiosis.rds")
eps  <- readRDS("data/derived/episodes.rds")

mbx_ids <- colnames(mat)
cat("metabolomic samples:", length(mbx_ids), "\n")


# helper: pull a column if it exists, else NA
col <- function(d, name) {
  if (name %in% colnames(d)) d[[name]] else rep(NA, nrow(d))
}


# ============================================================
# 1. ONE ROW PER METABOLOMIC SAMPLE
# ============================================================

mm <- meta[meta$data_type == "metabolomics", ]
i  <- match(mbx_ids, mm[["External ID"]])

stopifnot(!any(is.na(i)))

s <- data.frame(
  sample_id   = mbx_ids,
  subject_id  = mm[["Participant ID"]][i],
  collection  = col(mm, "site_sub_coll")[i],
  week        = col(mm, "week_num")[i],
  visit_num   = col(mm, "visit_num")[i],
  diagnosis   = col(mm, "diagnosis")[i],
  site        = col(mm, "site_name")[i],
  sex         = col(mm, "sex")[i],
  age         = col(mm, "consent_age")[i],
  race        = col(mm, "race")[i],
  stringsAsFactors = FALSE
)

# ---- 1a. timeline in DAYS, not visit index ----------------
s$day <- s$week * 7

cat("missing week_num:", sum(is.na(s$day)), "\n")

s <- s[order(s$subject_id, s$day), ]

# position within the subject's own series
s$visit_index <- ave(s$day, s$subject_id, FUN = seq_along)
s$days_since_first <- ave(s$day, s$subject_id, FUN = function(d) d - min(d))
s$gap_from_prev <- ave(s$day, s$subject_id,
                       FUN = function(d) c(NA, diff(d)))

cat("subjects:", length(unique(s$subject_id)), "\n")


# ============================================================
# 2. DISEASE ACTIVITY MEASURES
# ============================================================
#
# Coverage is uneven - that is reported, not hidden.

s$fecalcal <- col(mm, "fecalcal")[i]
s$hbi      <- col(mm, "hbi")[i]
s$sccai    <- col(mm, "sccai")[i]

cat("\n=== ACTIVITY MEASURE COVERAGE ===\n")
cat("faecal calprotectin:", sum(!is.na(s$fecalcal)), "of", nrow(s), "\n")
cat("HBI (CD only):      ", sum(!is.na(s$hbi)), "\n")
cat("SCCAI (UC only):    ", sum(!is.na(s$sccai)), "\n")

cat("\nby diagnosis:\n")
print(data.frame(
  diagnosis = names(table(s$diagnosis)),
  n         = as.integer(table(s$diagnosis)),
  fecalcal  = as.integer(tapply(!is.na(s$fecalcal), s$diagnosis, sum)),
  hbi       = as.integer(tapply(!is.na(s$hbi),      s$diagnosis, sum)),
  sccai     = as.integer(tapply(!is.na(s$sccai),    s$diagnosis, sum))
))


# ============================================================
# 3. MEDICATION COVARIATES
# ============================================================
#
# These come from the two-week recall questionnaire paired with
# each stool sample, so they are sample-level and complete.
#
# NOTE: brand-level medication columns (Asacol, Pentasa, Lialda,
# Apriso, Prednisone, Entocort and the rest) are recorded only at
# the five quarterly CLINIC visits. Metabolomic samples are mailed
# stool kits, so those columns are empty for all 546 of them -
# verified, not assumed. They are deliberately not used.

s$antibiotics       <- col(mm, "Antibiotics")[i]
s$immunosuppressant <-
  col(mm, "Immunosuppressants (e.g. oral corticosteroids)")[i]
s$chemotherapy      <- col(mm, "Chemotherapy")[i]

cat("\n=== MEDICATION COVERAGE ===\n")
for (v in c("antibiotics", "immunosuppressant", "chemotherapy")) {
  cat(sprintf("  %-18s non-missing: %4d\n", v, sum(!is.na(s[[v]]))))
}

cat("\nantibiotics:\n");       print(table(s$antibiotics, useNA = "ifany"))
cat("\nimmunosuppressant:\n"); print(table(s$immunosuppressant, useNA = "ifany"))


# ============================================================
# 4. DYSBIOSIS, ALIGNED BY COLLECTION EVENT
# ============================================================
#
# IMPORTANT: External ID is unique PER ASSAY. A stool sample has
# one ID for metabolomics and a different one for metagenomics.
# They are linked by site_sub_coll, the collection event.

mg <- meta[meta$data_type == "metagenomics", ]
dys$collection <- mg$site_sub_coll[match(dys$sample_id, mg[["External ID"]])]

j <- match(s$collection, dys$collection)

s$dysbiosis_score <- dys$dysbiosis_score[j]
s$dysbiotic       <- dys$dysbiotic[j]

cat("\n=== DYSBIOSIS LINKAGE ===\n")
cat("metabolomic samples with a same-collection dysbiosis score:",
    sum(!is.na(s$dysbiosis_score)), "of", nrow(s), "\n")
cat("of those, dysbiotic:", sum(s$dysbiotic, na.rm = TRUE), "\n")


# ============================================================
# 5. EPISODE ALIGNMENT
# ============================================================
#
# For each metabolomic sample, where does it sit relative to its
# subject's dysbiotic episodes? These columns are what Q1, Q3
# and Q4 are actually built on.

s$in_episode            <- FALSE
s$episode_id            <- NA_integer_
s$days_to_next_onset    <- NA_real_
s$days_since_resolution <- NA_real_
s$episode_phase         <- NA_character_   # pre / during / post

eps$episode_id <- seq_len(nrow(eps))

for (k in seq_len(nrow(s))) {
  
  e <- eps[eps$subject_id == s$subject_id[k], ]
  if (nrow(e) == 0) next
  d <- s$day[k]
  
  # inside an episode?
  hit <- which(d >= e$onset_day & d <= e$last_dys_day)
  if (length(hit) > 0) {
    s$in_episode[k]    <- TRUE
    s$episode_id[k]    <- e$episode_id[hit[1]]
    s$episode_phase[k] <- "during"
  }
  
  # time to the next onset (for Q1 / Q3)
  future <- e$onset_day[e$onset_day > d]
  if (length(future) > 0) {
    s$days_to_next_onset[k] <- min(future) - d
    if (is.na(s$episode_phase[k])) s$episode_phase[k] <- "pre"
  }
  
  # time since the most recent resolution (for Q4)
  past <- e$resolution_day[!is.na(e$resolution_day) & e$resolution_day <= d]
  if (length(past) > 0) {
    s$days_since_resolution[k] <- d - max(past)
    if (is.na(s$episode_phase[k])) s$episode_phase[k] <- "post"
  }
}

cat("\n=== EPISODE ALIGNMENT ===\n")
cat("samples inside an episode:  ", sum(s$in_episode), "\n")
cat("samples with a future onset:", sum(!is.na(s$days_to_next_onset)), "\n")
cat("samples after a resolution: ", sum(!is.na(s$days_since_resolution)), "\n")
cat("\nphase:\n")
print(table(s$episode_phase, useNA = "ifany"))
cat("\n(NA = subject has no episodes at all)\n")


# ============================================================
# 6. SUBJECT-LEVEL TRAIN/TEST SPLIT
# ============================================================
#
# SUBJECT level, never sample level. Sample-level splitting on
# this cohort is a documented source of inflated AUCs.
# Stratified by diagnosis so all three groups appear in both.

set.seed(SEED)

subj <- unique(s[, c("subject_id", "diagnosis")])
subj <- subj[order(subj$subject_id), ]

subj$split <- NA_character_
for (dx in unique(subj$diagnosis)) {
  idx <- which(subj$diagnosis == dx)
  n_test <- max(1, round(TEST_FRACTION * length(idx)))
  test_i <- sample(idx, n_test)
  subj$split[idx]    <- "train"
  subj$split[test_i] <- "test"
}

s$split <- subj$split[match(s$subject_id, subj$subject_id)]

cat("\n=== SPLIT (seed", SEED, ") ===\n")
cat("subjects:\n");  print(table(subj$diagnosis, subj$split))
cat("\nsamples:\n"); print(table(s$diagnosis, s$split))

# no subject may appear in both
stopifnot(length(intersect(
  s$subject_id[s$split == "train"],
  s$subject_id[s$split == "test"])) == 0)
cat("\nleakage check passed: no subject in both sets\n")


# ============================================================
# 7. SAVE
# ============================================================

saveRDS(s, "data/derived/samples.rds")
write.csv(s, "output/samples.csv", row.names = FALSE)
write.csv(subj, "output/split.csv", row.names = FALSE)

cat("\nwrote output/samples.csv :", nrow(s), "rows x", ncol(s), "cols\n")


# ---- data dictionary --------------------------------------

dd <- data.frame(
  column = c(
    "sample_id","subject_id","collection","week","day","visit_num",
    "visit_index","days_since_first","gap_from_prev",
    "diagnosis","site","sex","age","race",
    "fecalcal","hbi","sccai",
    "antibiotics","immunosuppressant","chemotherapy",
    "dysbiosis_score","dysbiotic",
    "in_episode","episode_id","episode_phase",
    "days_to_next_onset","days_since_resolution","split"),
  description = c(
    "Metabolomics External ID; column name in the feature matrix",
    "Participant ID",
    "site_sub_coll: the collection event, links assays to each other",
    "week_num from enrolment",
    "week * 7; the time axis for all longitudinal analysis",
    "Clinic visit number where applicable",
    "Order of this sample within the subject's own series",
    "Days from the subject's first metabolomic sample",
    "Days since the subject's previous metabolomic sample",
    "CD / UC / nonIBD",
    "Recruitment site",
    "Sex",
    "Age at consent",
    "Race",
    "Faecal calprotectin",
    "Harvey-Bradshaw Index (CD)",
    "Simple Clinical Colitis Activity Index (UC)",
    "Antibiotics in the past two weeks (Yes/No)",
    "Immunosuppressants in the past two weeks (Yes/No)",
    "Chemotherapy in the past two weeks (Yes/No)",
    "Median Bray-Curtis to non-IBD, from the same collection event",
    "dysbiosis_score > 0.8427",
    "Sample falls inside a dysbiotic episode",
    "Which episode, if inside one",
    "pre / during / post relative to episodes; NA if subject has none",
    "Days until this subject's next episode onset (Q1, Q3)",
    "Days since this subject's most recent resolution (Q4)",
    "train / test, assigned at SUBJECT level, seed 607"),
  stringsAsFactors = FALSE
)

stopifnot(all(dd$column %in% colnames(s)))
stopifnot(all(colnames(s) %in% dd$column))

write.csv(dd, "output/data_dictionary.csv", row.names = FALSE)
cat("wrote output/data_dictionary.csv (", nrow(dd), "columns documented )\n")


# ============================================================
# 8. COHORT FLOW FIGURE
# ============================================================

flow <- data.frame(
  stage = c("subjects enrolled",
            "with metabolomics",
            "with >=3 timepoints",
            "with >=1 dysbiotic episode"),
  n = c(length(unique(meta[["Participant ID"]])),
        length(unique(s$subject_id)),
        sum(table(s$subject_id) >= 3),
        length(unique(eps$subject_id[eps$subject_id %in% s$subject_id]))),
  stringsAsFactors = FALSE
)

cat("\n=== COHORT FLOW ===\n")
print(flow)

cat("\nsample-level: ", sum(!is.na(s$dysbiosis_score)), "of", nrow(s),
    "metabolomic samples have a same-collection dysbiosis score\n")

write.csv(flow, "output/cohort_flow.csv", row.names = FALSE)

png("output/cohort_flow.png", width = 900, height = 550)
par(mar = c(5, 16, 4, 6))
bp <- barplot(rev(flow$n), horiz = TRUE, names.arg = rev(flow$stage),
              las = 1, xlab = "subjects", col = "lightgreen",
              main = "Cohort flow", xlim = c(0, max(flow$n) * 1.15))
text(rev(flow$n), bp, labels = rev(flow$n), pos = 4, xpd = TRUE)
dev.off()

cat("\nsaved output/cohort_flow.png\n")
cat("\n=== DONE ===\n")
cat("samples.csv is the backbone - every downstream script joins to it.\n")


table(readRDS("data/derived/samples.rds")$chemotherapy)
