# DH607 Multi-Omics IBD Project: Transcriptomics & Metabolomics

**Dataset:** IBDMDB / HMP2 Longitudinal Cohort ([ibdmdb.org](https://ibdmdb.org))  
**Reference:** Lloyd-Price et al. (2019), *Nature* 569, 655–662 (HMP2 Multi-omics Study)

---

## 1. Project Overview

This repository investigates Inflammatory Bowel Disease (**Crohn's Disease [CD]** and **Ulcerative Colitis [UC]** vs. **non-IBD controls**) across two multi-omics branches:

1. **`transcriptomics/`** — Split into **`microbial/`** (stool metatranscriptomics [MTX] + paired metagenomic [MGX] DNA-adjustment) and **`htx/`** (host mucosal biopsy transcriptomics [HTX]).
2. **`metabolomics/`** — Longitudinal stool metabolomics (MBX) examining pre-dysbiosis volatility, personal-baseline deviation vs. absolute abundance, lead-time prediction, and post-flare recovery.

---

## 2. Repository Structure

```text
DH607_Project/
├── README.md                  # Unified project overview & index
├── .gitignore                 # Excludes raw/processed binaries, archives, logs & local docs
│
├── transcriptomics/
│   ├── microbial/             # Microbial MTX & Paired MGX Analysis
│   │   ├── scripts/           # R pipeline scripts (01–04)
│   │   ├── figures/           # Cohort QC, volcano, QQ & DNA-adjusted plots
│   │   └── results/           # QC summaries & MGX-adjusted EC tables
│   │
│   └── htx/                   # Host Mucosal Biopsy Transcriptomics (HTX)
│       ├── scripts/           # Host biopsy QC & exploration script (h1)
│       ├── figures/           # Library depth QC, PCA, heatmap & variable-gene plots
│       └── results/           # Host transcriptomics QC report
│
└── metabolomics/              # Longitudinal Stool Metabolomics (MBX)
    ├── scripts/               # R analysis scripts (step_0, personA/B/C, Q1–Q4)
    ├── csv files/             # Summary tables, test statistics & trajectories
    └── figures/               # Volatility, deviation, lead-time & recovery plots
```

*(Note: Local `data/raw/`, `data/processed/`, and `data/metadata/` folders are auto-created when running the download scripts and excluded via `.gitignore`.)*

---

## 3. Transcriptomics Summary (`transcriptomics/`)

### Microbial (`transcriptomics/microbial/`)
| Step | Script | Scope & Cohort | Key Outputs & Findings |
|---|---|---|---|
| **01** | `01_download_and_audit.R` | `MTX_genefamilies.tsv.gz` (2,164,739 rows × 736 cols) & metadata (5,533 × 490) | Verified checksums; discovered MTX values are **relative abundance (0–1)** rather than RPK; split 917,164 unstratified vs 1,247,575 taxon-stratified features. |
| **02** | `02_metadata_and_cohort_QC.R` | **835 MTX samples** across **109 participants** (CD=394/52, UC=238/30, nonIBD=203/27; median 7 visits) | Characterized confounders (age, site, antibiotics: CD 15.0%, UC 4.6%, nonIBD 3.0%); built same-visit MTX–MGX pairing table & 4 QC plots. |
| **03** | `03_MTX_differential_analysis.R` | **79,563** filtered unstratified UniRef90 gene families (`≥10%` prevalence, mean `≥1e-7`) across 706 samples | CD: 31,885 FDR < 0.05; UC: 28,979 FDR < 0.05 (98.6% overlap at 20% prevalence sensitivity). Sparse singular LMM fits / `lm()` fallback inflated QQ $\lambda$ (`8.04` CD, `6.50` UC), motivating Script 04. |
| **04** | `04_MGX_adjustment_and_transcriptional_activity.R` | **1,107** shared ECs across **735 paired MTX+MGX samples** (CD=337, UC=211, nonIBD=187) | `log2(RNA) ~ log2(DNA) + diagnosis + covariates`:<br>• **CD vs nonIBD:** **79** Active (7.1%), **82** Passive (7.4%), 946 NS (85.5%)<br>• **UC vs nonIBD:** **181** Active (16.4%), **126** Passive (11.4%), 800 NS (72.3%) |

### Host (`transcriptomics/htx/`)
| Step | Script | Scope & Cohort | Key Outputs & Findings |
|---|---|---|---|
| **h1** | `h1_download_and_visualise_HTX.R` | **249** host biopsies, **21,768** expressed genes (median 39.3M reads; 88/109 MTX overlap) | PCA: **PC1 (35.7%) = tissue site** (Ileum vs Colon/Rectum), **PC2 (23.1%) = inflammation**; highlights lipid/bile-acid (`APOA1/4`, `APOB`, `FABP6`) & antimicrobial (`DEFA5/6`, `REG3A`) genes. |

---

## 4. Metabolomics Summary (`metabolomics/`)

| Component | Contents | Focus |
|---|---|---|
| `metabolomics/scripts/` | `step_0.R`, `personA.R`, `PersonB.R`, `personC.R`, `Q1_volatility.R`, `Q2_NEW.R`, `Q3.R`, `q4.R` | Cohort preprocessing and longitudinal Q1–Q4 modeling |
| `metabolomics/csv files/` | `Q1_*`–`Q4_*` result CSVs & method summaries | Pre-onset volatility (45–75d window), ABS vs DEV CV, lead-time bins, and post-flare distance/trajectories |
| `metabolomics/figures/` | `Q1_*`–`Q4_*`, `Figure4_storyline.png`, `deviation_checks.png`, `volatility_checks.png` | Publication/presentation figures for Q1–Q4 metabolomics findings |

---

## 5. Reproducibility

* **R Version:** ≥ 4.2 (`set.seed(42)` in transcriptomics scripts).
* **Large Data Policy:** Raw `.tsv.gz` matrices, processed `.rds` objects, `.zip` archives, and 15 MB gene-family CSV dumps are excluded via `.gitignore`.
