# IBD Microbial Transcriptional Activity & Multi-Omics Candidate Prioritisation

**Project Title:** Prioritisation of Microbial Protein Candidates Associated with Host Responses in Inflammatory Bowel Disease Using Metatranscriptomic and Multi-omics Data  
**Dataset:** IBDMDB / HMP2 Longitudinal Cohort ([ibdmdb.org](https://ibdmdb.org))  
**Reference:** Lloyd-Price et al. (2019), *Nature* 569, 655–662 (HMP2 Multi-omics Study)

---

## 1. Research Objective

This repository implements a reproducible R analysis pipeline investigating microbial gene-family and enzyme-level transcription in **Crohn's Disease (CD)** and **Ulcerative Colitis (UC)** compared to **non-IBD controls**, while distinguishing **active transcriptional regulation** (RNA shifts per gene copy) from **passive abundance shifts** (DNA copy-number changes) using paired metatranscriptomics (MTX) and metagenomics (MGX), alongside host mucosal transcriptomics (HTX) and upcoming metabolomic integration.

---

## 2. Repository Structure

```text
IBD_MTX_project/
├── README.md                          # Project overview, pipeline guide & results index
├── .gitignore                         # Excludes large raw/processed binaries & OS artifacts
│
├── scripts/                           # Numbered, sequential R analysis scripts
│   ├── 01_download_and_audit.R        # Download HMP2 MTX/metadata files & audit structure/units
│   ├── 02_metadata_and_cohort_QC.R    # Cohort QC, confounder audit & MTX–MGX visit pairing
│   ├── 03_MTX_differential_analysis.R # Unstratified UniRef90 differential transcription + sensitivity
│   ├── 04_MGX_adjustment_and_transcriptional_activity.R
│   │                                  # Paired RNA+DNA (EC-level) active vs passive modeling
│   └── h1_download_and_visualise_HTX.R# Host biopsy transcriptomics QC, PCA & top variable genes
│
├── metabolomics/                      # [Reserved] Empty folder for co-researcher's metabolomic analysis
│
├── figures/                           # All QC, differential, QQ, DNA-adjusted & HTX plots
│   ├── qc_*.pdf                       # Script 02 cohort & multi-omics overlap plots
│   ├── volcano_CD_vs_nonIBD.png       # Script 03 gene-family volcano (CD vs nonIBD)
│   ├── volcano_UC_vs_nonIBD.png       # Script 03 gene-family volcano (UC vs nonIBD)
│   ├── qq_CD_vs_nonIBD.png            # Script 03 QQ diagnostic plot (CD, λ = 8.04)
│   ├── qq_UC_vs_nonIBD.png            # Script 03 QQ diagnostic plot (UC, λ = 6.50)
│   ├── rna_vs_dna_scatter_CD.png      # Script 04 unadjusted vs DNA-adjusted effect sizes (CD)
│   ├── volcano_MGX_adjusted_CD.png    # Script 04 DNA-adjusted EC volcano (CD vs nonIBD)
│   ├── volcano_MGX_adjusted_UC.png    # Script 04 DNA-adjusted EC volcano (UC vs nonIBD)
│   └── htx_*.png                      # Script h1 host mucosal biopsy QC, PCA, heatmap & gene plots
│
├── results/
│   ├── qc/                            # Cohort summary CSV & metadata/HTX QC text reports
│   └── differential/                  # Script 04 MGX-adjusted EC tables & sensitivity summary
│
└── logs/                              # Reproducibility audit logs, run logs & R sessionInfo()
```
*(Note: `data/raw/`, `data/processed/`, and `data/metadata/` are auto-created locally when running Script 01 and excluded via `.gitignore`.)*

---

## 3. Pipeline Summary & Key Findings (Scripts 01–04 + HTX)

| Step | Script | Scope & Cohort | Key Outputs & Findings |
|---|---|---|---|
| **01** | `scripts/01_download_and_audit.R` | Raw `MTX_genefamilies.tsv.gz` (2,164,739 rows × 736 cols) & metadata (5,533 rows × 490 cols) | Verified file checksums; identified that MTX values are **relative abundance (0–1)** rather than RPK despite `_Abundance-RPKs` headers; separated 917,164 unstratified vs 1,247,575 taxon-stratified features. |
| **02** | `scripts/02_metadata_and_cohort_QC.R` | **835 MTX samples** across **109 participants** (CD=394 / 52, UC=238 / 30, nonIBD=203 / 27; median 7 visits/person) | Characterized confounders (age, site, antibiotics: CD 15.0%, UC 4.6%, nonIBD 3.0%); built same-visit MTX–MGX pairing table (`MTX_MGX_same_visit_pairs.csv`) and 4 cohort QC plots. |
| **03** | `scripts/03_MTX_differential_analysis.R` | **79,563** filtered unstratified UniRef90 gene families (`≥10%` prevalence, mean `≥1e-7`) across 706 samples | Adjusted for age, sex, site & antibiotics (`CD`: 31,885 FDR < 0.05; `UC`: 28,979 FDR < 0.05; 98.6% overlap at 20% prevalence sensitivity). **Diagnostic finding:** 67.7% sparsity/singular fits + `lm()` fallback yielded inflated QQ $\lambda$ (`8.04` CD, `6.50` UC), motivating DNA-adjusted modeling in Script 04. |
| **04** | `scripts/04_MGX_adjustment_and_transcriptional_activity.R` | **1,107** shared Enzyme Commission (EC) features across **735 paired MTX+MGX samples** (CD=337, UC=211, nonIBD=187) | Modeled `log2(RNA) ~ log2(DNA) + diagnosis + age + sex + site + antibiotics` to isolate per-copy regulation:<br>• **CD vs nonIBD:** **79** Active (7.1%), **82** Passive (7.4%), 946 NS (85.5%)<br>• **UC vs nonIBD:** **181** Active (16.4%), **126** Passive (11.4%), 800 NS (72.3%) |
| **HTX** | `scripts/h1_download_and_visualise_HTX.R` | **249** host mucosal biopsies, **21,768** expressed genes (median depth 39.3M reads; overlaps 88/109 MTX participants) | PCA shows **PC1 (35.7%) = tissue location** (Ileum vs Colon/Rectum) and **PC2 (23.1%) = inflammation**; identifies host lipid/bile-acid transport (`APOA1/4`, `APOB`, `FABP6`) & antimicrobial (`DEFA5/6`, `REG3A`) anchors. |

---

## 4. Plots Index (`figures/`)

| Plot File | Generated By | Description |
|---|---|---|
| `qc_diagnosis_distribution.pdf` | Script 02 | MTX sample counts across CD, UC, and nonIBD |
| `qc_samples_per_participant.pdf` | Script 02 | Longitudinal visit counts per participant |
| `qc_multiomics_overlap.pdf` | Script 02 | Participant overlap across MTX, MGX, MPX, MBX, and HTX |
| `qc_CD_HBI_over_time.pdf` | Script 02 | Longitudinal Harvey–Bradshaw Index (HBI) trajectories in CD |
| `volcano_CD_vs_nonIBD.png` / `volcano_UC_vs_nonIBD.png` | Script 03 | Unadjusted-for-DNA gene-family differential transcription |
| `qq_CD_vs_nonIBD.png` / `qq_UC_vs_nonIBD.png` | Script 03 | QQ diagnostics showing $\lambda$ inflation from repeated-measures `lm()` fallback |
| `rna_vs_dna_scatter_CD.png` | Script 04 | Unadjusted vs DNA-adjusted EC effect sizes (Active vs Passive) |
| `volcano_MGX_adjusted_CD.png` / `volcano_MGX_adjusted_UC.png` | Script 04 | DNA-adjusted EC volcano plots separating Active (red) from Passive (blue) |
| `htx_qc_library_size.png` | Script h1 | Host biopsy sequencing depth across diagnosis & tissue site |
| `htx_pca_diagnosis_tissue.png` | Script h1 | Host transcriptome PCA (tissue site & inflammation axes) |
| `htx_heatmap_top50_genes.png` | Script h1 | Top 50 most variable host mucosal genes |
| `htx_barplot_top_variable_genes.png` | Script h1 | Mean `log2-CPM` (± SE) for top variable host epithelial/immune genes |

---

## 5. Co-Researcher Collaboration (`metabolomics/`)

The **`metabolomics/`** folder is intentionally initialized as an empty directory for integrating the **metabolomic (MBX) analysis** performed by a co-researcher (e.g., longitudinal deviation/pre-flare dynamics and metabolite signatures).

---

## 6. Reproducibility & Execution

* **R Version:** ≥ 4.2 (server runs logged in `logs/` used R 4.6.2; `set.seed(42)` in all scripts).
* **Key R Packages:** `data.table`, `dplyr`, `tidyr`, `ggplot2`, `scales`, `ggrepel`, `parallel`, `lme4`, `lmerTest`, `pheatmap`.
* **Large Data Policy:** Raw `.tsv.gz` matrices, processed `.rds` objects, and full 79k-row gene-family CSV tables are excluded via `.gitignore` to keep the repository lightweight. Run `scripts/01_download_and_audit.R` (and download `MTX_ecs.tsv.gz` / `MGX_ecs.tsv.gz` from IBDMDB into `data/raw/`) to regenerate local inputs, then execute Scripts `01` → `02` → `03` → `04` in order.
