# Technical Documentation: Multi-Omics Integrative Analysis Pipeline (`COL3A1` 2D vs 3D)

Comprehensive technical reference for the R bioinformatics pipeline designed to analyze, integrate, and biologically interpret **RNA-seq transcriptomics** and **mass spectrometry proteomics** datasets from **2D** and **3D** cellular models associated with the **`COL3A1`** gene.

---

## Table of Contents
- [Pipeline Architecture](#pipeline-architecture)
- [Prerequisites and Installation](#prerequisites-and-installation)
- [Input Data Structure](#input-data-structure)
- [Analytical Workflow](#analytical-workflow)
  - [1. Pre-processing and Normalization](#1-pre-processing-and-normalization)
  - [2. Unsupervised Multi-Omics Integration (MOFA2)](#2-unsupervised-multi-omics-integration-mofa2)
  - [3. Supervised Classification & Integration (sPLS-DA & DIABLO)](#3-supervised-classification--integration-spls-da--diablo)
  - [4. Functional Enrichment Analysis (GOEA & GSEA)](#4-functional-enrichment-analysis-goea--gsea)
  - [5. Automated Reporting and LaTeX Export](#5-automated-reporting-and-latex-export)
- [Core Custom Functions](#core-custom-functions)
- [Generated Outputs](#generated-outputs)
- [Reproducibility Notes](#reproducibility-notes)

---

## Pipeline Architecture

1. **Raw Inputs:** RNA-seq Counts (`2D` / `3D`) + Proteomics Intensities (`2D` / `3D`) + Clinical Metadata
2. **Pre-processing & Quality Control:**
   * *Transcriptomics:* `filterByExpr` (`edgeR`) -> `rlog` (`DESeq2`) -> Top 10k MAD -> Z-score
   * *Proteomics:* Group NA Filter -> `log2(x+1)` -> BPCA Imputation -> Z-score
3. **Unsupervised Exploratory & Latent Factor Analysis:**
   * Preliminary PCA (4 individual omics blocks)
   * MOFA2 (4 latent factors: Full cohort vs. Outlier `RR` exclusion)
4. **Supervised Modeling & Validation (`mixOmics`):**
   * *Single-Omics (`sPLS-DA`):* Dual-criterion LOO-CV (`BER`/`MCC`) + Stability + Permutation Test
   * *Multi-Omics (`DIABLO`):* 2-Phase Grid Search (Cross-Omics `2D`/`3D` & Cross-Dimension `T`/`P`)
5. **Biological Interpretation & Export:**
   * Gene Ontology Over-Representation Analysis (`GOEA`: `BP`, `CC`, `MF`)
   * Multi-Strategy GSEA (`PLS-DA` Loadings, `COL3A1` Pearson, `WGCNA` Network Screening)
   * Automated Exports: High-res PNGs, Formatted Excel Workbooks (`.xlsx`), LaTeX Tables (`.tex`)

---

## Prerequisites and Installation

### R Environment and Packages
The pipeline requires **R (>= 4.2)** and the following CRAN and Bioconductor libraries:

* **Data Wrangling & I/O:** `readr`, `readxl`, `tidyr`, `dplyr`, `stringr`, `purrr`, `writexl`, `openxlsx`, `rstudioapi`
* **Omics & Statistical Modeling:** `DESeq2`, `edgeR`, `pcaMethods`, `sva`, `mixOmics`, `WGCNA`, `MOFA2`, `reticulate`
* **Functional Enrichment & Annotation:** `clusterProfiler`, `org.Hs.eg.db`, `msigdbr`, `fgsea`, `enrichplot`
* **Visualization & Formatting:** `ggplot2`, `ggrepel`, `ggVennDiagram`, `plotly`, `patchwork`, `UpSetR`, `RColorBrewer`, `GGally`, `kableExtra`

**CRAN Installation Command:**
`install.packages(c("readr", "readxl", "tidyr", "dplyr", "ggplot2", "ggrepel", "reticulate", "ggVennDiagram", "plotly", "patchwork", "UpSetR", "stringr", "writexl", "openxlsx", "RColorBrewer", "GGally", "purrr", "kableExtra", "rstudioapi"))`

**Bioconductor Installation Command:**
`BiocManager::install(c("DESeq2", "edgeR", "pcaMethods", "sva", "mixOmics", "WGCNA", "MOFA2", "clusterProfiler", "org.Hs.eg.db", "msigdbr", "fgsea", "enrichplot"))`

### Python / Conda Environment for MOFA2
Training **MOFA2** models requires a Python environment (`mofa_env`, Python 3.10) with `mofapy2` installed (`pip install mofapy2`) and linked via `reticulate`.

---

## Input Data Structure

The pipeline aligns **14 biological samples** (10 experimental `EXP1` cases and 4 control `CTR1` subjects) across four molecular layers:

1. **`COL3A1_sampleAttributes_3D.xlsx`**: Sample metadata sheet containing sample identifiers and the experimental group column `Condition` (`CTR1` vs. `EXP1`).
2. **`COL3A1_rawCounts_mod_2D.txt`**: Raw RNA-seq count matrix for 2D cultures indexed by `ensembl_gene_id`.
3. **`COL3A1_rawCounts_mod_3D.tsv`**: Raw RNA-seq count matrix for 3D cultures indexed by `gene_id`.
4. **`proteomics_Original Matrix_Col3A1_2D_3D.xlsx`**: Raw protein intensity matrix indexed by `PG.Genes`, containing 2D samples in columns `1:14` and 3D samples in columns `15:28`.
5. *(Optional)* **Literature Reference Tables:** External GSEA benchmark workbooks (`Supplementary_Table_GSEA_2D_Correl.xlsx` and `Supplementary_Table_GSEA_2D_GSE27019.xlsx`) used for pathway overlap validation.

> **Sample Alignment:** At startup, explicit index vectors (`neworder2D`, `neworder3D`, `neworderMeta`) reorder all matrices so that columns match identically across all omics blocks (`19-1597` through `RR`).

---

## Analytical Workflow

### 1. Pre-processing and Normalization
* **Overdispersion Assessment:** Log-log mean-variance plots are generated against the Poisson expectation ($Var = \mu$) to verify overdispersion prior to count normalization.
* **Transcriptomics (2D & 3D):**
  1. Low-count genes are filtered using `edgeR::filterByExpr()` grouped by `Condition`, retaining **10,973 genes (2D)** and **11,376 genes (3D)** out of 33,694 initial features.
  2. Variance-stabilizing transformation is applied via `DESeq2::rlog(..., blind = TRUE)`.
  3. Common genes between 2D and 3D are intersected, ranked by **Median Absolute Deviation (MAD)** to select the **top 10,000 most variable genes**, and Z-score standardized per feature ($\mu = 0, \sigma = 1$).
* **Proteomics (2D & 3D):**
  1. String sanitization removes quotes/whitespace, standardizes decimal separators, strips secondary protein isoforms (`sub(";.*", "")`), and applies $\log_2(x + 1)$ transformation.
  2. Features are retained if detected in `> 6` `EXP1` samples or `> 2` `CTR1` samples.
  3. Missing values (`NA`) are imputed using **Bayesian PCA** (`pcaMethods::pca(..., method = "bpca", nPcs = 5)`), followed by 2D/3D intersection and feature-wise Z-score scaling.

### 2. Unsupervised Multi-Omics Integration (MOFA2)
* Trains **Multi-Omics Factor Analysis (MOFA2)** models with `num_factors = 4`, structuring samples into two groups (`groups = c("2D", "3D")`).
* **Outlier Sensitivity Analysis:** Compares the full 14-patient model (`MOFA_model_COL3A1`) against a model excluding sample **`RR`** (`MOFA_meno_RR_model_COL3A1`), which behaves as an outlier in unsupervised latent space.
* Evaluates variance explained ($R^2$) per view and group, factor orthogonality, sample distributions across `CTR1` and `EXP1`, and Reactome pathway enrichment (`msigdbr` C2:CP:REACTOME) on latent factor weights.

### 3. Supervised Classification & Integration (sPLS-DA & DIABLO)
* **Preliminary PCA:** Evaluates unsupervised sample clustering across 10 principal components for all four scaled datasets.
* **Single-Omics sPLS-DA (`Trascr2D`, `Prot2D`, `Trascr3D`, `Prot3D`):**
  * Optimizes sparsity (`keepX`) via Leave-One-Out Cross-Validation (LOO-CV) using Mahalanobis distance.
  * Balances **Balanced Error Rate (BER)**, **Matthews Correlation Coefficient (MCC)**, False Negative limits (`max_FN`), and LOO feature selection stability ($\ge 80\%$).
  * Assesses statistical significance via empirical/exhaustive **Permutation Tests** (1,001 runs).
* **Multi-Omics DIABLO (`block.splsda`):**
  * Pairwise `spls` confirms strong inter-block correlation ($\rho > 0.98$), justifying integrative modeling with discriminative design weights (`0.01` to `0.4`).
  * **Two-Phase Grid Search:** Coarse parameter sweep (`grid_search_splsda_advanced`) followed by local neighborhood refinement (`prepare_param_second_phase_advanced`).
  * Evaluates four integrative architectures:
    1. *Cross-Omics 2D:* `Transcriptomics 2D` + `Proteomics 2D`
    2. *Cross-Omics 3D:* `Transcriptomics 3D` + `Proteomics 3D`
    3. *Cross-Dimension Transcriptomics:* `Transcriptomics 2D` + `Transcriptomics 3D`
    4. *Cross-Dimension Proteomics:* `Proteomics 2D` + `Proteomics 3D`

### 4. Functional Enrichment Analysis (GOEA & GSEA)
* **Gene Ontology Over-Representation Analysis (GOEA):** Compares enriched `BP`, `CC`, and `MF` terms (`clusterProfiler::compareCluster`) between single-omics and integrated proteomics signatures across multiple Benjamini-Hochberg thresholds ($\alpha \in \{0.05, 0.10, 0.25\}$).
* **Multi-Strategy GSEA (`fgsea` on MSigDB C5 GO & C2 KEGG/Reactome):**
  1. **Multivariate PLS-DA Ranking:** Ranks features by Component 1 loadings in single-omics `plsda` and multi-block `block.plsda` (supporting both block-specific and Z-score consensus fused rankings).
  2. **`COL3A1` Co-expression Ranking:** Ranks all genes/proteins by Pearson correlation against the `COL3A1` expression profile.
  3. **WGCNA Network Screening Ranking:** Computes oriented topological network scores:
     $$GS_i = |\text{cor}(x_i, Y)|^\gamma, \quad a_{ij} = |\text{cor}(x_i, x_j)|^\alpha \ (j \neq i), \quad q_i = GS_i \sum_{j \neq i} a_{ij} GS_j, \quad r(g_i) = \text{sign}(\text{cor}(x_i, Y)) \cdot q_i$$
     using automated soft-threshold power estimation (`pickSoftThreshold`) and Z-score consensus merging across blocks.

### 5. Automated Reporting and LaTeX Export
* Filters significant GSEA pathways ($p_{\text{adj}} < 0.20, |\text{NES}| \ge 1.0$ for Excel reports; $p_{\text{adj}} \le 0.05$ for thesis appendices).
* Generates publication-ready `.tex` tables (`kableExtra` with `longtable` and `booktabs`) listing enriched pathways and selected sPLS-DA / DIABLO biomarkers with their loading weights.

---

## Core Custom Functions

| Function | Purpose |
| :--- | :--- |
| `tune_and_run_dual_criterion()` | Performs LOO-CV grid tuning for single-omics `splsda`, optimizing BER, MCC, error counts (FP/FN), and feature stability. |
| `permutation_test_singola_omica()` / `permutation_test_multi_omica()` | Evaluates exact (when $\binom{n}{n_{\text{pos}}} \le 5000$) or Monte Carlo null distributions to compute component-wise empirical p-values for MCC and BER. |
| `grid_search_splsda_advanced()` | Multi-block grid search for DIABLO evaluating LOO majority-vote classification independently for Component 1 and Component 2. |
| `prepare_param_second_phase_advanced()` | Builds a refined local hyperparameter grid around the best Phase 1 `design` and `keepX` values. |
| `run_splsda_pipeline()` | Fits the final DIABLO model, computes confusion matrix metrics (`TP`, `TN`, `FP`, `FN`, `SENS`, `SPEC`, `BER`, `B_ACC`, `MCC`), extracts biomarkers, and runs LOO stability analysis. |
| `generate_splsda_report_optimized()` | Exports an 8-plot high-resolution diagnostic suite (`0_Correlation_Plot.png` through `8_Circos.png`) highlighting misclassified samples and cross-omics networks. |
| `run_gsea_on_block_plsda()` | Runs `fgsea` on `block.plsda` Component 1 loadings either per individual block or via Z-standardized block fusion (`fuse_blocks = TRUE`). |
| `get_wgcna_network_rank()` / `fondi_ranking_wgcna()` | Calculates weighted topological network screening scores via `WGCNA` and merges Z-standardized rankings across omics/dimensions. |

---

## Generated Outputs

Running the pipeline produces the following files and directories:

1. **Saved MOFA2 Models:**
   * `MOFA_model_COL3A1.hdf5` & `MOFA_model_COL3A1.RData`
   * `MOFA_meno_RR_model_COL3A1.hdf5` & `MOFA_meno_RR_model_COL3A1.RData`
2. **Diagnostic Plot Directories (`Report_<Dataset>/`):**
   * `0_Correlation_Plot.png`, `1_Sample_Separation.png`, `2_Correlation_Circle.png`, `3_Boxplot_*.png`, `4_Loadings_*.png`, `5_Network.png`, `6_CIM_Heatmap.png`, `7_Stability_*.png`, `8_Circos.png`.
   * Comparative multi-omics figures (`Figure_sPLSDA_All_Omics_Misclassified_Filled.png`, `Figure_Boxplot_Scores_Comp1.png`, `Figure_PlotLoadings_Comp1_ggplot.png`, `metriche_*.png`).
   * Enrichment visualizations (`GO_Emapplot_*.png`, `GO_Cnetplot_*.png`, `plot_GSEA/GSEA_*.png`).
3. **Biomarker Gene Lists (`.txt`):**
   * `geni_trascr2D_per_GO.txt`, `prot_prot2D_per_GO.txt`, `geni_per_GO_3D.txt`, `prot_per_GO_3D.txt`, `GO_geni_Trascr_2D.txt`, `GO_prot_Prot_2D.txt`.
4. **Styled Excel Workbooks (`.xlsx`):**
   * `GO_Prot_singola_e_integrazione.xlsx`
   * `Risultati_GSEA_MultiOmica.xlsx` (including `4_Intersezioni_Core` for cross-validated core pathways)
   * `Risultati_GSEA_Correlazione_COL3A1.xlsx`
   * `Risultati_GSEA_WGCNA.xlsx`
5. **LaTeX Appendix Tables (`Tabelle_LaTeX/*.tex`):**
   * `Biomarcatori_Proteomica_2D_3D.tex`
   * `Biomarcatori_DIABLO_Prot2D_Prot3D.tex`
   * Significant GSEA pathway tables formatted per workbook and sheet.

---

## Reproducibility Notes

Before running the script on a different system, update the local paths configured in the script:
* **Conda Environment:** `use_condaenv("/Users/.../envs/mofa_env", required = TRUE)`
* **Input Files:** Paths passed to `read_excel()` and `read_table()` for transcriptomics, proteomics, and reference benchmark tables.
* **Output Directories:** Variables `cartella_output`, `dir_plot_gsea`, and `percorso_cartella` for Excel, PNG, and LaTeX exports.
