# Multi-Omics Integrative Analysis Pipeline (`COL3A1` 2D vs 3D Models)

![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3?logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-Omics-87B13F)
![License](https://img.shields.io/badge/License-MIT-blue.svg)

End-to-end bioinformatics pipeline written in **R** for exploratory analysis, predictive modeling, and multi-omics data integration (**RNA-seq Transcriptomics** and **Proteomics**) across **2D** and **3D** cellular models investigating **`COL3A1`**.

## Overview

The workflow combines unsupervised latent factor modeling with supervised discriminant feature selection and multi-level functional enrichment:

* **Data Pre-processing & QC:** Overdispersion diagnostics, `edgeR` expression filtering, `DESeq2` regularized log (`rlog`) transformation, MAD feature selection, $\log_2$ transformation, and Bayesian PCA (`BPCA`) missing value imputation.
* **Unsupervised Integration (`MOFA2`):** Latent factor extraction across 2D and 3D sample groups, sensitivity analysis (full cohort vs. outlier `RR` exclusion), and Reactome pathway enrichment on factor weights.
* **Supervised Classification (`mixOmics`):**
  * *Single-Omics (`sPLS-DA`):* Dual-criterion Leave-One-Out Cross-Validation (LOO-CV) tuning balancing Balanced Error Rate (BER), Matthews Correlation Coefficient (MCC), False Negative constraints, and feature stability ($\ge 80\%$), validated via empirical/exhaustive permutation tests.
  * *Multi-Omics Integration (`DIABLO` / `block.splsda`):* Two-phase grid search across cross-omics (`Transcriptomics + Proteomics`) and cross-dimensional (`2D + 3D`) designs.
* **Functional Enrichment (`GOEA` & `GSEA`):** Over-Representation Analysis (`clusterProfiler`) and `fgsea` across three ranking strategies: multivariate PLS-DA loadings, Pearson correlation with `COL3A1`, and oriented **WGCNA Network Screening**.
* **Automated Reporting:** High-resolution diagnostic plots, multi-sheet formatted Excel workbooks (`.xlsx`), and publication-ready LaTeX tables (`.tex`).

## Quick Start

1. **Install dependencies:** Ensure R ($\ge 4.2$) and a Conda environment (`mofa_env` with `mofapy2`) are installed. See full instructions in [`DOCUMENTATION.md`](DOCUMENTATION.md#prerequisites-and-installation).
2. **Configure file paths:** Update the input/output directory paths and Conda environment path at the top of the main R script.
3. **Run the pipeline:** Execute the script sequentially from RStudio or via `Rscript`.

## Full Documentation

For detailed methodology, mathematical formulations, function references, input schemas, and output descriptions, please refer to the full technical manual:

👉 **[Read the Complete Technical Documentation (`DOCUMENTATION.md`)](DOCUMENTATION.md)**
