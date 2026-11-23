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

```text
[Raw Inputs: RNA-seq Counts (2D/3D) + Proteomics Intensities (2D/3D) + Clinical Metadata]
       │
       ├──► 1. Pre-processing & Quality Control
       │       ├── Transcriptomics: filterByExpr (edgeR) -> rlog (DESeq2) -> Top 10k MAD -> Z-score
       │       └── Proteomics: Group NA Filter -> log2(x+1) -> BPCA Imputation -> Z-score
       │
       ├──► 2. Unsupervised Exploratory & Latent Factor Analysis
       │       ├── Preliminary PCA (4 individual omics blocks)
       │       └── MOFA2 (4 latent factors: Full cohort vs. Outlier "RR" exclusion)
       │
       ├──► 3. Supervised Modeling & Validation (mixOmics)
       │       ├── Single-Omics (sPLS-DA): Dual-criterion LOO-CV (BER/MCC) + Stability + Permutation Test
       │       └── Multi-Omics (DIABLO): 2-Phase Grid Search (Cross-Omics 2D/3D & Cross-Dimension T/P)
       │
       └──► 4. Biological Interpretation & Export
               ├── Gene Ontology Over-Representation Analysis (GOEA: BP, CC, MF)
               ├── Multi-Strategy GSEA (PLS-DA Loadings, COL3A1 Pearson, WGCNA Network Screening)
               └── Automated Exports: High-res PNGs, Formatted Excel Workbooks (.xlsx), LaTeX Tables (.tex)
