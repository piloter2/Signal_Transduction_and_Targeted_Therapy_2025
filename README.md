### Supplementary Figure 2.  Single-Cell RNA-seq Analysis Pipeline for Oligodendrocyte Lineage

This repository contains the R analysis pipeline for integrating and analyzing single-cell RNA sequencing (scRNA-seq) data from multiple sources to study oligodendrocyte lineage cells in various mouse models (Alzheimer's Disease, Cuprizone model, etc.).

#### Datasets Integrated
The following datasets are processed and integrated:
- **GSE160512**: PS2APP mouse model.
- **GSE153895**: P301L mouse model.
- **GSE148676**: Cuprizone-induced demyelination model.


### Figure 2 & Supplementary Figure 6

This repository contains the full computational pipeline for analyzing single-cell RNA-sequencing (scRNA-seq) data, focusing on Oligodendrocyte lineage cells in 5xFAD vs. WT mouse models.

#### Pipeline Overview

The analysis follows a standard single-cell workflow with specific customizations for doublet detection and pathway enrichment.

1.  **Data Loading & QC**: Loads Indrop count matrices and removes doublets using `scDblFinder` (Parallelized).
2.  **Integration**: Merges samples and corrects batch effects using `Harmony` and `SCTransform`.
3.  **Clustering & Annotation**: Unsupervised clustering and manual annotation of cell types (OPC, NFOL, MOL, etc.).
4.  **Scoring**: Calculates signature scores for Inflammation, Myelination, and NFkB using `UCell`.
5.  **Differential Expression**: Identifies DEGs using the `MAST` test.
6.  **Enrichment Analysis**: Performs GSEA for GO terms and Reactome pathways.

#### Prerequisites

#### R Environment
This pipeline uses **R version 4.x**. Ensure the following packages are installed:

```r
* **CRAN**: `Seurat`, `tidyverse`, `Matrix`, `harmony`, `parallel`, `doParallel`
* **Bioconductor**: `SingleCellExperiment`, `scDblFinder`, `clusterProfiler`, `ReactomePA`, `EnhancedVolcano`, `UCell`
* **Custom/GitHub**: `indRop`, `Nebulosa`
install.packages(c("Seurat", "tidyverse", "Matrix", "cowplot", "ggpubr", "harmony", "BiocManager"))
BiocManager::install(c("biomaRt", "SingleCellExperiment", "scDblFinder"))
# Note: 'indRop' is a custom package required for reading Indrop data.
