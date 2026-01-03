## Oligodendrocyte Lineage Single-Cell Integration Analysis

This repository contains the R analysis pipeline for integrating and analyzing single-cell RNA sequencing (scRNA-seq) data from multiple sources to study oligodendrocyte lineage cells in various mouse models (Alzheimer's Disease, Cuprizone model, etc.).

#### Datasets Integrated
The following datasets are processed and integrated:
- **GSE160512**: PS2APP mouse model.
- **GSE153895**: P301L mouse model.
- **GSE148676**: Cuprizone-induced demyelination model.

#### Prerequisites

#### R Environment
This pipeline uses **R version 4.x**. Ensure the following packages are installed:

```r
install.packages(c("Seurat", "tidyverse", "Matrix", "cowplot", "ggpubr", "harmony", "BiocManager"))
BiocManager::install(c("biomaRt", "SingleCellExperiment", "scDblFinder"))
# Note: 'indRop' is a custom package required for reading Indrop data.
