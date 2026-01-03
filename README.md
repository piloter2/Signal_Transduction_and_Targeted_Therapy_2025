## Oligodendrocyte precursor cells-microglia cross talk via BMP4 signaling drives microglial neuroprotective response and mitigates Alzheimer's disease progression.
This repository contains the full computational pipeline for analyzing single-cell RNA-sequencing (scRNA-seq) data, focusing on Oligodendrocyte lineage cells in 5xFAD vs. WT mouse models.

### 📌 Project Overview
#### Study Title: OPC-Microglia Crosstalk via BMP4 Signaling in Alzheimer's Disease
- Key Biological Focus: BMP4 signaling, OPCs, Microglia, Alzheimer's Disease (5xFAD model).
- Author: Jaemyung Jang (Korea Brain Research Institute)
- Contact: piloter2@kbri.re.kr

#### 📂 Directory Structure
- The code assumes the following directory structure (paths must be configured in the script):
/path/to/base_dir/
├── data/                  # Raw counts (.counts.tsv) and metrics (.metrics.tsv)
├── scripts/               # R scripts and custom functions
├── results/               # Output tables and figures
└── final_Rdata/           # Processed Seurat objects (.RDS)

### Supplementary Figure 2.  Single-Cell RNA-seq Analysis Pipeline for Oligodendrocyte Lineage
This repository contains the R analysis pipeline for integrating and analyzing single-cell RNA sequencing (scRNA-seq) data from multiple sources to study oligodendrocyte lineage cells in various mouse models (Alzheimer's Disease, Cuprizone model, etc.).

#### Datasets Integrated
The following datasets are processed and integrated:
- **GSE160512**: PS2APP mouse model.
- **GSE153895**: P301L mouse model.
- **GSE148676**: Cuprizone-induced demyelination model.

### Figure 2 & Supplementary Figure 6, 8
This repository contains the R analysis pipeline used to investigate the role of BMP4 signaling in Oligodendrocyte Precursor Cells (OPCs) and their crosstalk with Microglia in Alzheimer's disease (AD) progression. The analysis integrates in-house scRNA-seq datasets (inDrops platform) with public datasets (SmartSeq2) to identify neuroprotective responses.

### 🚀 Analysis Workflow
The pipeline consists of two main stages:

##### 1. Preprocessing & Integration
- Data Loading: Parallel loading of raw counts and metrics from the inDrops platform.
- Quality Control (QC):
  - Filtering based on UMI counts (UMIFM >= 1000).
  - Doublet Removal: Using scDblFinder.
  - Ambient RNA Correction: Using SoupX to remove background contamination.
  - Mitochondrial and Hemoglobin gene filtering.
- Integration:
  - Merging in-house data with public dataset GSE71585 (Reference).
  - Normalization via SCTransform (regressing out cell cycle, mito %, etc.).
- Batch correction using Harmony.
- Clustering and Silhouette score assessment.

##### 2. Downstream Analysis & Visualization
- Annotation: Mapping clusters to biological identities (e.g., Microglia, OPC, COP, MOL) using a predefined cell_type_map.
- Visualization:
  - DimPlot: UMAP visualization of cell types.
  - DotPlot: Expression of marker genes (e.g., Pdgfa, Slc1a2, Snap25).
  - Differential Expression (DE):
  - Identifying DEGs between genotypes (e.g., BMP4-KO vs. 5xFAD Control) using the MAST test.
- Comparison performed per cell type.
- Enrichment Analysis (GSEA):
  - Gene Ontology (GO): BP, CC, MF.
  - Reactome Pathway Analysis.
  - WikiPathways.

##### Enviroments
''r 
Core R Packages
Seurat (v4.0+)
Harmony (Integration)
scDblFinder (Doublet removal)
SoupX (Ambient RNA correction)
UCell (Scoring)
Tidyverse (dplyr, ggplot2, etc.)
Bio-conductor Packages
SingleCellExperiment
clusterProfiler, ReactomePA, DOSE
org.Mm.eg.db, AnnotationDbi
Custom Packages
indRop (In-house package for processing inDrops data)

