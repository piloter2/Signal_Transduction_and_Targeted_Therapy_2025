# ==============================================================================
# Title: Single-Cell RNA-seq Analysis of Oligodendrocyte Lineage in 5xFAD Mice
# Description: Integration and analysis of in-house and public scRNA-seq datasets 
#              to investigate oligodendrocyte changes in Alzheimer's disease models.
# Author: Jaemyung Jang / Korea Brain Research Institute
# Date: 2025-12-07
# ==============================================================================

# ==============================================================================
# 1. Setup and Libraries
# ==============================================================================

# Suppress warnings
options(warn = -1)

# Increase memory limits for large datasets
# options(future.globals.maxSize = 1000 * 1024^2) # 1GB
unix::rlimit_as(4 * 1024^4) # Uncomment if using 'unix' package for resource limits

# Python configuration for UMAP (if needed)
# library(reticulate)
# use_python("/usr/bin/python3") 
# py_install("numpy")
# py_install("setuptools")
# py_install("umap-learn")
# py_install("scikit-learn")
# py_install("pandas")
# py_install("numexpr")
# py_install("louvain")
# py_install("leidenalg")
# py_install("MAST")
# py_config()
# py_discover_config()
# cat(names(reticulate::import("umap")), sep="\n")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(harmony)
  library(readxl)
  library(stringr)
  library(tibble)
  library(ggsci)
  library(ggpubr)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(UCell)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ReactomePA)
  library(DOSE)
  library(biomaRt)
  library(openxlsx)
  # library(indRop) # Custom or specific package for reading inDrop data
})

# Load custom functions
source("./functions/single_cell_function_mouse.r")

# Define Base Directory (Adjust this for your environment)
BASE_DIR <- "/path/to/" 
OUTPUT_DIR <- file.path(BASE_DIR, "OL_sorted", "results")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Set seed for reproducibility
set.seed(1234)

# ==============================================================================
# 2. Data Loading and Preprocessing Functions
# ==============================================================================

# Function to read 10X data from a directory (Wrapper for Read10X or custom loading)
# Note: You were using 'readIndrop' and 'Read10X_GEO'. 
# Ensure these functions are defined or loaded.
# Here is a generic placeholder structure.

load_seurat_object <- function(data_path, sample_name, project_name, min_cells = 0) {
  # Logic to read data based on format (mtx, txt, etc.)
  # This part relies on your specific 'readIndrop' or 'Read10X' usage
  # For example:
  # counts <- Read10X(data_path) 
  # or
  # counts <- readIndrop(data_path)
  
  # For this refactor, assuming counts are loaded into 'counts_matrix'
  # seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = project_name, min.cells = min_cells)
  # return(seurat_obj)
  return(NULL) # Placeholder
}

# Function for basic QC and Filtering
preprocess_seurat <- function(obj, mt_pattern = "^mt\\.", hb_pattern = "^Hb[ab]") {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = hb_pattern)
  
  # Add filtering logic here if needed (e.g., subset based on features/counts)
  
  return(obj)
}

# ==============================================================================
# 3. Loading Datasets
# ==============================================================================

# --- 3.1 Load Public Datasets - GSE278199
cat("Loading GSE278199...\n")
gse278199_path <- file.path(BASE_DIR, "/data/GSE278199/GSE278199_OPC_cells.rds")
if (file.exists(gse278199_path)) {
  counts_278199 <- readRDS(gse278199_path)
  rownames(counts_278199) <- gsub("-", "\\.", rownames(counts_278199))
  
  GSE278199.obj <- CreateSeuratObject(counts = counts_278199, min.cells = 0)
  GSE278199.obj$basic <- "GSE278199"
  
  # Assign Metadata (simplify logic)
  sample_ids <- str_sub(colnames(GSE278199.obj), 1, 1)
  GSE278199.obj$id <- case_when(
    sample_ids == names(table(sample_ids))[1] ~ "naïve",
    sample_ids == names(table(sample_ids))[2] ~ "saline_CPZ3",
    sample_ids == names(table(sample_ids))[3] ~ "XPro1595_CPZ3",
    TRUE ~ "Unknown"
  )
  GSE278199.obj <- preprocess_seurat(GSE278199.obj)
}

# --- 3.2 Load Public Datasets - GSE160512

path_GSE160512 <- file.path(BASE_DIR, "data","GSE160512","raw")
fileListCNTs <- list.files(path_GSE160512, pattern = ".txt.gz")
path_in_counts.raw <- paste0(path_GSE160512,"/",fileListCNTs)
names(fileListCNTs) <- path_in_counts.raw

GSE160512.meta <- read.table(paste0(BASE_DIR,"/data/GSE160512/GSE160512_PS2APPTsneCoordinatesAndClusterAssignments.txt.gz"), header =T, sep="\t")
GSE160512.meta <- GSE160512.meta %>% tibble::column_to_rownames('cellID')

require(biomaRt)
annot <- c()
GSE160512.obj <- list()

ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
for(x in path_in_counts.raw){
    require(indRop)
    temp <- readIndrop(x)
    annot <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'),filters = 'ensembl_gene_id',values = colnames(temp),mart = ensembl)
    annot <- annot[-c(which((annot$mgi_symbol == '') | duplicated(annot$mgi_symbol) | duplicated(annot$ensembl_gene_id))),]
    temp <- temp[,annot$ensembl_gene_id]
    colnames(temp) <- gsub("-","\\.", annot$mgi_symbol)
    GSE160512.obj[[fileListCNTs[x]]] <- CreateSeuratObject(counts = t(temp), min.cells = 0)               
}

for(x in path_in_counts.raw){
    y <- fileListCNTs[x]
    GSE160512.obj[[y]]@meta.data$id <- fileListCNTs[x]
    GSE160512.obj[[y]]@meta.data$basic <- "GSE160512"
    GSE160512.obj[[y]]@meta.data$percent.mt <- PercentageFeatureSet(GSE160512.obj[[y]], pattern = paste0(c(mt.genes,"^mt."), collapse = "|"))
    GSE160512.obj[[y]]@meta.data$percent.hb <- PercentageFeatureSet(GSE160512.obj[[y]], pattern = paste0(c('Hbq1a','Hbb.y','Hbb.bh1','Hbb.bs','Hba.x','Hba.a2','Hba.a1','Hbq1b','Hbb.bt','Hbb.bh2','Hbb.bh3','Hba.ps4','Hbb.bh0'), collapse = "|"))

    GSE160512.obj[[y]]@meta.data$celltype <- GSE160512.meta[colnames(GSE160512.obj[[y]]),]$allCells.cluster.interpretation
    GSE160512.obj[[y]]@meta.data$type <- GSE160512.meta[colnames(GSE160512.obj[[y]]),]$genotype
}

GSE160512.oligo <- list()

for(y in seq_len(length(GSE160512.obj))){
    GSE160512.oligo[[y]] <- GSE160512.obj[[y]][, grep("^oligo|^OPC|^COP", GSE160512.obj[[y]]@meta.data$celltype)]
}

# --- 3.3 Load Public Datasets -  GSE153895

path_GSE153895 <- file.path(BASE_DIR, "data" ,"reference","GSE153895","raw")
fileListCNTs <- list.files(path_GSE153895, pattern = ".txt.gz")
path_in_counts.raw <- paste0(path_GSE153895,"/",fileListCNTs)
names(fileListCNTs) <- path_in_counts.raw

GSE153895.meta <- read.table(paste0(BASE_DIR,"/data//GSE153895/GSE153895_TsneCoordinatesAndClusterAssignments.txt.gz"), header =T, sep="\t")
GSE153895.meta <- GSE153895.meta %>% tibble::column_to_rownames('cellID')

require(biomaRt)
annot <- c()
GSE153895.obj <- list()
ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

for(x in path_in_counts.raw){
    require(indRop)
    temp <- readIndrop(x)
    annot <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'),filters = 'ensembl_gene_id',values = colnames(temp),mart = ensembl)
    annot <- annot[-c(which((annot$mgi_symbol == '') | duplicated(annot$mgi_symbol) | duplicated(annot$ensembl_gene_id))),]
    temp <- temp[,annot$ensembl_gene_id]
    colnames(temp) <- gsub("-","\\.", annot$mgi_symbol)

    GSE153895.obj[[fileListCNTs[x]]] <- CreateSeuratObject(counts = t(temp), min.cells = 0)               
}


for(x in path_in_counts.raw){
    y <- fileListCNTs[x]
    GSE153895.obj[[y]]@meta.data$celltype <- GSE153895.meta[colnames(GSE153895.obj[[y]]),]$allCells.cluster.interpretation
    GSE153895.obj[[y]]@meta.data$type <- GSE153895.meta[colnames(GSE153895.obj[[y]]),]$genotype

    GSE153895.obj[[y]]@meta.data$id <- fileListCNTs[x]
    GSE153895.obj[[y]]@meta.data$basic <- "GSE153895"
    GSE153895.obj[[y]]@meta.data$percent.mt <- PercentageFeatureSet(GSE153895.obj[[y]], pattern = paste0(c(mt.genes,"^mt."), collapse = "|"))
    GSE153895.obj[[y]]@meta.data$percent.hb <- PercentageFeatureSet(GSE153895.obj[[y]], pattern = paste0(c('Hbq1a','Hbb.y','Hbb.bh1','Hbb.bs','Hba.x','Hba.a2','Hba.a1','Hbq1b','Hbb.bt','Hbb.bh2','Hbb.bh3','Hba.ps4','Hbb.bh0'), collapse = "|"))
    # GSE153895.obj.OL[[y]] <- GSE153895.obj[[y]][, grep("^oligo|^OPC|^COP|C15", GSE153895.obj[[y]]$GSE153895.celltype)]

}


GSE153895.oligo <- list()

for(y in seq_len(length(GSE153895.obj))){
    GSE153895.oligo[[y]] <- GSE153895.obj[[y]][, grep("^oligo|^OPC|^COP|C15", GSE153895.obj[[y]]@meta.data$celltype)]
}


# ==============================================================================
# 4. Integration and Dimensionality Reduction
# ==============================================================================

# Merge Datasets
OL.merge <- merge(GSE278199.obj, y = c(GSE153895.oligo, GSE160512.oligo), project = 'merged')

# Normalization and Regression
OL.merge <- OL.merge %>%
                  SCTransform(assay = 'RNA',      
                              new.assay.name = 'SCT',  
                              vars.to.regress = c("percent.mt","percent.hb","nFeature_RNA"),
                              verbose = T) 
# Integration Features

OL.features <- SelectIntegrationFeatures(object.list = list(GSE278199.obj, GSE153895.oligo, GSE160512.oligo), nfeatures = 3000)
OL.feature <- OL.features[-grep("^mt.|Rik$|^Hba.|^Hbb.",OL.features)]

# Dimensionality Reduction
pcs <- 1:25
OL.merged <- OL.merge %>%
                          RunPCA(verbose = FALSE, assay = "SCT", features = OL.feature ) %>%
                          RunHarmony(group.by=c("basic"), assay.use = "SCT", dims=pcs, max.iter.harmony = 100, max.iter.cluster=250) %>% 
                          RunUMAP(reduction = "harmony", umap.method = "umap-learn", assay = "SCT",  dims=pcs) %>%
                          FindNeighbors(reduction = "harmony", dims=pcs, assay = "SCT") %>%
                          FindClusters(resolution = seq(0.1,1.6,0.1), algorithm = 2)

# ==============================================================================
# 5. Cell Type Annotation and Markers
# ==============================================================================

# Set Identity
Idents(OL.merged) <- "SCT_snn_res.0.4"

# Visualization
p_umap <- DimPlot(OL.merged, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
print(p_umap)

# DotPlot for Marker Genes
markers_to_plot <- c("Pdgfra", "Cspg4", "Olig1", "Olig2", # OPC
                     "Fyn", "Sirt2", "Enpp6", "Bcas1",    # COP/NFOL
                     "Plp1", "Mbp", "Mog", "Mag",         # MOL
                     "Cnp", "Mal") 

DotPlot(OL.merged, features = markers_to_plot, assay = "SCT", cols = "RdBu") + 
  RotatedAxis() + coord_flip()
