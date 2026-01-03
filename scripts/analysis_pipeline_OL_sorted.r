# ==============================================================================
# Title : Oligodendrocyte precusor cells-microglia cross talk via BMP4 signaling drives microglial neuroprotective response and mitigates Alzheimer's diseas progression
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
options(warn=-1)

# Base Directories (USER MUST UPDATE THESE PATHS)
BASE_DIR <- "/path/to/data" 
OUT_DIR  <- "/path/to/results"

# Define Marker Genes
OL_MARKERS <- c(
  "Pdgfra","Vcan","Olig1","Ptprz1","Olig2", # OPC
  "Plp1","Cldn11","Cnp","Tspan2","Mobp",    # Oligodendrocyte
  "Trf","Apod","Ermn","Car2",
  "Ptgds","Plp1","Mal","Cryab","Car2", 
  "Fyn","Sirt2","Tuba1a","Bcas1","Enpp6",   # COP
  "Hexb","Ctss","C1qb","C1qc","Cx3cr1",     # Microglia
  "Lyz2","Mrc1","Pf4","Ctsb","F13a1",       # Macrophage
  "Ly6c1","Ly6a","Flt1","Itm2a","Slco1a4",  # Endothelial
  "Rgs5","Acta2","Myl9","Tagln","Cald1",    # Mural
  "Dcn","Col1a2","Igf2","Igfbp2","Pcolce",  # VLMC
  "Ccdc153","Rarres2","Tmem212","Nnat","Rsph1", # Ependyma    
  "Cpe","Apoe","Cst3","Mt1","Dbi","Aldoc","Gja1","Slc1a2","Slc1a3","Clu", # Astrocyte
  "Ccnd2","Tubb5","Dlx1","Tmsb10","Tuba1a", # NSC
  "Hmgb2","Top2a","Hmgn2","Tubb5","H2afz",  # NSC    
  "Sst","Npy","Resp18","Nap1l5","Nos1",     # Neuron
  "Lypd1","Meg3","Atp1b1","Snhg11","Olfm1", # Neuron
  "Penk","Arpp21","Ppp3ca","Atp2b1","Hpca"  # Neuron
)

# Mitochondrial & Hemoglobin Patterns
mt.genes.pattern <- "^mt."
hb.genes.pattern <- paste0(c('Hbp1','Hbq1a','Hbb.y','Hbb.bh1','Hbb.bs','Hba.x','Hba.a2','Hba.a1','Hbq1b','Hbb.bt','Hbb.bh2','Hbb.bh3','Hba.ps4','Hbb.bh0'), collapse = "|")

# Initialize Ensembl Mart
ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

# Helper Function for Gene Conversion and Object Creation
create_seurat_from_indrop <- function(file_path, mart_obj) {
  temp <- readIndrop(file_path)
  
  # Convert Ensembl IDs to Gene Symbols
  annot <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'),
                 filters = 'ensembl_gene_id',
                 values = colnames(temp),
                 mart = mart_obj)
  
  # Remove duplicates and empty symbols
  annot <- annot[-c(which((annot$mgi_symbol == '') | duplicated(annot$mgi_symbol) | duplicated(annot$ensembl_gene_id))),]
  
  temp <- temp[, annot$ensembl_gene_id]
  colnames(temp) <- gsub("-", "\\.", annot$mgi_symbol)
  
  sobj <- CreateSeuratObject(counts = t(temp), min.cells = 0)
  return(sobj)
}

# ------------------------------------------------------------------------------
# 2. Data Loading & Preprocessing
# ------------------------------------------------------------------------------

### 2.1 Load GSE160512
cat("Processing GSE160512...\n")
path_gse160512 <- file.path(BASE_DIR, "GSE160512", "raw")
files_160512 <- list.files(path_gse160512, pattern = ".txt.gz", full.names = TRUE)
meta_160512 <- read.table(file.path(BASE_DIR, "GSE160512", "GSE160512_PS2APPTsneCoordinatesAndClusterAssignments.txt.gz"), header =T, sep="\t") %>% 
  tibble::column_to_rownames('cellID')

obj_list_160512 <- list()
for(f in files_160512) {
  fname <- basename(f)
  sobj <- create_seurat_from_indrop(f, ensembl)
  
  sobj$id <- fname
  sobj$basic <- "GSE160512"
  sobj$percent.mt <- PercentageFeatureSet(sobj, pattern = paste0(c(mt.genes.pattern), collapse = "|"))
  sobj$percent.hb <- PercentageFeatureSet(sobj, pattern = hb.genes.pattern)
  
  # Metadata mapping
  common_cells <- intersect(colnames(sobj), rownames(meta_160512))
  sobj <- subset(sobj, cells = common_cells) # Ensure matching cells
  sobj$GSE160512.celltype <- meta_160512[colnames(sobj),]$allCells.cluster.interpretation
  sobj$GSE160512.genotype <- meta_160512[colnames(sobj),]$genotype
  
  obj_list_160512[[fname]] <- sobj
}

# 2.2 Load GSE153895
cat("Processing GSE153895...\n")
path_gse153895 <- file.path(BASE_DIR, "GSE153895", "raw")
files_153895 <- list.files(path_gse153895, pattern = ".txt.gz", full.names = TRUE)
meta_153895 <- read.table(file.path(BASE_DIR, "GSE153895", "GSE153895_TsneCoordinatesAndClusterAssignments.txt.gz"), header =T, sep="\t") %>% 
  tibble::column_to_rownames('cellID')

obj_list_153895 <- list()
for(f in files_153895) {
  fname <- basename(f)
  sobj <- create_seurat_from_indrop(f, ensembl)
  
  sobj$id <- fname
  sobj$basic <- "GSE153895"
  sobj$percent.mt <- PercentageFeatureSet(sobj, pattern = paste0(c(mt.genes.pattern), collapse = "|"))
  sobj$percent.hb <- PercentageFeatureSet(sobj, pattern = hb.genes.pattern)
  
  common_cells <- intersect(colnames(sobj), rownames(meta_153895))
  sobj <- subset(sobj, cells = common_cells)
  sobj$GSE153895.celltype <- meta_153895[colnames(sobj),]$allCells.cluster.interpretation
  sobj$GSE153895.genotype <- meta_153895[colnames(sobj),]$genotype
  
  obj_list_153895[[fname]] <- sobj
}
# Select specific subsets as per original code
obj_list_153895_sub <- obj_list_153895[c(6,7,9,11,12)] 

# 2.3 Load GSE148676
cat("Processing GSE148676...\n")
path_gse148676 <- file.path(BASE_DIR, "GSE148676")
files_148676 <- list.files(path_gse148676, pattern = ".txt.gz", full.names = TRUE)

# Metadata for GSE148676
type_map <- data.frame(
  type = c("WT_Baseline_rep1", "WT_Baseline_rep2", "WT_Cuprizone_4w_rep1", 
           "WT_Cuprizone_4w_rep2", "WT_Baseline_rep3", "WT_Cuprizone_4w_rep3", "WT_Cuprizone_4w_rep4"),
  row.names = c("GSM4476524", "GSM4476526","GSM4476528","GSM4476532","GSM4476536", "GSM4476538", "GSM4476540")
)

obj_list_148676 <- list()
for(f in files_148676) {
  fname <- basename(f)
  sobj <- create_seurat_from_indrop(f, ensembl)
  
  sobj$id <- fname
  sobj$basic <- "GSE148676"
  gsm_id <- strsplit(fname, "_")[[1]][1]
  sobj$type <- type_map[gsm_id, "type"]
  sobj$percent.mt <- PercentageFeatureSet(sobj, pattern = paste0(c(mt.genes.pattern), collapse = "|"))
  sobj$percent.hb <- PercentageFeatureSet(sobj, pattern = hb.genes.pattern)
  
  obj_list_148676[[fname]] <- sobj
}

# ------------------------------------------------------------------------------
# 3. Integration & Normalization
# ------------------------------------------------------------------------------
cat("Merging datasets...\n")
ref.merge <- merge(
  x = obj_list_160512[[1]], 
  y = c(obj_list_160512[-1], obj_list_153895_sub, obj_list_148676), 
  project = 'merged', 
  merge.data = TRUE
)

cat("Running SCTransform...\n")
ref.merge <- SCTransform(ref.merge, assay = 'RNA', new.assay.name = 'SCT', 
                         vars.to.regress = c("percent.mt", "nFeature_RNA"), verbose = TRUE)
ref.merge <- FindVariableFeatures(ref.merge, assay = "SCT", nfeatures = 2000)

cat("Running PCA & Harmony...\n")
pcs <- 1:30
DefaultAssay(ref.merge) <- "SCT"
ref.merge <- RunPCA(ref.merge, verbose = FALSE, assay = "SCT", features = VariableFeatures(ref.merge))
ref.merge <- RunHarmony(ref.merge, group.by.vars = "id", assay.use = "SCT", dims = pcs)

cat("Clustering & UMAP...\n")
ref.merge <- RunUMAP(ref.merge, reduction = "harmony", umap.method = "umap-learn", assay = "SCT", dims = pcs)
ref.merge <- FindNeighbors(ref.merge, reduction = "harmony", dims = pcs, assay = "SCT")
ref.merge <- FindClusters(ref.merge, resolution = seq(0.1, 1, 0.1), algorithm = 2)

# ------------------------------------------------------------------------------
# 4. Visualization
# ------------------------------------------------------------------------------
options(repr.plot.width = 25, repr.plot.height = 35)

p1 <- list()
p1[[1]] <- DimPlot(ref.merge, reduction = "umap", group.by = "id", label = T) + NoLegend()
p1[[2]] <- DimPlot(ref.merge, group.by = "id", label = T) + NoLegend()
p1[[3]] <- DimPlot(ref.merge, group.by = "SCT_snn_res.0.5", label = T)
p1[[4]] <- DimPlot(ref.merge, group.by = "basic", label = T)
p1[[5]] <- DotPlot(ref.merge, group.by = "SCT_snn_res.0.5", features = unique(OL_MARKERS), assay="SCT") + coord_flip()

# Save plots
# ggsave(file.path(OUT_DIR, "QC_Plots.pdf"), plot = plot_grid(plotlist = p1, ncol=2), width=25, height=35)

# ------------------------------------------------------------------------------
# 5. Annotation & Subsetting
# ------------------------------------------------------------------------------
# NOTE: Cluster IDs ("0", "1", etc.) may change between runs. Verify before assigning labels!
celltype_map <- c(
  "1" = "Oligodendrocyte", "3" = "OPC", "4" = "Oligodendrocyte", 
  "5" = "Oligodendrocyte", "6" = "Oligodendrocyte", "14" = "Oligodendrocyte", 
  "15" = "Oligodendrocyte", "21" = "Oligodendrocyte", "24" = "COP", 
  "25" = "Oligodendrocyte", "26" = "OPC", "32" = "OPC", "33" = "Oligodendrocyte"
)

Idents(ref.merge) <- "SCT_snn_res.0.5"
ref.merge <- RenameIdents(ref.merge, celltype_map)
ref.merge$ref.celltype <- Idents(ref.merge) # Note: Clusters not in map become NA/Identity

# Genotype Normalization
ref.merge$ref.genetype <- NA
ref.merge$ref.genetype[grep("WT_Baseline", ref.merge$type)] <- "WT"
ref.merge$ref.genetype[grep("WT_Cuprizone", ref.merge$type)] <- "MS"
ref.merge$ref.genetype[grep("Non-Tg", ref.merge$GSE160512.genotype)] <- "WT"
ref.merge$ref.genetype[grep("PS2APP", ref.merge$GSE160512.genotype)] <- "PS2APP"
ref.merge$ref.genetype[grep("NonTg", ref.merge$GSE153895.genotype)] <- "WT"
ref.merge$ref.genetype[grep("P301L", ref.merge$GSE153895.genotype)] <- "P301L"

# Subset Oligodendrocyte Lineage
target_clusters <- c("1", "4", "5", "6", "14", "15", "21", "25", "33", "3", "24", "26", "32")
ref.merge.oligo <- subset(ref.merge, idents = target_clusters) # Simplified subsetting logic if Idents are reset, otherwise use cell IDs

# Final Plot check
p_final_1 <- DimPlot(ref.merge.oligo, reduction = "umap", group.by = "ref.celltype", label = T)
p_final_2 <- DimPlot(ref.merge.oligo, reduction = "umap", group.by = "ref.genetype", label = T)
plot_grid(p_final_1, p_final_2)

# ------------------------------------------------------------------------------
# 6. Save Data
# ------------------------------------------------------------------------------
saveRDS(ref.merge, file.path(OUT_DIR, "ref.merge.rds"))
saveRDS(ref.merge.oligo, file.path(OUT_DIR, "ref.merge.oligo.rds"))

cat("Analysis Complete.\n")
