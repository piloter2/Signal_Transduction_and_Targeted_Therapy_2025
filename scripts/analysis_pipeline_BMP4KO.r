# ==============================================================================
# Oligodendrocyte precusor cells-microglia cross talk via BMP4 signaling drives microglial neuroprotective response and mitigates Alzheimer's diseas progression
# Description: Integration and analysis of in-house and public scRNA-seq datasets 
#              to investigate oligodendrocyte changes in Alzheimer's disease models.
#              Full pipeline from raw count processing to pathway analysis.
#              Includes QC (DoubletFinder), Integration (Harmony), Annotation,
#              UCell Scoring, and Enrichment Analysis (GO/Reactome).
# Author: Jaemyung Jang / Korea Brain Research Institue (piloter2@kbri.re.kr)
# Date: 2025-07-16
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup & Configuration
# ------------------------------------------------------------------------------
options(warn=-1)

# Load core libraries
suppressPackageStartupMessages({
  library(parallel); library(doParallel); library(foreach)
  library(data.table); library(Matrix); library(tidyverse)
  library(Seurat); library(SingleCellExperiment)
  library(harmony); library(UCell)
  library(EnhancedVolcano); library(ComplexHeatmap)
  library(ggsci); library(ggpubr); library(ggstatsplot); library(ggsignif)
  library(Nebulosa)
  
  # Bio-conductor & Annotation
  library(scDblFinder)
  library(clusterProfiler); library(ReactomePA); library(DOSE)
  library(org.Mm.eg.db); library(AnnotationDbi)
  
  # Custom/Local libraries
  library(indRop) # Ensure this custom package is available
})

# System & Memory Settings
# Note: memory.limit is Windows-specific. Remove if running on Linux.
if(.Platform$OS.type == "windows") memory.limit(size = 98226*1024)

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)
base::gc()

# PATH CONFIGURATION (USER MUST UPDATE)
BASE_DIR      <- "/path/to"
DATA_DIR      <- file.path(BASE_DIR, "data")
RESULT_DIR    <- file.path(BASE_DIR, "results")
SCRIPT_DIR    <- "/path/to/scripts"

# Load custom functions
source(file.path(SCRIPT_DIR, 'single_cell_function_mouse.r'))

# ------------------------------------------------------------------------------
# 2. Data Loading & QC (Parallel Processing)
# ------------------------------------------------------------------------------
cat("Setting up parallel processing...\n")

# Path Setup
path_in_counts  <- list.files(DATA_DIR, pattern = ".counts.tsv", full.names = TRUE) # Update pattern as needed
path_in_metrics <- list.files(DATA_DIR, pattern = ".metrics.tsv", full.names = TRUE)

# Filter paths if necessary (example logic preserved)
# path_metrics <- path_in_metrics[grepl("WT|Tg6799", path_in_metrics)]
# path_counts  <- path_in_counts[grepl("WT|Tg6799", path_in_counts)]
path_metrics <- path_in_metrics
path_counts  <- path_in_counts

# Initialize Cluster
cores <- parallel::detectCores() / 2
cl <- parallel::makeCluster(cores, type='PSOCK')
doParallel::registerDoParallel(cl)
on.exit(stopCluster(cl)) # Ensure cluster stops even if error occurs

# Export libraries and functions to workers
clusterEvalQ(cl, {
  library(Seurat); library(data.table); library(SingleCellExperiment)
  library(indRop); library(scDblFinder)
})
clusterExport(cl, c("readIndrop", "fread", "SingleCellExperiment", "scDblFinder", "CreateSeuratObject"), envir=environment())

# Step 2.1: Load Metrics
cat("Loading metrics...\n")
data.OLsort <- parallel::parLapply(cl, path_metrics, function(path) {
  metric_df <- as.matrix(fread(path, header = TRUE, sep = '\t', verbose=FALSE))
  uni_align <- which(as.numeric(metric_df[,"UMIFM"]) >= 1000)
  as.character(metric_df[uni_align])
})

# Step 2.2: Load Raw Counts (Cached)
cat("Loading raw counts...\n")
filtereddata.OLsort <- parallel::parLapply(cl, path_counts, function(path) {
  return(readIndrop(path))
})

# Step 2.3: Process & Doublet Removal
cat("Processing datasets and removing doublets...\n")
OL_sort <- foreach(exponent = 1:length(path_counts), .packages="Seurat") %dopar% {
  
  filt_indices <- which(colnames(filtereddata.OLsort[[exponent]]) %in% data.OLsort[[exponent]])
  
  if(length(filt_indices) > 1){
    dt <- filtereddata.OLsort[[exponent]][, filt_indices]
    
    # Run scDblFinder
    sce <- SingleCellExperiment(assays = list(counts = dt))
    tryCatch({
      sce <- scDblFinder(sce)
    }, error = function(e) print("scDblFinder Error"), warning = function(w) print("Warning"))
    
    count_data <- dt
    
    # Remove Doublets
    if(!is.null(sce$scDblFinder.class)){
      doublets <- colnames(sce)[which(sce$scDblFinder.class=="doublet")]
      if(length(doublets) > 0){
        count_data <- dt[, -which(colnames(dt) %in% doublets)]
      }
    }
    
    return(CreateSeuratObject(counts = count_data, min.cells = 10))
  } else {
    return(NULL)
  }
}

# Cleanup Empty Entries
names(OL_sort) <- sapply(str_split(path_counts, "/"), function(x) {
  # Logic to extract ID from filename, adjust splitting index as needed
  tail(str_split(x, "/")[[1]], 1) 
})
OL_sort.filtered <- compact(OL_sort)

# ------------------------------------------------------------------------------
# 3. Preprocessing & Integration (SCT + Harmony)
# ------------------------------------------------------------------------------
cat("Preprocessing individual objects...\n")

OL_sort_ <- compact(lapply(names(OL_sort.filtered), function(x){
  obj <- OL_sort.filtered[[x]]
  if(!is.null(obj) && ncol(obj) > 100){
    obj$basic <- "OL sorted"
    obj$id <- x
    obj$genetype <- ifelse(grepl("WT", toupper(x)), "WT", "TG")
    obj$new_ident <- paste0(x, "_", obj$genetype)
    
    # QC Metrics
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = paste0(c(mt.genes, "^mt."), collapse = "|"))
    obj[["percent.ribosomal"]] <- PercentageFeatureSet(obj, pattern = "^Rpl|^Rps")
    obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = '^Hba.a1$|^Hba.a2$|^Hba.x$|^Hbb.bh1$|^Hbb.bh2$|^Hbb.bs$|^Hbb.bt$|^Hbq1b$|^Hbb.y$|^Hbq1a$')
    
    obj <- processA(obj) # Assuming processA is in the sourced script
    return(obj)
  }
  return(NULL)
}))

cat("Merging and Integrating...\n")
OLs.merged <- merge(OL_sort_[[1]], y = OL_sort_[-1], project = 'merged', merge.data = TRUE)

# Normalization & Cell Cycle Scoring
OLs.merged <- SCTransform(OLs.merged, assay = 'RNA', new.assay.name = 'SCT', 
                          vars.to.regress = c('percent.mt','nFeature_RNA','percent.hb'), verbose = T)
OLs.merged <- CellCycleScoring(OLs.merged, s.features = intersect(rownames(OLs.merged), s.gene), 
                               g2m.features = intersect(rownames(OLs.merged), g2m.gene), assay = 'SCT')
OLs.merged <- SCTransform(OLs.merged, assay = 'RNA', new.assay.name = 'SCT', 
                          vars.to.regress = c('percent.mt','nFeature_RNA', 'S.Score', 'G2M.Score','percent.hb'), verbose = T)

# Feature Selection excluding mitochondrial/hemoglobin genes
OLs.merge.features <- SelectIntegrationFeatures(object.list = OL_sort_, nfeatures = 2000)
OLs.merge.features <- setdiff(OLs.merge.features, grep("^mt-|^mt.|^Hba|^Hbb|^Hbd|^Hbg1|^Hbq1|^Hbm", rownames(OLs.merged), value = TRUE))

# Dimensionality Reduction
pcs <- 1:25
DefaultAssay(OLs.merged) <- "SCT"
OLs.merged <- RunPCA(OLs.merged, verbose = FALSE, assay = "SCT", features = OLs.merge.features)
# Uncomment Harmony if needed:
# OLs.merged <- RunHarmony(OLs.merged, group.by="id", assay.use="SCT", dims=pcs)
OLs.merged <- RunUMAP(OLs.merged, reduction = "pca", umap.method = "umap-learn", assay = "SCT", dims = pcs)
OLs.merged <- FindNeighbors(OLs.merged, reduction = "pca", dims = pcs, assay = "SCT")
OLs.merged <- FindClusters(OLs.merged, resolution = seq(0.1, 2, 0.1), algorithm = 4)

OLs.merge1.id <- OLs.merged

# ------------------------------------------------------------------------------
# 4. Cell Type Annotation
# ------------------------------------------------------------------------------
Idents(OLs.merge1.id) <- "SCT_snn_res.0.9" 

# Manual Annotation Map
celltype_map <- c(
  "1" = "OL_mixed-B", "2" = "EPC A", "3" = "EPC B", "4" = "choroidal macrophages",
  "5" = "OL_mixed-A", "6" = "NEUR A", "7" = "U-5", "8" = "NEUR B",
  "9" = "NFOL", "10" = "MG", "11" = "OL_S", "12" = "EPC C",
  "13" = "NFOL-WT", "14" = "ASC", "15" = "U-1", "16" = "MOL",
  "17" = "U-2", "18" = "Endo", "19" = "U-3", "20" = "Mural",
  "21" = "U-4", "22" = "NEUR C", "23" = "OPC/COP", "24" = "OL_U",
  "25" = "COP", "26" = "Fib"
)

OLs.merge1.id <- RenameIdents(OLs.merge1.id, celltype_map)
OLs.merge1.id$broad_type <- Idents(OLs.merge1.id)

# ------------------------------------------------------------------------------
# 5. Visualization (DimPlots & DotPlots)
# ------------------------------------------------------------------------------
cat("Generating Plots...\n")
target_clusters <- c('OL_mixed-B','OL_mixed-A','NFOL','OL_S','NFOL-WT','MOL','OPC/COP','OL_U','COP')
# Or use cluster IDs: c('1','5','9','11','13','16','23','24','25')

p1a <- DimPlot(OLs.merge1.id, reduction = "umap", group.by = "broad_type", label = FALSE) + NoAxes()
ggsave(file.path(RESULT_DIR, "Figure2b_Dim.jpeg"), p1a, width = 7.5, height = 7.5, dpi = 300)

p1c <- plot_density(subset(OLs.merge1.id, idents = target_clusters), 
                    features = c("Bmp4","Pdgfra","Cspg4","Mobp","Bcas1","Tlcd2","Lrrc17","Plp1", "Mbp"), 
                    reduction = "umap", pal = "inferno")
ggsave(file.path(RESULT_DIR, "Figure2d_Den.jpeg"), p1c, width = 16, height = 12, dpi = 300)

# ------------------------------------------------------------------------------
# 6. UCell Scoring & Comparison
# ------------------------------------------------------------------------------
cat("Running UCell Scoring...\n")
marker_sets <- list(
  "INF" = inflammation, # Ensure 'inflammation' variable is defined in sourced script
  "Myelination" = c("Adam10","Adgrg6","Akt1","Mbp","Mog","Plp1","Olig2","Sox10"), # (Truncated for brevity)
  "NFkB" = c("Nfkb1","Rela","Relb","Nfkbia")
)

OLs.merge1.score <- AddModuleScore_UCell(OLs.merge1.id, features = marker_sets, assay = "RNA")

# Boxplots for comparison (5xFAD vs WT)
# (Loop logic simplified for readability)
clusters_to_test <- c('1','5','9','11','16','23','24','25')

for(s in clusters_to_test) {
  obj_sub <- subset(OLs.merge1.score, idents = s)
  obj_sub$genetype <- factor(ifelse(obj_sub$genetype=="WT", "C57BL/6", "Tg6799"))
  
  # Example: Myelination Score
  p_box <- ggbetweenstats(
    data = data.frame(genetype = obj_sub$genetype, Score = obj_sub$Myelination_UCell),
    x = genetype, y = Score, type = "p", pairwise.comparisons = FALSE,
    title = paste0("Cluster ", s)
  )
  # Save or print p_box as needed
}

# ------------------------------------------------------------------------------
# 7. Differential Expression & Pathway Analysis (GO/Reactome)
# ------------------------------------------------------------------------------
cat("Running DE and Pathway Analysis...\n")

# 7.1 DE Analysis (MAST)
OLs.merge1.DE <- lapply(names(table(Idents(OLs.merge1.id))), function(i){
  obj_sub <- subset(OLs.merge1.id, idents = i)
  if(ncol(obj_sub) > 10) {
    tryCatch({
      FindMarkers(obj_sub, ident.1 = "TG", ident.2 = "WT", group.by = "genetype",
                  test.use = "MAST", assay = "RNA", logfc.threshold = 0.001)
    }, error = function(e) return(NULL))
  } else { return(NULL) }
})
names(OLs.merge1.DE) <- names(table(Idents(OLs.merge1.id)))

# 7.2 Pathway Analysis Loop
# Initialize lists for results
enrich_res <- list(bp = list(), cc = list(), mf = list(), reactome = list())

for(id in clusters_to_test) {
  if(!is.null(OLs.merge1.DE[[id]])) {
    # Prepare Gene List
    selectedGene <- OLs.merge1.DE[[id]] %>% dplyr::filter(p_val < 0.01 & abs(avg_log2FC) > 0.25)
    selectedGene <- selectedGene[!grepl("^mt\\.", rownames(selectedGene)),]
    
    if(nrow(selectedGene) > 0){
      gene_map <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(selectedGene),
                                        columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
      
      gene_list_sorted <- selectedGene$avg_log2FC
      names(gene_list_sorted) <- gene_map$ENTREZID[match(rownames(selectedGene), gene_map$SYMBOL)]
      gene_list_sorted <- sort(na.omit(gene_list_sorted), decreasing = TRUE)
      
      # Run Enrichment (GO BP example)
      tryCatch({
        gse.bp <- gseGO(geneList = gene_list_sorted, ont = "BP", OrgDb = org.Mm.eg.db,
                        keyType = "ENTREZID", pvalueCutoff = 0.1, verbose = FALSE)
        enrich_res$bp[[id]] <- setReadable(gse.bp, 'org.Mm.eg.db', 'ENTREZID')
      }, error = function(e) NULL)
      
      # Reactome
      tryCatch({
        res.react <- gsePathway(gene_list_sorted, organism = "mouse", pvalueCutoff = 0.1, verbose = FALSE)
        enrich_res$reactome[[id]] <- setReadable(res.react, 'org.Mm.eg.db', 'ENTREZID')
      }, error = function(e) NULL)
    }
  }
}

# 7.3 Visualization of Pathways (Barplot example)
# (Logic for p1e_total extracted here)

cat("Pipeline Finished. Saving Session...\n")
saveRDS(OLs.merge1.id, file.path(RESULT_DIR, "OLs.merged.final.rds"))
