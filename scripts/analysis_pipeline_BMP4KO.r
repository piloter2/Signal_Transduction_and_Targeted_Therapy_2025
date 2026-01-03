# ==============================================================================
# Oligodendrocyte precusor cells-microglia cross talk via BMP4 signaling drives microglial neuroprotective response and mitigates Alzheimer's diseas progression
# Description: Integration and analysis of in-house and public scRNA-seq datasets 
#              to investigate oligodendrocyte changes in Alzheimer's disease models.
#              Full pipeline from raw count processing to pathway analysis.
#              Includes QC (DoubletFinder), Integration (Harmony), Annotation,
#              UCell Scoring, and Enrichment Analysis (GO/Reactome).
# Author: Jaemyung Jang / Korea Brain Research Institue (piloter2@kbri.re.kr)
# Date: 2023-07-06
# ==============================================================================
###### UMIFM = 500, decontamn = 0.001, BMP4KO_final_2023-07-06 15-58-10.Rdata

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
cor_colagen <- foreach(exponent = 1:length(path_counts), .packages="Seurat") %dopar% {
  
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

# # Cleanup Empty Entries
# names(cor_colagen) <- sapply(str_split(path_counts, "/"), function(x) {
#   # Logic to extract ID from filename, adjust splitting index as needed
#   tail(str_split(x, "/")[[1]], 1) 
# })
# cor_colagen <- compact(cor_colagen)

genetype <- c(rep("Pdgfa-cre;BMP4fl.fl;5xFAD",5),
              rep("5xFAD (C57BL/6)",7))
filename <- unlist(lapply(fileListCNTs,function(x){strsplit(x,".counts.tsv.gz")[[1]][1]}))
date <- c()
date[grep("210706", filename)] <- "210706"
date[grep("210910", filename)] <- "210910"
date[c(1:2,6:7)] <- "200820"

for(i in 1:length(cor_colagen)) {
  cor_colagen[[i]]$genetype <- genetype[i]
  cor_colagen[[i]]$date <- date[i]
  cor_colagen[[i]]$platform <- "inDrops"
  cor_colagen[[i]]$id <- filename[i]
  # cor_colagen[[i]]$basic <- "cor_colagen"
  cor_colagen[[i]]$new_ident <- paste0(filename[i],"_", genetype[i],"_",date[i])  
    
}

# # soupX
if (!require("SoupX", quietly = TRUE))
    install.packages("SoupX")

require(SoupX)

scNoDrop <-list()
for(i in 1:length(cor_colagen)) {
  toc <- GetAssayData(cor_colagen[[i]]) 
  # tod <-  cor_colagen_[[i]]

  scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)

  soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
  scNoDrop[[i]] = setSoupProfile(scNoDrops, soupProf)
}

dreopGenes <- list()
for(i in 1:length(cor_colagen)) {
  sc <- scNoDrop[[i]] 
  #sc = setContaminationFraction(sc, 0.1)
  #head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
  
  dreopGenes[[i]]<-sc$soupProfile[which(sc$soupProfile$est>0.001), ] %>% rownames()
}


 dropGenes.filtered <- lapply(dreopGenes, function(gene) {
   return(setdiff(gene, grep("^mt.|^Hba.|^Hbb.", gene, value = TRUE)))
 })
 

cor_colagen.filtered<-list()
for(i in 1:length(cor_colagen)) {
  cor_colagen.filtered[[i]] <- cor_colagen[[i]][-which(rownames(cor_colagen[[i]]) %in%  dropGenes.filtered[[i]]),]
}

cor_colagen_ <- lapply(cor_colagen.filtered, processA)
             
 ## Reference GSE71585
GSE71585.RefSeq.counts <- read.table(file = "/home/choelab/working/public_dataset/GSE71585/GSE71585_RefSeq_counts.csv", sep =",", header =T)
GSE71585.meta <- read.table(file = "/home/choelab/working/public_dataset/GSE71585/GSE71585_Clustering_Results.csv", sep =",", header =T)
GSE71585.meta <- GSE71585.meta %>% tibble::column_to_rownames('sample_title')
rownames(GSE71585.meta) <- gsub("\\-", ".", rownames(GSE71585.meta))
GSE71585.RefSeq.counts <- GSE71585.RefSeq.counts %>% tibble::column_to_rownames('gene')

rownames(GSE71585.RefSeq.counts)[grep("Epb4", rownames(GSE71585.RefSeq.counts))] <- c('Epb41', 'Epb41l1', 'Epb41l2', 'Epb41l3', 'Epb41l4a', 'Epb41l4b', 'Epb41l5', 'Epb42', 'Epb49')
rownames(GSE71585.RefSeq.counts) <- gsub("\\-", ".", rownames(GSE71585.RefSeq.counts))

GSE71585.RefSeq.count <- CreateSeuratObject(counts =  GSE71585.RefSeq.counts, project = "GSE71585.RefSeq", min.cells = 0,  meta.data = GSE71585.meta)
# GSE71585.RefSeq.counts_ <- CreateSeuratObject(counts =  GSE71585.RefSeq.counts %>% tibble::column_to_rownames('gene'), project = "GSE71585.RefSeq", min.cells = 10)

base::gc()
GSE71585.RefSeq.counts_ <- processB(GSE71585.RefSeq.count)
GSE71585.RefSeq.counts_$basic <- GSE71585.RefSeq.counts_$genetype <- GSE71585.RefSeq.counts_$id <- "GSE71858"
GSE71585.RefSeq.counts_$platform <- "SmartSeq2"
             
# ------------------------------------------------------------------------------
# 3. Preprocessing & Integration (SCT + Harmony)
# ------------------------------------------------------------------------------

##### PROCESS 1 - integration & SCT & Clustering

cor.71585.merge <- merge(GSE71585.RefSeq.counts_, y= cor_colagen_ , project = 'bmp4ko', merge.data = FALSE)
cor.71585.merge.features <- SelectIntegrationFeatures(object.list = c(cor_colagen_ , GSE71585.RefSeq.counts_), nfeatures = 5000)
cor.71585.merge.features <- setdiff(cor.71585.merge.features, grep(c("^mt-|^mt.|^Hba|^Hbb|^Hbd|^Hbg1|^Hbq1|^Hbm"), cor.71585.merge.features, value =T))

base::gc()

cor.71585.merge <- cor.71585.merge  %>%
                  SCTransform(assay = 'RNA',      
                              new.assay.name = 'SCT', 
                              vars.to.regress = c('percent.mt', 'nFeature_RNA','percent.Hb'),
                              verbose = T)  %>% # normalize data with SCTransform()
                  CellCycleScoring(s.features = intersect(rownames(cor.71585.merge),s.gene),  # Perform cell cycle analysis
                                   g2m.features = intersect(rownames(cor.71585.merge),g2m.gene),
                                   assay = 'SCT',
                                   set.ident = TRUE) %>%
                  SCTransform(assay = 'RNA',
                              new.assay.name = 'SCT',
                              vars.to.regress = c('percent.mt', 'nFeature_RNA', 'S.Score', 'G2M.Score','percent.Hb'),
                              verbose = T) 


pcs <- 50
set.seed(1234)
if (!require("harmony", quietly = TRUE))
    install.packages("harmony")
require(harmony)
cor.71585.merge.genetype <-cor.71585.merge %>%
                          RunPCA(verbose = FALSE, assay = "SCT", features = cor.71585.merge.features) %>%
                          RunHarmony(group.by=c("platform", "id"), assay.use = "SCT", dims=1:pcs, max.iter.harmony = 100,max.iter.cluster = 200) %>% 
                          RunUMAP(reduction = "harmony", assay = "SCT", umap.method = "umap-learn", dims=1:pcs)  %>% 
                          FindNeighbors(reduction = "harmony", assay = "SCT", dims=1:pcs, graph.name = "RNA_snn") %>%
                          FindClusters(resolution = seq(0.1,2,0.1), algorithm = 4, graph.name = "RNA_snn")

Idents(cor.71585.merge.genetype) <- "RNA_snn_res.1.2" # original

# measuring silhouette score
cor.71585.merge.genetype.sub <- subset(cor.71585.merge.genetype, subset= platform == "inDrops")
pcs <- 1:50
for(res in seq(0.1,2,0.1)) {
  
  require(cluster)
  distance_matrix <- dist(Embeddings(cor.71585.merge.genetype.sub[['pca']])[, pcs])
  # distance_matrix <- dist(Embeddings(hipp_merged[['pca']]))
  clusters <- eval(parse(text=sprintf("cor.71585.merge.genetype.sub@meta.data$RNA_snn_res.%s",res)))
  
  #dat.combined@meta.data$seurat_clusters
  silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
  cor.71585.merge.genetype.sub@meta.data$silhouette_score <- silhouette[,3]
  
  mean_silhouette_score <- mean(cor.71585.merge.genetype.sub@meta.data$silhouette_score)
  median_silhouette_score <- median(cor.71585.merge.genetype.sub@meta.data$silhouette_score)
  mad_silhouette_score <- mad(cor.71585.merge.genetype.sub@meta.data$silhouette_score)
  
  df <- data.frame('resolution'=res, 'mean'=mean_silhouette_score,'median'=median_silhouette_score,'mad'=mad_silhouette_score)
  
  print(df)
  
}         

saveRDS(cor_colagen, "/path/to/data/final_Rdata/bmp4ko_rawdata_UMIFM500.RDS")
saveRDS(cor.71585.merge.genetype.sub, file.path(RESULT_DIR, "bmp4ko_workingOBJ_20230703.RDS"))
saveRDS(cor.71585.merge.genetype, file.path(RESULT_DIR, "bmp4ko_merged_20230704.RDS"))
