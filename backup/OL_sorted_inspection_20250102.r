
options(warn=-1)

require(reticulate)
use_python("/usr/bin/python3")
# use_miniconda("r-reticulate")
# miniconda_update()
# py_install("scikit-learn")
# py_install("pandas")
# py_install("numexpr")
# py_install("louvain")
# py_install("leidenalg")
# py_install("umap-learn==0.5.3")
# # py_config()
py_discover_config()

cat(names(reticulate::import("umap")), sep="\n")



# Pass TRUE if you want to see progress output on some of Monocle 3's operations
DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)
options(future.globals.maxSize= 1000e6)

base::gc()
# base::rm(list = ls()) # Clear the environment

setwd('~/R')
require(readxl)
source('~/R/single_cell_function_20230601.R')

################# path - OL sorted

# pathIMZ <- file.path("/mnt","workingA","OL_sorted","reference")
# fileListCNTs <- list.files(pathIMZ, pattern = ".counts.tsv")
# fileListMETRICS <- list.files(pathIMZ, pattern = ".metrics.tsv")

# path_in_counts.raw <- paste0(pathIMZ,"/",fileListCNTs)
# path_in_metrics.raw <- paste0(pathIMZ,"/",fileListMETRICS)


# path_metrics <- c(path_in_metrics.raw[grepl("WT",path_in_metrics.raw)],
#                   path_in_metrics.raw[grepl("Tg6799",path_in_metrics.raw)])

# path_counts<-c(path_in_counts.raw[grepl("WT",path_in_counts.raw)],
#                path_in_counts.raw[grepl("Tg6799",path_in_counts.raw)])

################# read data
# library(scDblFinder)
# library(indRop)
# require(doParallel) 
# require(foreach)
# library(data.table)

# system <- Sys.info()['sysname']
# cores <- parallel::makeCluster(detectCores()/2, type='PSOCK')
# cl <- NULL

# cl <- parallel::makeCluster(getOption('cl.cores', cores))
# doParallel::registerDoParallel(cl)
#   ###################################################################
#   clusterEvalQ(cl, c( library(Seurat), library(data.table),library(SingleCellExperiment) ,library(indRop),library(scDblFinder)))
#   clusterExport(cl, c("readIndrop","fread","SingleCellExperiment","scDblFinder","CreateSeuratObject"),
#                 envir=environment())
  
#   data.OLsort <- parallel::parLapply(cl,path_metrics,  function(path) {
#     metric_OLsort      <- as.matrix(fread(path, header = TRUE, sep = '\t',verbose=FALSE))
#     uni_align_OLsort <- which(as.numeric(metric_OLsort[,"UMIFM"]) >= 1000)
#     as.character(metric_OLsort[uni_align_OLsort])
#   })
  
  
#   cacheParallel.OLsort <- function() {
#     parallel::parLapply(cl,path_counts, function(path) {
#       return(readIndrop(path))
#     })
#   }
  
#   elapse <- base::system.time(filtereddata.OLsort<- cacheParallel.OLsort())
  
  
#   OL_sort<-foreach::foreach(exponent = 1:length(path_counts), .packages="Seurat")  %dopar% {
#     dt<-c()
#     runData<-c()
#     filt <- which(colnames(filtereddata.OLsort[[exponent]]) %in% data.OLsort[[exponent]])
#     if(length(filt) > 1){
#         dt<-filtereddata.OLsort[[exponent]][,filt]
#         print(dim(dt))
#         runData<-dt
#         print(dim(runData))
#         runData.sce <- SingleCellExperiment(assays = list(counts = runData))
#         tryCatch(runData.sce <- scDblFinder(runData.sce),
#                  error = function(e) print("Error"),
#                  warning = function(w) print("Warning"),
#                  finally = NULL)
        

#         count <- runData
        
#         if(!is.null(runData.sce$scDblFinder.class)){
#             doubletlists <- colnames(runData.sce)[which(runData.sce$scDblFinder.class=="doublet")]
#             if(length(doubletlists) > 0){
#                 count <- runData[,-which(colnames(dt) %in% doubletlists)]
#             }
#         }
        
        
#         # count
#         CreateSeuratObject(counts = count, min.cells = 10)


#     } else{
        
#     }
#   }
  
  
  

################################################################
# on.exit(stopCluster(cl))
# stopCluster(cl)
# detach("package:indRop", unload=TRUE)
# detach("package:scDblFinder", unload=TRUE)

# names(OL_sort) <- sapply(str_split(path_counts, "/"), function(x){str_split(x[7], "\\.")[[1]][1]})

# saveRDS(OL_sort, "/home/choelab/working/OL_sorted/OL_sorted_rawdata_UMIFM500.RDS")
# OL_sort.filtered <- compact(OL_sort)


################# OL sorted
# names(OL_sort) <- sapply(str_split(path_counts, "/"), function(x){str_split(x[7], "\\.")[[1]][1]})
# OL_sort_ <- compact(sapply(names(OL_sort.filtered), function(x){
#     if(!is.null(OL_sort.filtered[[x]]) && dim(OL_sort.filtered[[x]])[2] > 100){
#         OL_sort.filtered[[x]]$basic <- "OL sorted"
#         OL_sort.filtered[[x]]$id <- x
#         OL_sort.filtered[[x]]$new_ident <- paste0(x, "_", ifelse(grepl("WT",toupper(x)),"WT","TG"))
#         OL_sort.filtered[[x]]$genetype <- ifelse(grepl("WT",toupper(x)),"WT","TG")
#         OL_sort.filtered[[x]][["percent.mt"]] <- PercentageFeatureSet(OL_sort.filtered[[x]], pattern = paste0(c(mt.genes,"^mt."), collapse = "|"))
#         OL_sort.filtered[[x]][["percent.ribosomal"]] <- PercentageFeatureSet(OL_sort.filtered[[x]], pattern = "^Rpl|^Rps")
#         OL_sort.filtered[[x]][["percent.hb"]] <- PercentageFeatureSet(OL_sort.filtered[[x]], pattern = paste0('^Hba.a1$|^Hba.a2$|^Hba.x$|^Hbb.bh1$|^Hbb.bh2$|^Hbb.bs$|^Hbb.bt$|^Hbq1b$|^Hbb.y$|^Hbq1a$'))
#         # paste0(unique(unlist(lapply(OL_sort, function(x){grep("^Hb" ,rownames(x), value = T)}))), collapse = "|" )

#         OL_sort.filtered[[x]]$log10GenesPerUMI <- log10(OL_sort.filtered[[x]]$nFeature_RNA) / log10(OL_sort.filtered[[x]]$nCount_RNA)
#         OL_sort.filtered[[x]] <- processA(OL_sort.filtered[[x]])
#         # if(!is.null(OL_sort[[x]]) && dim(OL_sort[[x]])[2] > 100){
#              OL_sort.filtered[[x]]
#         # } else{}
#         # OL_sort[[x]]
#     }  
    
# }))

# PROCESS 1
# OLs.merge1 <- merge(OL_sort_[[1]], y= c(OL_sort_[2:length(OL_sort_)]),project = 'merged',merge.data=T)

# OLs.merged <- OLs.merge1 %>%
#                   SCTransform(assay = 'RNA',      
#                               new.assay.name = 'SCT',  
#                               vars.to.regress = c('percent.mt','nFeature_RNA','percent.hb'),
#                               verbose = T)  %>% # normalize data with SCTransform()
#                   CellCycleScoring(s.features = intersect(rownames(OLs.merge1),s.gene),  # Perform cell cycle analysis
#                                    g2m.features = intersect(rownames(OLs.merge1),g2m.gene),
#                                    assay = 'SCT',
#                                    set.ident = TRUE) %>%
#                   SCTransform(assay = 'RNA',
#                               new.assay.name = 'SCT',
#                               vars.to.regress = c('percent.mt','nFeature_RNA', 'S.Score', 'G2M.Score','percent.hb'),
#                               verbose = T)      

# OLs.merge.features <- SelectIntegrationFeatures(object.list = OL_sort_, nfeatures = 2000)
# OLs.merge.features <- setdiff(OLs.merge.features, grep(c("^mt-|^mt.|^Hba|^Hbb|^Hbd|^Hbg1|^Hbq1|^Hbm"), rownames(OLs.merged), value =T))



# pcs <- 1:25
# set.seed(1234)
# library(harmony)
# DefaultAssay(OLs.merged) <- "SCT"
# OLs.merge1.id <- OLs.merged %>%
#                           RunPCA(verbose = FALSE, assay = "SCT", features = OLs.merge.features  ) %>%
#                         #   RunHarmony(group.by=c("id"), assay.use = "SCT", dims=pcs) %>% 
#                           RunUMAP(reduction = "pca", umap.method = "umap-learn", assay = "SCT",  dims=pcs) %>%
#                           FindNeighbors(reduction = "pca", dims=pcs, assay = "SCT") %>%
#                           FindClusters(resolution = seq(0.1,2,0.1), algorithm = 4)


# # saveRDS(OLs.merge1.id, "/home/choelab/working/OL_sorted/final_Rdata/OL_sorted_filtered_20230720.RDS")

# DimPlot(OLs.merge1.id, reduction = "umap", group.by = "id", label =T, label.size = 5, pt.size = 0.01,raster=FALSE) + 
#            theme(plot.title = element_text(size = 12, face = "bold",hjust = 0),  axis.line = element_blank())
# DimPlot(OLs.merge1.id, group.by = "genetype", label =T, label.size = 5, pt.size = 0.01,raster=FALSE) + 
#            theme(plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank())


# for(res in seq(0.1,2,0.1)) {
  
#   require(cluster)
#   distance_matrix <- dist(Embeddings(OLs.merge1.id[['pca']])[, pcs])
#   # distance_matrix <- dist(Embeddings(hipp_merged[['pca']]))
#   clusters <- eval(parse(text=sprintf("OLs.merge1.id@meta.data$SCT_snn_res.%s",res)))
  
#   #dat.combined@meta.data$seurat_clusters
#   silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
#   OLs.merge1.id@meta.data$silhouette_score <- silhouette[,3]
  
#   mean_silhouette_score <- mean(OLs.merge1.id@meta.data$silhouette_score)
#   median_silhouette_score <- median(OLs.merge1.id@meta.data$silhouette_score)
#   mad_silhouette_score <- mad(OLs.merge1.id@meta.data$silhouette_score)
  
#   df <- data.frame('resolution'=res, 'mean'=mean_silhouette_score,'median'=median_silhouette_score,'mad'=mad_silhouette_score)
  
#   print(df)
  
# }     

# Idents(OLs.merge1.id) <- "SCT_snn_res.0.9" 
# OLs.merge1.id.marker <- FindAllMarkers(OLs.merge1.id,  assay = "RNA", logfc.threshold = 0.25)

# marker <- OLs.merge1.id.marker %>% dplyr::filter(p_val_adj < 0.1 & avg_log2FC > 0.25)

# intersect(marker %>% dplyr::filter(cluster=='23') %>% pull(gene),marker %>% dplyr::filter(cluster=='25') %>% pull(gene))

# OLs.merge1.id.marker[grep("Sox10",rownames(OLs.merge1.id.marker)),] %>% dplyr::filter(p_val_adj < 0.1 & avg_log2FC > 0.25)


# OLs.merge1.id.submarker <- FindAllMarkers(subset(OLs.merge1.id, idents=c('1','5','9','11','13','16','23','24','25')),  
#                                             test.use = "MAST",
#                                             assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE)


# OLs.merge1.id.submarker %>% dplyr::filter(p_val_adj < 0.1 & avg_log2FC > 0.5 & cluster=='9')

# marker.sub <- OLs.merge1.id.submarker %>% dplyr::filter(p_val_adj < 0.1 & avg_log2FC > 1)#%>% top_n(n=10, wt=avg_log2FC)

# Idents(OLs.merge1.id) <- "SCT_snn_res.0.9" 
# celltype <- c()
# celltype[1] <- "OL_mixed-B" # 1 
# celltype[2] <- "EPC A" # 2 
# celltype[3] <- "EPC B"# 3
# celltype[4] <- "choroidal macrophages"#4 
# celltype[5] <- "OL_mixed-A"# 5 
# celltype[6] <- "NEUR A"# 6
# celltype[7] <- "U-5"# 7
# celltype[8] <- "NEUR B"# 7
# celltype[9] <- "NFOL"# 8
# celltype[10] <- "MG"# 9
# celltype[11] <- "OL_S"# 10
# celltype[12] <- "EPC C"
# celltype[13] <- "NFOL-WT"
# celltype[14] <- "ASC"
# celltype[15] <- "U-1"
# celltype[16] <- "MOL"
# celltype[17] <- "U-2"
# celltype[18] <- "Endo"
# celltype[19] <- "U-3"
# celltype[20] <- "Mural"
# celltype[21] <- "U-4"
# celltype[22] <- "NEUR C"
# celltype[23] <- "OPC/COP"
# celltype[24] <- "OL_U"
# celltype[25] <- "COP"
# celltype[26] <- "Fib"

# names(celltype) <- c(1:26)

# for(i in names(celltype)){
#     OLs.merge1.id@meta.data$broad_type[colnames(OLs.merge1.id) %in% WhichCells(subset(OLs.merge1.id, ident = i))] <- celltype[i]
# }

# table(OLs.merge1.id$broad_type) %>% length()



#GSE160512

path_GSE160512 <- file.path("/mnt","workingA","OL_sorted","reference","GSE160512","raw")
fileListCNTs <- list.files(path_GSE160512, pattern = ".txt.gz")
# fileListCNTs <- list.files(pathIMZ, pattern = "GSE180041_exonic_norm_umi.csv.gz")
path_in_counts.raw <- paste0(path_GSE160512,"/",fileListCNTs)
# raw <- read.table(path_in_counts.raw, header =T, sep=",")
names(fileListCNTs) <- path_in_counts.raw

GSE160512.meta <- read.table("/mnt/workingA/OL_sorted/reference/GSE160512/GSE160512_PS2APPTsneCoordinatesAndClusterAssignments.txt.gz", header =T, sep="\t")
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
#     if(names(table(GSE160512.obj[[y]]$GSE160512.genotype)) == "Non-Tg"){
#         GSE160512.obj.OL.wt[[y]] <- GSE160512.obj[[y]][, grepl("^oligo|^OPC|^COP", GSE160512.obj[[y]]$GSE160512.celltype)]
#         GSE160512.obj.OL.wt2[[y]] <- GSE160512.obj[[y]][, grepl("^oligo|^COP", GSE160512.obj[[y]]$GSE160512.celltype)]

#     }
    

    # GSE160512.obj.OL[[y]] <- GSE160512.obj[[y]][, grep("^oligo|^OPC|^COP", GSE160512.obj[[y]]$GSE160512.celltype)]
    # GSE160512.obj.OL2[[y]] <- GSE160512.obj[[y]][, grep("^oligo|^COP", GSE160512.obj[[y]]$GSE160512.celltype)]
    # GSE160512.obj.OL3[[y]] <- merge(GSE160512.obj.OL2[[y]], GSE160512.obj.OL[[y]][,GSE160512.obj.OL[[y]]$GSE160512.celltype == "OPC"][,1:200])
}

GSE160512.oligo <- list()

for(y in seq_len(length(GSE160512.obj))){
    GSE160512.oligo[[y]] <- GSE160512.obj[[y]][, grep("^oligo|^OPC|^COP", GSE160512.obj[[y]]@meta.data$celltype)]
}


#GSE153895

path_GSE153895 <- file.path("/mnt","workingA","OL_sorted","reference","GSE153895","raw")
fileListCNTs <- list.files(path_GSE153895, pattern = ".txt.gz")
# fileListCNTs <- list.files(pathIMZ, pattern = "GSE180041_exonic_norm_umi.csv.gz")
path_in_counts.raw <- paste0(path_GSE153895,"/",fileListCNTs)
# raw <- read.table(path_in_counts.raw, header =T, sep=",")
names(fileListCNTs) <- path_in_counts.raw

GSE153895.meta <- read.table("/mnt/workingA/OL_sorted/reference/GSE153895/GSE153895_TsneCoordinatesAndClusterAssignments.txt.gz", header =T, sep="\t")
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



#GSE148676

# GSM4476524	WT_Baseline_rep1
# GSM4476526	WT_Baseline_rep2
# GSM4476528	WT_Cuprizone_4w_rep1
# GSM4476532	WT_Cuprizone_4w_rep2
# GSM4476536	WT_Baseline_rep3
# GSM4476538	WT_Cuprizone_4w_rep3
# GSM4476540	WT_Cuprizone_4w_rep4


type <- as.data.frame(c("WT_Baseline_rep1", "WT_Baseline_rep2", "WT_Cuprizone_4w_rep1", "WT_Cuprizone_4w_rep2", "WT_Baseline_rep3", "WT_Cuprizone_4w_rep3", "WT_Cuprizone_4w_rep4")) 
              
rownames(type) <- c("GSM4476524", "GSM4476526","GSM4476528","GSM4476532","GSM4476536", "GSM4476538", "GSM4476540")
colnames(type) <- "type"


path_GSE148676 <- file.path("/mnt","workingA","OL_sorted","reference","GSE148676")
fileListCNTs <- list.files(path_GSE148676, pattern = ".txt.gz")
# fileListCNTs <- list.files(pathIMZ, pattern = "GSE180041_exonic_norm_umi.csv.gz")
path_in_counts.raw <- paste0(path_GSE148676,"/",fileListCNTs)
# raw <- read.table(path_in_counts.raw, header =T, sep=",")
names(fileListCNTs) <- path_in_counts.raw


require(biomaRt)
annot <- c()
GSE148676.obj <- list()

ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
for(x in path_in_counts.raw){
    require(indRop)
    temp <- readIndrop(x)
    annot <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'),filters = 'ensembl_gene_id',values = colnames(temp),mart = ensembl)
    annot <- annot[-c(which((annot$mgi_symbol == '') | duplicated(annot$mgi_symbol) | duplicated(annot$ensembl_gene_id))),]
    temp <- temp[,annot$ensembl_gene_id]
    colnames(temp) <- gsub("-","\\.", annot$mgi_symbol)
    GSE148676.obj[[fileListCNTs[x]]] <- CreateSeuratObject(counts = t(temp), min.cells = 0)               
}


# GSE148676.obj.OL <- list()


for(x in path_in_counts.raw){
    y <- fileListCNTs[x]
    GSE148676.obj[[y]]@meta.data$id <- y
    GSE148676.obj[[y]]@meta.data$basic <- "GSE148676"
    GSE148676.obj[[y]]@meta.data$type <- type[strsplit(y,"_")[[1]][1],]
    GSE148676.obj[[y]]@meta.data$percent.mt <- PercentageFeatureSet(GSE148676.obj[[y]], pattern = paste0(c(mt.genes,"^mt."), collapse = "|"))
    GSE148676.obj[[y]]@meta.data$percent.hb <- PercentageFeatureSet(GSE148676.obj[[y]], pattern = paste0(c('Hbp1','Hbq1a','Hbb.y','Hbb.bh1','Hbb.bs','Hba.x','Hba.a2','Hba.a1','Hbq1b','Hbb.bt','Hbb.bh2','Hbb.bh3','Hba.ps4','Hbb.bh0'), collapse = "|"))

    # GSE148676.obj[[y]]@meta.data$celltype <- GSE160512.meta[colnames(GSE148676.obj[[y]]),]$allCells.cluster.interpretation
    # GSE148676.obj[[y]]@meta.data$GSE148676.genotype <- GSE160512.meta[colnames(GSE148676.obj[[y]]),]$genotype
#     if(names(table(GSE160512.obj[[y]]$GSE160512.genotype)) == "Non-Tg"){
#         GSE160512.obj.OL.wt[[y]] <- GSE160512.obj[[y]][, grepl("^oligo|^OPC|^COP", GSE160512.obj[[y]]$GSE160512.celltype)]
#         GSE160512.obj.OL.wt2[[y]] <- GSE160512.obj[[y]][, grepl("^oligo|^COP", GSE160512.obj[[y]]$GSE160512.celltype)]

#     }
    

    # GSE160512.obj.OL[[y]] <- GSE160512.obj[[y]][, grep("^oligo|^OPC|^COP", GSE160512.obj[[y]]$GSE160512.celltype)]
    # GSE160512.obj.OL2[[y]] <- GSE160512.obj[[y]][, grep("^oligo|^COP", GSE160512.obj[[y]]$GSE160512.celltype)]
    # GSE160512.obj.OL3[[y]] <- merge(GSE160512.obj.OL2[[y]], GSE160512.obj.OL[[y]][,GSE160512.obj.OL[[y]]$GSE160512.celltype == "OPC"][,1:200])
}


table(GSE148676.obj[[1]]@meta.data$celltype)
table(GSE148676.obj[[3]]$type)

# PROCESS 1
GSE148676.merge <- merge(GSE148676.obj[[1]], y= c(GSE148676.obj[2:length(GSE148676.obj)]), project = 'merged',merge.data=FALSE)

GSE148676.merge <- GSE148676.merge %>%
                  SCTransform(assay = 'RNA',      
                              new.assay.name = 'SCT',  
                              vars.to.regress = c("percent.mt","percent.hb","nFeature_RNA"),
                              verbose = T)   # normalize data with SCTransform()

GSE148676.features <- SelectIntegrationFeatures(object.list = GSE148676.obj, nfeatures = 5000)

GSE148676.feature <- GSE148676.features[-grep("^mt.|Rik$|^Hba.|^Hbb.",GSE148676.features)]

table(GSE148676.merge$id)
pcs <- 1:50
set.seed(1234)
require(harmony)

GSE148676.merged <- GSE148676.merge %>%
                          RunPCA(verbose = FALSE, assay = "SCT", features = GSE148676.feature ) %>%
                          RunHarmony(group.by=c("id"), assay.use = "SCT", dims=pcs) %>% 
                          RunUMAP(reduction = "harmony", umap.method = "umap-learn", assay = "SCT",  dims=pcs) %>%
                          FindNeighbors(reduction = "harmony", dims=pcs, assay = "SCT") %>%
                          FindClusters(resolution = seq(0.1,1,0.1), algorithm = 2)

Idents(GSE148676.merged) <- "SCT_snn_res.0.4"

DimPlot(GSE148676.merged, 
        reduction = "umap", #group.by = "basic", 
        label = T, repel = T, 
        label.size = 6, pt.size = 0.01) + 
  theme( axis.title = element_blank(), 
        plot.title = element_text(size = 14, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
  ggtitle(paste0("GSE148676 - cells : ",comma(length(colnames(GSE148676.merged)),format="d")," cells"))

GSE118918.marker <- c(
"Pdgfra","Vcan","Olig1","Ptprz1","Olig2", # opc
"Plp1","Cldn11","Cnp","Tspan2","Mobp", 
"Trf","Apod","Ermn","Car2",
"Ptgds","Plp1","Mal","Cryab","Car2", # 
"Fyn","Sirt2","Tuba1a","Bcas1","Enpp6", # cop
"Hexb","Ctss","C1qb","C1qc","Cx3cr1", # micro
"Lyz2","Mrc1","Pf4","Ctsb","F13a1", # macrophage
"Ly6c1","Ly6a","Flt1","Itm2a","Slco1a4", # Endothelial
"Rgs5","Acta2","Myl9","Tagln","Cald1", # Mural
"Dcn","Col1a2","Igf2","Igfbp2","Pcolce",#VLMC
"Ccdc153","Rarres2","Tmem212","Nnat","Rsph1", # Ependyma    
"Cpe","Apoe","Cst3","Mt1","Dbi","Aldoc","Gja1","Slc1a2","Slc1a3","Clu", 
"Ccnd2","Tubb5","Dlx1","Tmsb10","Tuba1a", #NSC
"Hmgb2","Top2a","Hmgn2","Tubb5","H2afz", # NSC    
"Sst","Npy","Resp18","Nap1l5","Nos1", # NEURON
"Lypd1","Meg3","Atp1b1","Snhg11","Olfm1",# NEURON
"Penk","Arpp21","Ppp3ca","Atp2b1","Hpca" #NEURON
)

oligo.mark <- c("Sirt2","Bmp4","Tnr","Fyn", "Enpp6", "Neu4",  # COP 
                  "C1ql1", "Olig1", "Pdgfra", "Cspg4",# OPC
                  "Apod", "Trf", "Plp1", "Mbp", "Fth1","Mag", "Mog","Cldn11","Mobp", # MOL
                  "Tmeff2","Idh1","Neu4", "Cdh20", "Tmem163", "Nfasc", "Mpzl1", # NFOL
                  "Mal","Opalin", "Ptgds","Bfsp2" #MFOL
                 )

DotPlot(GSE148676.merged, 
        # group.by = "SCT_snn_res.0.3", 
        idents = names(table(Idents(GSE148676.merged)))[table(Idents(GSE148676.merged)) > 100],
        features = unique(GSE118918.marker), assay="SCT" , cols ="RdBu") + coord_flip() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


GSE148676.oligo <- subset(GSE148676.merged, idents = c(1, 3, 6, 7, 17, 19, 30))

GSE148676.oligo.subset <- list()

for(i in seq_len(length(GSE148676.obj))) {
    GSE148676.oligo.subset[[i]] <- subset(GSE148676.obj[[i]], cells = intersect(colnames(GSE148676.obj[[i]]), colnames(GSE148676.oligo)))
}

# GSE148676.merged <- GSE148676.merge %>%
#   RunPCA(verbose = TRUE, features = GSE148676.feature) %>%  
#   RunHarmony(reduction.use = "pca", group.by.vars="Donor", max.iter.harmony = 100, max.iter.cluster=250) %>% #"platform",
#   RunUMAP(reduction = "harmony", umap.method = "umap-learn", dims=1:50) %>%
#   FindNeighbors(reduction = "harmony", dims=1:50) 
#GSE266687

data <- lapply(str_sub(list.files("/mnt/workingA/OL_sorted/reference/GSE266687/", pattern = "_matrix.mtx.gz"), 1, str_locate(list.files("/mnt/workingA/OL_sorted/reference/GSE266687/", pattern = "_matrix.mtx.gz"),"_matrix.mtx.gz")[,1]), function(x) {
    Read10X_GEO("/mnt/workingA/OL_sorted/reference/GSE266687/",x)
})

GSE266687 <- unlist(data)

# require(biomaRt)
# annot <- c()
GSE266687.obj <- list()

ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
for(x in names(GSE266687)){
    temp <- GSE266687[[x]]
    # annot <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'),filters = 'ensembl_gene_id',values = rownames(temp),mart = ensembl)
    # annot <- annot[-c(which((annot$mgi_symbol == '') | duplicated(annot$mgi_symbol) | duplicated(annot$ensembl_gene_id))),]
    # temp <- temp[annot$ensembl_gene_id,]
    rownames(temp) <- gsub("-","\\.", rownames(temp))
    GSE266687.obj[[x]] <- CreateSeuratObject(counts = temp, min.cells = 0)               
}


for(i in seq_len(length(GSE266687.obj))){
    GSE266687.obj[[i]]@meta.data$id <- str_sub(list.files("/mnt/workingA/OL_sorted/reference/GSE266687/", pattern = "_matrix.mtx.gz"), 1, 10)[i]
    GSE266687.obj[[i]]@meta.data$basic <- "GSE266687"
    GSE266687.obj[[i]]@meta.data$type <- str_sub(list.files("/mnt/workingA/OL_sorted/reference/GSE266687/", pattern = "_matrix.mtx.gz"), 12, str_locate(list.files("/mnt/workingA/OL_sorted/reference/GSE266687/", pattern = "_matrix.mtx.gz"),"_matrix.mtx.gz")[,1]-1)[i]
    GSE266687.obj[[i]]@meta.data$percent.mt <- PercentageFeatureSet(GSE266687.obj[[i]], pattern = paste0(c(mt.genes,"^mt."), collapse = "|"))
    GSE266687.obj[[i]]@meta.data$percent.hb <- PercentageFeatureSet(GSE266687.obj[[i]], pattern = paste0(c('Hbp1','Hbq1a','Hbb.y','Hbb.bh1','Hbb.bs','Hba.x','Hba.a2','Hba.a1','Hbq1b','Hbb.bt','Hbb.bh2','Hbb.bh3','Hba.ps4','Hbb.bh0'), collapse = "|"))
}

table( GSE266687.obj[[1]]@meta.data$id)
table( GSE266687.obj[[1]]@meta.data$type)

GSE266687.merge <- merge(GSE266687.obj[[1]], y= c(GSE266687.obj[2:length(GSE266687.obj)]), project = 'merged',merge.data=FALSE)

GSE266687.merge <- GSE266687.merge %>%
                  SCTransform(assay = 'RNA',      
                              new.assay.name = 'SCT',  
                              vars.to.regress = c("percent.mt","percent.hb","nFeature_RNA"),
                              verbose = T)   # normalize data with SCTransform()

GSE266687.features <- SelectIntegrationFeatures(object.list = GSE266687.obj, nfeatures = 5000)

GSE266687.feature <- GSE266687.features[-grep("^mt.|Rik$|^Hba.|^Hbb.",GSE266687.features)]

table(GSE266687.merge$id)
pcs <- 1:50
set.seed(1234)
require(harmony)

GSE266687.merged <- GSE266687.merge %>%
                          RunPCA(verbose = FALSE, assay = "SCT", features = GSE266687.feature ) %>%
                          RunHarmony(group.by=c("id"), assay.use = "SCT", dims=pcs) %>% 
                          RunUMAP(reduction = "harmony", umap.method = "umap-learn", assay = "SCT",  dims=pcs) %>%
                          FindNeighbors(reduction = "harmony", dims=pcs, assay = "SCT") %>%
                          FindClusters(resolution = seq(0.1,1,0.1), algorithm = 2)


Idents(GSE266687.merged) <- "SCT_snn_res.0.8"

DimPlot(GSE266687.merged, 
        reduction = "umap", #group.by = "basic", 
        label = T, repel = T, 
        label.size = 6, pt.size = 0.01) + 
  theme( axis.title = element_blank(), 
        plot.title = element_text(size = 14, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
  ggtitle(paste0("GSE266687 - cells : ",comma(length(colnames(GSE266687.merged)),format="d")," cells"))



DotPlot(GSE266687.merged, 
        # group.by = "SCT_snn_res.0.3", 
        # idents = names(table(Idents(GSE266687.merged)))[table(Idents(GSE266687.merged)) > 100],
        features = unique(GSE118918.marker), assay="SCT" , cols ="RdBu") + coord_flip() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


GSE266687.oligo <- subset(GSE266687.merged, idents = c(3, 15, 16, 18))

intersect(colnames(GSE266687.obj[[4]]), str_sub(colnames(GSE266687.oligo), 1, 18)[which(str_sub(colnames(GSE266687.oligo), 20, 25) == names(table(str_sub(colnames(GSE266687.oligo), 20, 25)))[5])])

# table(GSE266687.oligo$type)

GSE266687.oligo.subset <- list()

GSE266687.oligo.subset[[1]] <- subset(GSE266687.obj[[1]], cells = intersect(colnames(GSE266687.obj[[1]]), str_sub(colnames(GSE266687.oligo), 1, 18)[which(str_sub(colnames(GSE266687.oligo), 20, 25) == names(table(str_sub(colnames(GSE266687.oligo), 20, 25)))[1])]))
GSE266687.oligo.subset[[2]] <- subset(GSE266687.obj[[10]], cells = intersect(colnames(GSE266687.obj[[10]]), str_sub(colnames(GSE266687.oligo), 1, 18)[which(str_sub(colnames(GSE266687.oligo), 20, 25) == names(table(str_sub(colnames(GSE266687.oligo), 20, 25)))[2])]))
for(i in c(2:9)) {
    GSE266687.oligo.subset[[i+1]] <- subset(GSE266687.obj[[i]], cells = intersect(colnames(GSE266687.obj[[i]]), str_sub(colnames(GSE266687.oligo), 1, 18)[which(str_sub(colnames(GSE266687.oligo), 20, 25) == names(table(str_sub(colnames(GSE266687.oligo), 20, 25)))[i+1])]))
}

################## GSE278199

GSE278199 <- readRDS("/mnt/workingA/OL_sorted/reference/GSE278199/GSE278199_OPC_cells.rds")

    # annot <- getBM(attributes = c('mgi_symbol','ensembl_gene_id'),filters = 'ensembl_gene_id',values = rownames(temp),mart = ensembl)
    # annot <- annot[-c(which((annot$mgi_symbol == '') | duplicated(annot$mgi_symbol) | duplicated(annot$ensembl_gene_id))),]
    # temp <- temp[annot$ensembl_gene_id,]
rownames(GSE278199) <- gsub("-","\\.", rownames(GSE278199))
GSE278199.obj <- CreateSeuratObject(counts = GSE278199, min.cells = 0)               

GSE278199.obj@meta.data$id[which(str_sub(colnames(GSE278199.obj), 1, 1) == names(table(str_sub(colnames(GSE278199.obj), 1, 1)))[1])] <- "naïve"
GSE278199.obj@meta.data$id[which(str_sub(colnames(GSE278199.obj), 1, 1) == names(table(str_sub(colnames(GSE278199.obj), 1, 1)))[2])] <- "saline_CPZ3"
GSE278199.obj@meta.data$id[which(str_sub(colnames(GSE278199.obj), 1, 1) == names(table(str_sub(colnames(GSE278199.obj), 1, 1)))[3])] <- "XPro1595_CPZ3"
GSE278199.obj@meta.data$basic <- "GSE278199"
GSE278199.obj@meta.data$percent.mt <- PercentageFeatureSet(GSE278199.obj, pattern = paste0(c(mt.genes,"^mt."), collapse = "|"))
GSE278199.obj@meta.data$percent.hb <- PercentageFeatureSet(GSE278199.obj, pattern = paste0(c('Hbp1','Hbq1a','Hbb.y','Hbb.bh1','Hbb.bs','Hba.x','Hba.a2','Hba.a1','Hbq1b','Hbb.bt','Hbb.bh2','Hbb.bh3','Hba.ps4','Hbb.bh0'), collapse = "|"))

#######################################################################


OL.merge <- merge(GSE278199.obj, y= c(GSE266687.oligo, GSE153895.oligo, GSE160512.oligo), project = 'merged',merge.data=FALSE) #GSE148676.oligo.subset, 

OL.merge <- OL.merge %>%
                  SCTransform(assay = 'RNA',      
                              new.assay.name = 'SCT',  
                              vars.to.regress = c("percent.mt","percent.hb","nFeature_RNA"),
                              verbose = T) 


OL.features <- SelectIntegrationFeatures(object.list = c(GSE278199.obj, GSE266687.oligo, GSE153895.oligo, GSE160512.oligo), nfeatures = 3000) #GSE148676.oligo.subset, 
OL.feature <- OL.features[-grep("^mt.|Rik$|^Hba.|^Hbb.",OL.features)]


sum(table(OL.merge$basic))
pcs <- 1:25
set.seed(1234)
require(harmony)

OL.merged <- OL.merge %>%
                          RunPCA(verbose = FALSE, assay = "SCT", features = OL.feature ) %>%
                          RunHarmony(group.by=c("basic"), assay.use = "SCT", dims=pcs, max.iter.harmony = 100, max.iter.cluster=250) %>% 
                          RunUMAP(reduction = "harmony", umap.method = "umap-learn", assay = "SCT",  dims=pcs) %>%
                          FindNeighbors(reduction = "harmony", dims=pcs, assay = "SCT") %>%
                          FindClusters(resolution = seq(0.1,1,0.1), algorithm = 2)

Idents(OL.merged) <- "SCT_snn_res.0.8"

 names(table(OL.merged$id))[8]
# c(7:9, 13:24)
table(subset(OL.merged, subset = id == names(table(OL.merged$id))[15])$celltype)


table(OL.merged$type)

require(ggsci)
p1s<-DimPlot(subset(OL.merged, subset = id == names(table(OL.merged$id))[7] | 
                           id == names(table(OL.merged$id))[8] | 
                           id == names(table(OL.merged$id))[9] | 
                           id == names(table(OL.merged$id))[13] | 
                           id == names(table(OL.merged$id))[14] | 
                           id == names(table(OL.merged$id))[15] | 
                           id == names(table(OL.merged$id))[16] | 
                           id == names(table(OL.merged$id))[17] | 
                           id == names(table(OL.merged$id))[18] | 
                           id == names(table(OL.merged$id))[19] | 
                           id == names(table(OL.merged$id))[20] | 
                           id == names(table(OL.merged$id))[21] | 
                           id == names(table(OL.merged$id))[22] | 
                           id == names(table(OL.merged$id))[23] | 
                           id == names(table(OL.merged$id))[24] )
, 
        cols = pal_npg("nrc")(10),
        group.by = "new_cell_type", 
        reduction = "umap", #group.by = "basic", 
        label = T, repel = T, 
        label.size = 6, pt.size = 0.01) + 
  theme( axis.title = element_blank(), 
        plot.title = element_text(size = 14, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
  ggtitle(paste0("OL - cells : ",comma(length(colnames(subset(OL.merged, subset = id == names(table(OL.merged$id))[7] | 
                           id == names(table(OL.merged$id))[8] | 
                           id == names(table(OL.merged$id))[9] | 
                           id == names(table(OL.merged$id))[13] | 
                           id == names(table(OL.merged$id))[14] | 
                           id == names(table(OL.merged$id))[15] | 
                           id == names(table(OL.merged$id))[16] | 
                           id == names(table(OL.merged$id))[17] | 
                           id == names(table(OL.merged$id))[18] | 
                           id == names(table(OL.merged$id))[19] | 
                           id == names(table(OL.merged$id))[20] | 
                           id == names(table(OL.merged$id))[21] | 
                           id == names(table(OL.merged$id))[22] | 
                           id == names(table(OL.merged$id))[23] | 
                           id == names(table(OL.merged$id))[24] ))),format="d")," cells")) + 
  NoLegend()

table(OL.merged$basic)

p1s_1<-DimPlot(subset(OL.merged, subset = id == names(table(OL.merged$id))[7] | 
                           id == names(table(OL.merged$id))[8] | 
                           id == names(table(OL.merged$id))[9] | 
                           id == names(table(OL.merged$id))[13] | 
                           id == names(table(OL.merged$id))[14] | 
                           id == names(table(OL.merged$id))[15] | 
                           id == names(table(OL.merged$id))[16] | 
                           id == names(table(OL.merged$id))[17] | 
                           id == names(table(OL.merged$id))[18] | 
                           id == names(table(OL.merged$id))[19] | 
                           id == names(table(OL.merged$id))[20] | 
                           id == names(table(OL.merged$id))[21] | 
                           id == names(table(OL.merged$id))[22] | 
                           id == names(table(OL.merged$id))[23] | 
                           id == names(table(OL.merged$id))[24] )
, 
        # cols = pal_npg("nrc")(10),
        group.by = "basic", 
        reduction = "umap", #group.by = "basic", 
        label = T, repel = T, 
        label.size = 6, pt.size = 0.01) + 
  theme( axis.title = element_blank(), 
        plot.title = element_text(size = 14, face = "bold",hjust = 0), axis.line = element_blank(), 
        axis.text = element_blank(),  axis.ticks = element_blank()) + 
  ggtitle(paste0("OL - cells : ",comma(length(colnames(subset(OL.merged, subset = id == names(table(OL.merged$id))[7] | 
                           id == names(table(OL.merged$id))[8] | 
                           id == names(table(OL.merged$id))[9] | 
                           id == names(table(OL.merged$id))[13] | 
                           id == names(table(OL.merged$id))[14] | 
                           id == names(table(OL.merged$id))[15] | 
                           id == names(table(OL.merged$id))[16] | 
                           id == names(table(OL.merged$id))[17] | 
                           id == names(table(OL.merged$id))[18] | 
                           id == names(table(OL.merged$id))[19] | 
                           id == names(table(OL.merged$id))[20] | 
                           id == names(table(OL.merged$id))[21] | 
                           id == names(table(OL.merged$id))[22] | 
                           id == names(table(OL.merged$id))[23] | 
                           id == names(table(OL.merged$id))[24] ))),format="d")," cells")) + 
  NoLegend()

PS2APP <- subset(OL.merged, subset = id == names(table(OL.merged$id))[7] | 
                           id == names(table(OL.merged$id))[8] | 
                           id == names(table(OL.merged$id))[9] | 
                           id == names(table(OL.merged$id))[13] | 
                           id == names(table(OL.merged$id))[14] | 
                           id == names(table(OL.merged$id))[15] | 
                           id == names(table(OL.merged$id))[16] | 
                           id == names(table(OL.merged$id))[17] | 
                           id == names(table(OL.merged$id))[18] | 
                           id == names(table(OL.merged$id))[19] | 
                           id == names(table(OL.merged$id))[20] | 
                           id == names(table(OL.merged$id))[21] | 
                           id == names(table(OL.merged$id))[22] | 
                           id == names(table(OL.merged$id))[23] | 
                           id == names(table(OL.merged$id))[24] )

PS2APP <- PrepSCTFindMarkers(PS2APP)

COP.markers <- FindMarkers(PS2APP,
            group.by = "celltype", 
            test.use = "MAST",
            ident.1 = "COP")

write.xlsx(COP.markers, "/mnt/workingA/OL_sorted/results/fig1S_b.xlsx", rowNames = TRUE)

selGenes <- COP.markers %>% dplyr::filter(p_val_adj < 1e-300 & avg_log2FC > 2) %>% rownames()

setdiff(selGenes, grep("^Gm|Rik$", selGenes, value = TRUE))

require(EnhancedVolcano)
ggsave("/mnt/workingA/OL_sorted/results/fig1S_a.pdf",p1s, width = 7.5, height = 7.5, units= "in" )
pdf("/mnt/workingA/OL_sorted/results/fig1S_b.pdf", width = 5, height = 7.5, family = "ArialMT")
    EnhancedVolcano(COP.markers,
                          lab = rownames(COP.markers),
                          x = 'avg_log2FC',
                          y = 'p_val_adj'   ,
                        #   xlim = c(-max(abs(COP.markers$avg_log2FC))-0.1, max(abs(COP.markers$avg_log2FC))+0.1),
                        #   ylim = c(-1, max(-log10(COP.markers$p_val_adj))+0.1),
                          title = "COP vs OPC/Oligodendrocytes",
                          subtitle = paste0(" COP : ",comma(length(colnames(subset(PS2APP, subset = celltype == "COP"))),format = "d")," cells -  Others : ", comma(length(colnames(subset(PS2APP, subset = celltype != "COP"))),format = "d"), " cells"),
                          xlab = bquote(~Log[2]~ 'fold change'),
                          ylab = bquote(~Log[10]~ 'Adjust P-value'),
                          pCutoff = 1e-300,
                          FCcutoff = 2,
                          pointSize = 0.25,
                          labSize = 3.2,
                        #   selectLab = setdiff(selGenes, grep("^Gm|Rik$", selGenes, value = TRUE)),
                        selectLab = rownames(dfs.combined[-16,c(2:4)]), 
                          # shapeCustom = keyvals.shape,
                          #colCustom = keyvals.colour,
                          # col=c("grey30", "forestgreen", "royalblue", "red2"),
                          #                colCustom = keyvals,
                          colAlpha = 0.5,
                          legendPosition = 'bottom',
                          drawConnectors = TRUE,
                          widthConnectors = 0.05,
                          colConnectors = 'black',
                          border = 'full',
                          borderWidth = 0.25,
                          # labCol = "black",
                          # legendPosition = 'none',
                          axisLabSize = 10,
                          titleLabSize = 12,
                          subtitleLabSize = 10,
                          captionLabSize = 10,
                          maxoverlapsConnectors = 50)
dev.off()

# DotPlot(OL.merged, 
#         # group.by = "SCT_snn_res.0.3", 
#         # idents = names(table(Idents(OL.merged)))[table(Idents(OL.merged)) > 100],
#         features = unique(oligo.mark), 
#         assay="SCT" , cols ="RdBu") + 
#         coord_flip() + 
#            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

OL.merged <- PrepSCTFindMarkers(OL.merged)
OL.COP <- subset(OL.merged, idents = c(11,29))
# OL.COP <- PrepSCTFindMarkers(OL.COP)

OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "0"))] <- "OPC"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "13"))] <- "OPC"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "14"))] <- "OPC"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "17"))] <- "OPC"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "27"))] <- "OPC"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "1"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "2"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "3"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "4"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "5"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "6"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "7"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "8"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "9"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "10"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "12"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "15"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "16"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "18"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "19"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "20"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "21"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "22"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "23"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "24"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "25"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "26"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "28"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "30"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "31"))] <- "Oligodendrocytes"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "11"))] <- "COP"
OL.merged@meta.data$new_cell_type[colnames(OL.merged) %in% colnames(subset(OL.merged, idents = "29"))] <- "COP"



# OL.COP@meta.data$stimulation <- rep("xx", length(colnames(OL.COP))
# OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[1])] <- "Mertk-WT_Baseline"
# OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[2])] <- "Mertk-WT_Baseline"
# OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[3])] <- "Mertk-WT_Cuprizone_4w"
# OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[4])] <- "Mertk-WT_Cuprizone_4w"
# OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[5])] <- "Mertk-WT_Baseline"
# OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[6])] <- "Mertk-WT_Cuprizone_4w"
# OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[7])] <- "Mertk-WT_Cuprizone_4w"
# OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[8])] <- "Mertk-KO_Cuprizone_4+3w"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[1])] <- "PS2/APP/P301L"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[2])] <- "PS2/APP/P301L/TREM2KO"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[3])] <- "PS2/APP/P301L"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[4])] <- "PS2/APP/P301L/TREM2KO"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[5])] <- "PS2/APP/P301L"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[6])] <- "P301L"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[7])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[8])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[9])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[10])] <- "PS2/APP/P301L/TREM2KO"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[11])] <- "P301L"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[12])] <- "P301L"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[13])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[14])] <- "PS2/APP"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[15])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[16])] <- "PS2/APP"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[17])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[18])] <- "PS2/APP"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[19])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[20])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[21])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[22])] <- "PS2/APP"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[23])] <- "WT"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[24])] <- "PS2/APP"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[25])] <- "WT/Ctrl_0wks"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[26])] <- "WT/Ctrl_3wks"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[27])] <- "WT/CupRap_3wks"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[28])] <- "WT/CupRap_3wks"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[29])] <- "NesCre/Ctrl_3wks"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[30])] <- "NesCre/Ctrl_3wks"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[31])] <- "NesCre/CupRap_3wks"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[32])] <- "NesCre/CupRap_3wks"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[33])] <- "WT/CupRap_0wk"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[34])] <- "WT/CupRap_0wk"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[35])] <- "WT/CTRL"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[36])] <- "WT/CPZ3"
OL.COP@meta.data$genetype[which(OL.COP@meta.data$id == names(table(OL.COP@meta.data$id))[37])] <- "WT/XPro1595_CPZ3"

# save.image("/mnt/workingA/OL_sorted/final_Rdata/OL_sorted_new_20240103.RData")
# load("/mnt/workingA/OL_sorted/final_Rdata/OL_sorted_new_20240103.RData")
names(table(OL.COP@meta.data$genetype))

table(OL.COP@meta.data$basic)



OL.COP.DEP.DM_with_CupRap_0wks <- FindMarkers(OL.COP, 
                                        ident.1 = names(table(OL.COP@meta.data$genetype))[12],
                                        ident.2 = names(table(OL.COP@meta.data$genetype))[10],
                                        group.by = "genetype", 
                                        test.use = "MAST",
                                        recorrect_umi = FALSE)

OL.COP.DEP.DM_with_CupRap_3wks <- FindMarkers(OL.COP, 
                                        ident.1 = names(table(OL.COP@meta.data$genetype))[13],
                                        ident.2 = names(table(OL.COP@meta.data$genetype))[11],
                                        group.by = "genetype", 
                                        test.use = "MAST",
                                        recorrect_umi = FALSE)

OL.COP.DEP.DM_with_CPZ3 <- FindMarkers(OL.COP, 
                                        ident.1 = names(table(OL.COP@meta.data$genetype))[8],
                                        ident.2 = names(table(OL.COP@meta.data$genetype))[9],
                                        group.by = "genetype", 
                                        test.use = "MAST",
                                        recorrect_umi = FALSE)

# OL.COP.DEP.DM_on_Mertk_WT <- FindMarkers(OL.COP, 
#                                         ident.1 = names(table(OL.COP@meta.data$genetype))[3],
#                                         ident.2 = names(table(OL.COP@meta.data$genetype))[2],
#                                         group.by = "genetype", 
#                                         test.use = "MAST",
#                                         recorrect_umi = FALSE)

OL.COP.DEP.PS2_APP <- FindMarkers(OL.COP, 
                                        ident.1 = names(table(OL.COP@meta.data$genetype))[4],
                                        ident.2 = names(table(OL.COP@meta.data$genetype))[7],
                                        group.by = "genetype", 
                                        test.use = "MAST",
                                        recorrect_umi = FALSE)

OL.COP.DEP.tau <- FindMarkers(OL.COP, 
                                        ident.1 = names(table(OL.COP@meta.data$genetype))[5],
                                        ident.2 = names(table(OL.COP@meta.data$genetype))[7],
                                        group.by = "genetype", 
                                        test.use = "MAST",
                                        recorrect_umi = FALSE)

table(subset(OL.COP, subset = genetype == names(table(OL.COP@meta.data$genetype))[4])$basic)
table(subset(OL.COP, subset = genetype == names(table(OL.COP@meta.data$genetype))[5])$basic)

df <- merge(OL.COP.DEP.DM_with_CPZ3[grep(paste0(setdiff(selGenes, grep("^Gm|Rik$", selGenes, value = TRUE)), collapse = "$|^"),rownames(OL.COP.DEP.DM_with_CPZ3)),]  %>% dplyr::select("avg_log2FC"), 
          OL.COP.DEP.PS2_APP[grep(paste0(setdiff(selGenes, grep("^Gm|Rik$", selGenes, value = TRUE)), collapse = "$|^"),rownames(OL.COP.DEP.PS2_APP)),] %>% dplyr::select("avg_log2FC"), by = 0)

OL.COP.DEP.DM_with_CPZ3[grep("Bmp4",rownames(OL.COP.DEP.DM_with_CPZ3)),]
OL.COP.DEP.PS2_APP[grep("Bmp4",rownames(OL.COP.DEP.PS2_APP)),]
OL.COP.DEP.tau[grep("Bmp4",rownames(OL.COP.DEP.tau)),]

dfs <- df[,c(2:3)]
rownames(dfs) <- df$Row.names

dfs.combined <- merge(dfs, OL.COP.DEP.tau[grep(paste0(setdiff(selGenes, grep("^Gm|Rik$", selGenes, value = TRUE)), collapse = "$|^"),rownames(OL.COP.DEP.tau)),] %>% dplyr::select("avg_log2FC"), by = 0)
colnames(dfs.combined) <- c("Genes","GSE278199 vs GSE278199", "GSE160512 vs GSE153895 and GSE160512", "GSE153895 vs  GSE153895 and GSE160512")

rownames(dfs.combined) <-  dfs.combined$Genes
dim(dfs.combined)

write.xlsx(OL.COP.DEP.DM_with_CPZ3, "/mnt/workingA/OL_sorted/results/fig1S_column_1_rev.xlsx", rowNames = TRUE)
write.xlsx(OL.COP.DEP.PS2_APP, "/mnt/workingA/OL_sorted/results/fig1S_column_2.xlsx", rowNames = TRUE)
write.xlsx(OL.COP.DEP.tau, "/mnt/workingA/OL_sorted/results/fig1S_column_3.xlsx", rowNames = TRUE)

# paletteLength <- 11
# myColor <- colorRampPalette(c("royalblue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
# get.breaks <- c(seq(min(dfs.combined[,c(2:4)]), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(dfs.combined[,c(2:4)])/paletteLength, max(dfs.combined[,c(2:4)]), length.out=floor(paletteLength/2)))


# nbreaks <- round(max(abs(dfs.combined[,c(2:4)])))*2+1 # as many breaks as you want
# get.breaks <- seq(from = range(dfs.combined[,c(2:4)])[1],
#                   to = range(dfs.combined[,c(2:4)])[2],
#                   length.out = nbreaks)

# get.breaks <- seq(from = -round(max(abs(dfs.combined[,c(2:4)])))/2,
#                   to = round(max(abs(dfs.combined[,c(2:4)])))/2, 
#                                      length.out = nbreaks)
# require(RColorBrewer)
# require(pheatmap)
# pheatmap(dfs.combined[,c(2:4)],  
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     color = rev(brewer.pal(n = 13, name = "RdYlBu"))[c(1:5, 7, 9, 11)],
#     legend_breaks =  get.breaks)

require(ComplexHeatmap)

val = round(max(abs(dfs.combined[,c(2:4)])))/2

pdf("/mnt/workingA/OL_sorted/results/fig1S_c.pdf", width = 5, height = 15, family = "ArialMT")
ComplexHeatmap::pheatmap(
  dfs.combined[-16,c(2:4)], #25,24
  name ="Average Fold Change",
  cellwidth = 25,
  cellheight = 25,
#   fontsize = 10,
#   angle_col = 270,
  cluster_cols =F,
   breaks = c(-val,0, val),
  cluster_rows = F,
  annotation_names_col =T,
  border_color="white",
#    legend_breaks =  seq(-val,val,2*val/5), 
#    legend_labels = seq(-val,val,2*val/5),                                    
  main = "Average expression")
dev.off()

DimPlot(OL.COP, 
           reduction = "umap", 

        #     cols = c(brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Set3")),
            label =T, label.size = 4, pt.size = 0.01,raster=FALSE)+
            theme(legend.position = "none",  axis.title = element_blank(), 
                  plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
                  axis.text = element_blank(),  axis.ticks = element_blank()) + 
                  ggtitle("") + NoAxes()

DimPlot(OL.merged, 
            reduction = "umap", 
            
            # cells.highlight = WhichCells(subset(OLs.merge1.id, idents=c('1','5','9','11','13','16','23','24','25'))),
            # cols.highlight = "#DE2D26",
            # sizes.highlight = 0.5,
             label =F, label.size = 4, pt.size = 0.01,raster=FALSE)+
            theme(legend.position = "none",  axis.title = element_blank(), 
                  plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
                  axis.text = element_blank(),  axis.ticks = element_blank()) + 
                ggtitle(paste0(comma(length(colnames(OL.merged)),format = "d")," cells"))+
                  NoAxes()

##################################################################################        
p1b<-DimPlot(subset(OLs.merge1.id, idents=c('1','5','9','11','13','16','23','24','25')), 
           reduction = "umap", 
           group.by = "broad_type", 
           cols = c(brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Set3"))[c(1,5,9,11,13,16,23,24,25)],
            label =T, label.size = 4, pt.size = 1,raster=FALSE)+
            theme(legend.position = "none",  axis.title = element_blank(), 
                  plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
                  axis.text = element_blank(),  axis.ticks = element_blank()) +
            ggtitle(paste0(comma(length(colnames(subset(OLs.merge1.id, idents=c('1','5','9','11','13','16','23','24','25')))),format = "d")," cells"))+
                  NoAxes()


DotPlot(OLs.merge1.id, 
        group.by = "broad_type",
        idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
        # idents = c('1','9','11','13','16','23','24','25'), #,'24','25'        
        features = unique(c(GSE118918.marker,"Epyc")), assay="SCT" , cols ="RdBu") + coord_flip() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pS1b<-DotPlot(OLs.merge1.id, 
        group.by = "broad_type",
        idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
        # idents = c('1','9','11','13','16','23','24','25'), #,'24','25'        
        features = c("Slc1a2","Sparcl1","Corin","Olfr807","Cldn5","Flt1","Epyc","Cst3","Plp1","Acta2","Ccdc186","Ttc9b","B3galt2","Prdm11","Creg1","Safb",
        "Tlcd2","Pdxk","Mobp","Fth1","Bcas1","Mbp","Bc1","Lrrc17","Cldn11",
        "Rpl10","Rps4x","Gng10","Eef1a1","Ppan","Atxn3","Prpf8"), 
        assay="SCT", cols = c("RdBu"), dot.scale = 8) + 
        ylab("") + xlab("") + 
                      RotatedAxis() + coord_flip()+
                      theme(legend.title = element_text(size = 12, hjust = 0),
                            legend.text = element_text(size = 12),
                            axis.title = element_text(size = 15, face = "bold"),
                            axis.text = element_text(size = 18))



DotPlot(OLs.merge1.id, 
        group.by = "broad_type",
        # idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
        idents = c('1','5','9','11','13','16','23','24','25'), #,'24','25'        
        features = c("Mobp","Plp1","Tlcd2","Pdxk","Fth1","Bcas1","Mbp","Bc1","Lrrc17","Cldn11","Fos","Cspg4","Pdgfra","Gpr17"), assay="SCT" , cols ="RdBu") + coord_flip() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


DotPlot(subset(OLs.merge1.id, idents=c('1','5','9','11','13','16','23','24','25')), 
        group.by = "broad_type",
        idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
        # idents = c('1','9','11','13','16','23','24','25'), #,'24','25'        
        # features = unique(setdiff(marker.sub$gene,grep("^Gm|Rik$|^mt.",marker.sub$gene, value = TRUE)))
        features = c("Plp1","Mobp","Mbp","Bcas1","Cstb","Cldn11","Cnp","Mal","Car2"), assay="SCT" , cols ="RdBu") + coord_flip() + 

           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

require(Nebulosa)
p1c<-plot_density(subset(OLs.merge1.id, idents=c('1','5','9','11','13','16','23','24','25')), 
              features = c("Bmp4","Pdgfra","Cspg4","Mobp","Bcas1","Tlcd2","Lrrc17","Plp1", "Mbp"
                          ), reduction = "umap" ,
               pal = "inferno", 
               method = c("wkde"), size = 0.5, adjust=1)

DotPlot(OLs.merge1.id, 
        group.by = "broad_type",
        idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
        # idents = c('1','9','11','13','16','23','24','25'), #,'24','25'        
        features = unique(setdiff(marker$gene,grep("^Gm|Rik$",marker$gene, value = TRUE)))[c(226:250)], assay="SCT" , cols ="RdBu") + coord_flip() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



DotPlot(OLs.merge1.id, 
        group.by = "broad_type",
        # idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
        idents = c('1','5','9','11','13','16','21','23'), #,'24','25'
        features = unique(oligo.mark), assay="SCT" , cols ="RdBu") + coord_flip() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(OLs.merge1.id, 
        group.by = "broad_type",
        split.by = "genetype",
        # idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
        idents = c('13','9'), #,'24','25'
        features = c("Cnp","Mag","Mal","Mog","Opalin","Erbb3","Myrf","Plp1","Plxnb3","Enpp6","Gpr62","Sox10","Bcas1","Prkcq","Gjb1","Gsn","Ninj2","Plekhm1","S1pr5","Ugt8a","Pla2g4a","Lingo1","S1pr4","Fa2h"
        ), assay="SCT" , cols ="RdBu") + coord_flip() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#four different gene sets for myelination with the online tool Harmonize
#https://maayanlab.cloud/Harmonizome/gene_set/myelination/GO+Biological+Process+Annotations
DotPlot(OLs.merge1.id, 
        group.by = "broad_type",
        # idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
        idents = c('1','5','9','11','13','16','21','23'), #,'24','25'
        features = c("Adam10","Adgrg6","Akt1","Arhgef10","Atrx","Bace1","Dag1","Ddx54",
                    "Dicer1","Edn2","Egr2","Ermn","Fam126a","Gh1","Gjb1","Igf1","Itga6","Itgb4",
                    "Kdr","Lingo1", "Mbp", "Med17","Mmp2","Mog","Mtmr2","Mtor","Ncam1","Nrg1",
                    "Olig2","Oprm1","Plp1","Pmp22","Pou3f1","Ptpn6","Pvrl1","Qki","Slc16a2",
                    "Smarca4","Sod1","Tcf7l2","Tyrobp","Vegfa","Yy1"
        ), assay="SCT" , cols ="RdBu") + coord_flip() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(OLs.merge1.id, 
        group.by = "broad_type", 
        # idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
                idents = c('1','5','9','11','13','16','23'), #,'24','25'
        features = unique(c("Sox10","Olig1","Oligo2","Plp1","Myrf","Mog","Mag",
                     "Slc9a3r2", "Thbs3", "Klk6", "Hopx","Fos",
                     "Ptprz1","Pdgfra","Vcan","Bmp4","Gpr17", "Rsf1", "Olfr78",
                     "Enpp6","Bcas1","Sox10","Olig2","Cspg4",
                     "Bcas1", "Rras2", "Prom1","Mbp","Mobp")), assay="SCT" , cols ="RdBu") + coord_flip() + 
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DotPlot(OLs.merge1.id, 
        idents = names(table(Idents(OLs.merge1.id)))[table(Idents(OLs.merge1.id)) > 100 ], 
        group.by = "broad_type",
                features =  c("Acta2", "Myh11", "Tagln", "Cnn1","Pln",
                    "Kcnj8", "Olfr78", "Rgs4","Rgs5", "Abcc9", "Vtn",
                    "Cxcl14","Cdh11","Wt1","Rgs2","Smoc2",
                    "Spon2","Aspg","Dpep1","Ngfr","Angptl7",
                    "Clec3b","Cd55","Fap","Cxcl16","Col14a1"), assay="SCT" , cols ="RdBu") + coord_flip() + 
                                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


ggsave(paste0(file.path("/home","choelab","working","OL_sorted","results"),"/Figure2b_Dim_20230720.jpeg"), p1a, width = 7.5, height = 7.5, units = "in", device = "jpeg", dpi = 300)
ggsave(paste0(file.path("/home","choelab","working","OL_sorted","results"),"/Figure2c_Dim_20230720.jpeg"), p1b, width = 7.5, height = 7.5, units = "in", device = "jpeg", dpi = 300)
ggsave(paste0(file.path("/home","choelab","working","OL_sorted","results"),"/Figure2d_Den_20230720.jpeg"), p1c, width = 16, height = 12, units = "in", device = "jpeg", dpi = 300)
ggsave(paste0(file.path("/home","choelab","working","OL_sorted","results"),"/SuppFig_1a_Dim_20230720.jpeg"), pS1a, width = 7.5, height = 7.5, units = "in", device = "jpeg", dpi = 300)
ggsave(paste0(file.path("/home","choelab","working","OL_sorted","results"),"/SuppFig_1b_Dot_20230720.jpeg"), pS1b, width = 15, height = 10, units = "in", device = "jpeg", dpi = 300)

require(UCell)

markers <- list("INF" = inflammation,
                "Myelination" = c("Adam10","Adgrg6","Akt1","Arhgef10","Atrx","Bace1","Dag1","Ddx54",
                    "Dicer1","Edn2","Egr2","Ermn","Fam126a","Gh1","Gjb1","Igf1","Itga6","Itgb4",
                    "Kdr","Lingo1", "Mbp", "Med17","Mmp2","Mog","Mtmr2","Mtor","Ncam1","Nrg1",
                    "Olig2","Oprm1","Plp1","Pmp22","Pou3f1","Ptpn6","Pvrl1","Qki","Slc16a2",
                    "Smarca4","Sod1","Tcf7l2","Tyrobp","Vegfa","Yy1"),
                "NFkB" = c("Chuk","Fadd","Ikbkb","Ikbkg","Il1a","Il1r1","Irak1","Map3k1","Map3k7","Map3k4","Map4k4",
                            "Myd88","Nfkb1","Nfkbia","Rela","Ripk1","Tab1","Tlr4","Tnf","Tnfrsf1a","Tradd","Traf6"
                            )) 
OLs.merge1.score <- AddModuleScore_UCell(OLs.merge1.id, features = markers, assay = "RNA")
signature.names <- paste0(names(markers), "_UCell")


require(ComplexHeatmap)
require(gtools)

# p_heat <- list()
mat <- list()

for(i in levels(Idents(OLs.merge1.score))){
  
  subdat <- c()
  meta <- c()
  TG <- c()
  KO <- c()
  subdat <- subset(OLs.merge1.score, ident = i)
  meta <- subdat@meta.data
  
  TG <- c(meta$genetype == names(table(OLs.merge1.score$genetype))[2])
  KO <- c(meta$genetype == names(table(OLs.merge1.score$genetype))[1])

  
  
  mat[[i]] <- rbind(
    mean(foldchange(meta$INF_UCell[KO][meta$INF_UCell[KO] != 0], meta$INF_UCell[TG][meta$INF_UCell[TG] != 0])),
    mean(foldchange(meta$Myelination_UCell[KO][meta$Myelination_UCell[KO] != 0], meta$Myelination_UCell[TG][meta$Myelination_UCell[TG] != 0])))
  
  rownames(mat[[i]]) <- c("Inflammation","Myelination")
  
}

arr=array(unlist(mat),dim=c(2,length(levels(Idents(OLs.merge1.score)))-1))
rownames(arr) <- c("Inflammation","Myelination")
colnames(arr) <- celltype[c(1:25)]
 celltype


val <- round(max(arr[,-13]),1)

ComplexHeatmap::pheatmap(
  arr[,c(23,9,5,1,16,11)], #25,24
  name ="Average Fold Change",
  cellwidth = 50,
  cellheight = 50,
#   fontsize = 10,
#   angle_col = 270,
  cluster_cols =F,
   breaks = c(-val,0, val),
  cluster_rows = F,
  annotation_names_col =T,
  border_color="white",
#    legend_breaks =  seq(-val,val,2*val/5), 
#    legend_labels = seq(-val,val,2*val/5),                                    
  main = "Gene Scoring : 5xFAD vs WT - O4-MACS")

pbox_INF <- list()
pbox_MYL <- list()
pbox_NFkB <- list()

for(s in c('1','5','9','11','16','23','24','25')) {

    obj <- subset(OLs.merge1.score, idents = s)
    obj$genetype[which(obj$genetype=="WT")]<-"C57BL/6"
    obj$genetype[which(obj$genetype=="TG")]<-"Tg6799"

    table(obj$broad_type)

    require(ggstatsplot)
    require(ggsignif)
    require(ggpubr)
    df <- compare_means(Inflammation ~ genetype, 
                    data = data.frame("genetype"=obj@meta.data$genetype, "Inflammation"=obj@meta.data$INF_UCell),
                    method = "t.test")  %>%
    dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
    dplyr::arrange(group1) 

    pbox_INF[[s]]<-ggbetweenstats(data.frame("genetype"=obj@meta.data$genetype, "Inflammation"=obj@meta.data$INF_UCell), x = genetype, y = Inflammation, 
                #  type = "np",
                pairwise.comparisons = FALSE,
                pairwise.display = "all", 
                    type = "p",
                    conf.level = 0.99,
                    xlab = "", ## label for the x-axis
                    ylab = "Inflammation Gene Score", ## label for the y-axis
                    package = "ggsci",
                    palette = "nrc_npg",
                plot.type = "box")+   
    geom_signif(
        comparisons = df$groups,
        map_signif_level = TRUE,
        annotations = df$p.signif,
        # y_position = c(5.5, 5.75, 6.0),
        test = NULL,
        na.rm = TRUE
    ) + 
    ggtitle(paste0(names(table(obj$broad_type)), " : ",comma(length(colnames(obj)),format = "d")," cells - 5xFAD : ", table(obj$genetype)[1], " cells, WT : ", comma(table(obj$genetype)[2],,format = "d"), " cells"))+
    theme(legend.position ="none", legend.title = element_blank(),
                                legend.text = element_text(size = 20),
                                axis.title = element_text(size = 28, face = "bold"),
                                axis.text = element_text(size = 24))+
    rotate_x_text(angle = 30)

    df <- compare_means(NFkB ~ genetype, 
                    data = data.frame("genetype"=obj@meta.data$genetype, "NFkB"=obj@meta.data$NFkB_UCell),
                    method = "t.test")  %>%
    dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
    dplyr::arrange(group1) 

    pbox_NFkB[[s]]<-ggbetweenstats(data.frame("genetype"=obj@meta.data$genetype, "NFkB"=obj@meta.data$NFkB_UCell), x = genetype, y = NFkB, 
                #  type = "np",
                pairwise.comparisons = FALSE,
                pairwise.display = "all", 
                    type = "p",
                    conf.level = 0.99,
                    xlab = "", ## label for the x-axis
                    ylab = "NF-kappaB Gene Score", ## label for the y-axis
                    package = "ggsci",
                    palette = "nrc_npg",
                plot.type = "box")+   
    geom_signif(
        comparisons = df$groups,
        map_signif_level = TRUE,
        annotations = df$p.signif,
        # y_position = c(5.5, 5.75, 6.0),
        test = NULL,
        na.rm = TRUE
    ) + 
    ggtitle(paste0(names(table(obj$broad_type)), " : ",comma(length(colnames(obj)),format = "d")," cells - 5xFAD : ", table(obj$genetype)[1], " cells, WT : ", comma(table(obj$genetype)[2],,format = "d"), " cells"))+
    theme(legend.position ="none", legend.title = element_blank(),
                                legend.text = element_text(size = 20),
                                axis.title = element_text(size = 28, face = "bold"),
                                axis.text = element_text(size = 24))+
    rotate_x_text(angle = 30)

    df <- compare_means(Myelination ~ genetype, 
                    data = data.frame("genetype"=obj@meta.data$genetype, "Myelination"=obj@meta.data$Myelination_UCell),
                    method = "t.test")  %>%
    dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
    dplyr::arrange(group1) 

    pbox_MYL[[s]]<-ggbetweenstats(data.frame("genetype"=obj@meta.data$genetype, "Myelination"=obj@meta.data$Myelination_UCell), x = genetype, y = Myelination, 
                #  type = "np",
                pairwise.comparisons = FALSE,
                pairwise.display = "all", 
                    type = "p",
                    conf.level = 0.99,
                    xlab = "", ## label for the x-axis
                    ylab = "Myelination Gene Score", ## label for the y-axis
                    package = "ggsci",
                    palette = "nrc_npg",
                plot.type = "box")+   
    geom_signif(
        comparisons = df$groups,
        map_signif_level = TRUE,
        annotations = df$p.signif,
        # y_position = c(5.5, 5.75, 6.0),
        test = NULL,
        na.rm = TRUE
    ) + 
        ggtitle(paste0(names(table(obj$broad_type)), " : ",comma(length(colnames(obj)),format = "d")," cells - 5xFAD : ", table(obj$genetype)[1], " cells, WT : ", comma(table(obj$genetype)[2],,format = "d"), " cells"))+
    theme(legend.position ="none", legend.title = element_blank(),
                                legend.text = element_text(size = 20),
                                axis.title = element_text(size = 28, face = "bold"),
                                axis.text = element_text(size = 24))+
    rotate_x_text(angle = 30)

  }

ggsave(paste0(file.path("/home","choelab","working","OL_sorted","results"),"/Figure2d_score_20230720.jpeg"), 
        pbox_MYL[[6]] + pbox_MYL[[3]] + ylab("") + pbox_MYL[[5]]+ ylab(""), width = 20, height = 7.5, units = "in", device = "jpeg", dpi = 300)

ggsave(paste0(file.path("/home","choelab","working","OL_sorted","results"),"/Figure2e_VC_20230720.jpeg"), 
        OLs.merge1.VC[[7]] + OLs.merge1.VC[[3]] + OLs.merge1.VC[[6]], width = 15, height = 10, units = "in", device = "jpeg", dpi = 300)


pbox_MYL[[6]]+pbox_NFkB[[6]] #OPC/COP
pbox_MYL[[3]]+pbox_NFkB[[3]] #NFOL
pbox_MYL[[5]]+pbox_NFkB[[5]] #MOL

pbox_MYL[[2]]+pbox_NFkB[[2]] #NFOL-to-MOL
pbox_MYL[[4]]+pbox_NFkB[[4]] 
 
pbox_MYL[[7]]+pbox_NFkB[[7]]
pbox_MYL[[8]]+pbox_NFkB[[8]]
# i<-'11'
# obj <- subset(OLs.merge1.id, idents = i)
# table(obj$genetype)
OLs.merge1.DE <- lapply(names(table(Idents(OLs.merge1.id))), function(i){

                obj <- subset(OLs.merge1.id, idents = i)
    
    if(dim(obj)[2] > 10) {
        tryCatch(
            { 
                res <- FindMarkers(obj,
                          ident.1 = names(table(obj$genetype))[1],
                          ident.2 = names(table(obj$genetype))[2], 
                          group.by = "genetype",
                          logfc.threshold = 0.00001, 
                          assay = "RNA",
                          test.use = "MAST",
                          only.pos = FALSE)},
       error = function(id) {
                    res <- c()
                 },

                        #  warning = function(w) print("Warning"),
                 finally = {


                return(res)

                 })
    } else {
        res <- c()
        return(res)
    }            

})

names(OLs.merge1.DE) <- names(table(Idents(OLs.merge1.id)))


###### pathway analysis

require(DOSE)
# library(plyr)
# library(dplyr)

require(AnnotationDbi)
#library(org.Hs.eg.db)
require(org.Mm.eg.db)
require(clusterProfiler)

Mm <- org.Mm.eg.db

wt_vs_tg6799.cc <- list()
wt_vs_tg6799.bp <- list()
wt_vs_tg6799.mf <- list()
wt_vs_tg6799.reactome <- list()
wt_vs_tg6799.kegg <- list()
wt_vs_tg6799.wiki <- list()

genenlist <- list()

for(id in c('1','5','9','11','13','16','23','24','25')) { #,'13',
#   if(is.null( neuron.05m[[id]])){
#     print(paste0(id," is NULL *****************"))
#   } else {
#     
    de <- c()
    selectedGenes <- c()
    gene_list <- c()
    gene_lists <- c()
    res.reactome <- c()
    res.wiki <- c()
    
    if(!is.null(OLs.merge1.DE[[id]])) {
    selectedGene <- OLs.merge1.DE[[id]] %>% dplyr::filter(p_val < 0.01  & abs(avg_log2FC) > 0.25)
    selectedGenes <- selectedGene[!grepl("^mt\\.", rownames(selectedGene)),]
    # selectedGenes <- hipp.DE.15m[['16']] %>% filter(p_val_adj < 0.1) %>% filter(abs(avg_log2FC) > 0.25)
    # selectedGene<-selectedGenes[-which(rownames(selectedGenes) %in% out_features(selectedGenes)),]
    
    gene_lists<- selectedGenes %>% dplyr::select(avg_log2FC)
    

    de<-AnnotationDbi::select(Mm,
                              keys = rownames(gene_lists),
                              columns = c("ENTREZID", "SYMBOL"),
                              keytype = "SYMBOL")
    
    geneLST<-data.frame("SYMBOL"=rownames(gene_lists), "avg_log2FC" = gene_lists$avg_log2FC)
    
    gene_df <- de %>% mutate(rank = rank(de$ENTREZID,  ties.method = "random")) # %>% 
      
    gene_dfs <-merge(gene_df, geneLST, by = "SYMBOL", all=TRUE) %>%   arrange(desc(rank))
    
    
    gene_dfs<-gene_dfs[!is.na(gene_dfs$ENTREZID),]
    
    if(length(gene_dfs$ENTREZID)>0) {
      
      gene.list<-gene_dfs$'avg_log2FC'
      names(gene.list)<-gene_dfs$ENTREZID
      
      
      gene.list = sort(gene.list, decreasing = TRUE)
      gene.list <- na.omit(gene.list)
      genenlist[[id]] <- gene.list
      gse.cc <- gseGO(geneList=gene.list, 
                      ont ="CC",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
      wt_vs_tg6799.cc[[id]] <- setReadable(gse.cc, 'org.Mm.eg.db', 'ENTREZID')

      gse.bp <- gseGO(geneList=gene.list, 
                      ont ="BP",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)

      wt_vs_tg6799.bp[[id]] <- setReadable(gse.bp, 'org.Mm.eg.db', 'ENTREZID')

      gse.mf <- gseGO(geneList=gene.list, 
                      ont ="MF",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
      
      wt_vs_tg6799.mf[[id]] <- setReadable(gse.mf, 'org.Mm.eg.db', 'ENTREZID')

      require(ReactomePA)
      tryCatch( 
        {
       res.reactome <- gsePathway(gene.list, 
                organism = "mouse",
                  pvalueCutoff = 0.1,
                  minGSSize = 3,
                  maxGSSize = 1000,
                  pAdjustMethod = "none", 
                  verbose = TRUE) %>% 
                  setReadable('org.Mm.eg.db', 'ENTREZID')
        },
        error = function(id) {
                res.reactome <- c()
     }, 
     finally = {

              wt_vs_tg6799.reactome[[id]]<- res.reactome

     } )

    #   tryCatch( { res.kegg <- gseKEGG(gene.list,
    #            organism     = 'mouse',
    #            keyType = "ENTREZID",
    #             minGSSize = 2,
    #               maxGSSize = 1000,
    #            pvalueCutoff = 0.1,
    #            verbose      = FALSE) 
    #   },
    #    error = function(id) {
    #             res.kegg <- c()
    #    }, finally = {
    #     bmp4ko_vs_tg6799.kegg[[id]] <- res.kegg

    #    }
    #   )
    tryCatch( 
        {
        res.wiki <- gseWP(gene.list,
               organism     = 'Mus musculus',       
                minGSSize = 3,
                  maxGSSize = 1000,
               pvalueCutoff = 0.1,
               verbose      = FALSE)
    },
        error = function(id) {
                res.wiki <- c()
     }, 
     finally = {

        wt_vs_tg6799.wiki[[id]] <- res.wiki

     } )
    }
    }
    print(id)
}

OLs.merge1.VC[[9]] #COP
OLs.merge1.VC[[7]] #OPC/COP
OLs.merge1.VC[[3]] #NFOL
OLs.merge1.VC[[2]] #between NFOL and MOL

OLs.merge1.VC[[6]] #MOL
OLs.merge1.VC[[8]] #OL-U

OLs.merge1.VC[[4]] #OL-stressed

require(ggsci)
p1e_total<- rbind(cbind(wt_vs_tg6799.bp[["23"]][grep("myelin|ensheath|membrane lipid|lipid metabolic|ipid catabolic|phosphatidylcholine metabolic|phosphorus metabolic|phosphate metabolic|immune|inflammation|NF|kappaB", wt_vs_tg6799.bp[["23"]]$Description),] %>% dplyr::filter(pvalue < 0.1), "group" = "OPC/COP"),
        cbind(wt_vs_tg6799.bp[["9"]][grep("myelin|ensheath|membrane lipid|lipid metabolic|ipid catabolic|phosphatidylcholine metabolic|phosphorus metabolic|phosphate metabolic|immune|inflammation|NF|kappaB", wt_vs_tg6799.bp[["9"]]$Description),] %>% dplyr::filter(pvalue < 0.1), "group" = "NFOL"),
        cbind(wt_vs_tg6799.bp[["16"]][grep("myelin|ensheath|membrane lipid|lipid metabolic|ipid catabolic|phosphatidylcholine metabolic|phosphorus metabolic|phosphate metabolic|immune|inflammation|NF|kappaB", wt_vs_tg6799.bp[["16"]]$Description),] %>% dplyr::filter(pvalue < 0.1), "group" = "MOL")) %>%
        ggplot(aes(NES, fct_reorder(Description, group), fill=group)) + 
        geom_col(orientation='y') + 
        scale_fill_npg() + 
        # scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
        theme_minimal() + xlab("Normalized Enrichment Score") + ylab(NULL) +
                      theme(legend.position = "bottom",
                            legend.title = element_text(size = 8, hjust = 0),
                            legend.text = element_text(size = 8),
                            #axis.title = element_text(size = 15, face = "bold"),
                            axis.text = element_text(size = 8))

ggsave(paste0(file.path("/home","choelab","working","OL_sorted","results"),"/Figure2f_GO-BP_20230720.jpeg"), 
        p1e_total, width = 12, height = 12, units = "in", device = "jpeg", dpi = 300)



wt_vs_tg6799.bp[[s]][grep("myelin|immune|inflammation|NF|kappaB", wt_vs_tg6799.bp[[s]]$Description),] %>% dplyr::filter(pvalue < 0.1)
wt_vs_tg6799.cc[[s]][grep("myelin|immune|inflammation|NF|kappaB", wt_vs_tg6799.cc[[s]]$Description),] %>% dplyr::filter(pvalue < 0.1)
wt_vs_tg6799.mf[[s]][grep("myelin|immune|inflammation|NF|kappaB", wt_vs_tg6799.mf[[s]]$Description),] %>% dplyr::filter(pvalue < 0.1)
wt_vs_tg6799.reactome[[s]][grep("myelin|immune|inflammation|NF|kB|ensheath|Myelin|Immune|Inflammation|NF|KappaB|Ensheath", wt_vs_tg6799.reactome[[s]]$Description),]  %>% dplyr::filter(pvalue < 0.1)

dotplot(wt_vs_tg6799.bp[[s]], showCategory = grep("myelin|immune|inflammation|NF|kappaB", wt_vs_tg6799.bp[[s]]$Description, value = TRUE))

wt_vs_tg6799.bp[[s]] %>% dplyr::filter(pvalue < 0.1 ) %>% view()
wt_vs_tg6799.mf[[s]] %>% dplyr::filter(pvalue < 0.1) %>% view()
wt_vs_tg6799.reactome[[s]]  %>% dplyr::filter(pvalue < 0.1) %>% view()

GeneSymbol <- list()
p1e <- list()

for(s in c('1','5','9','11','13','16','23','24','25')) { 
    x <- wt_vs_tg6799.bp[[s]]
    if(!is.null(x)) {
    y <- arrange(x, abs(NES)) %>% 
            group_by(sign(NES)) 

    GeneSymbol[[s]] <- unique(unlist(str_split(x$core_enrichment[grep("myelin|ensheath|lipid metabolic|ipid catabolic|phosphatidylcholine metabolic|phosphorus metabolic|phosphate metabolic|immune|inflammation|NF|kappaB",x$Description)],"\\/")))

    p1e[[s]] <- ggplot(y, aes(NES, fct_reorder(Description, NES), fill=p.adjust), showCategory=grep("myelin|ensheath|lipid metabolic|ipid catabolic|phosphatidylcholine metabolic|phosphorus metabolic|phosphate metabolic|immune|inflammation|NF|kappaB", wt_vs_tg6799.bp[[s]]$Description, value = TRUE)) + 
        geom_col(orientation='y') + 
        scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
        theme_minimal() + ylab(NULL)
    }
}



OLs.merge1.VC<-lapply(c('1','5','9','11','13','16','23','24','25'), function(s){

  DE <- OLs.merge1.DE[[s]]
  obj <- subset(OLs.merge1.id, idents = s)
  require(EnhancedVolcano)
  res <- EnhancedVolcano(DE,
                          lab = rownames(DE),
                          x = 'avg_log2FC',
                          y = 'p_val'   ,
                          xlim = c(-max(abs(DE$avg_log2FC))-0.1, max(abs(DE$avg_log2FC))+0.1),
                          ylim = c(-1, max(-log10(DE$p_val))+0.1),
                          title = paste0(names(table(obj$broad_type)), " : ",comma(length(colnames(obj)),format = "d")," cells - 5xFAD : ", table(obj$genetype)[1], " cells, WT : ", comma(table(obj$genetype)[2],,format = "d"), " cells"),
                          subtitle = "5xFAD versus C57BL/6 : differential expression using MAST",
                          xlab = bquote(~Log[2]~ 'fold change'),
                          ylab = bquote(~Log[10]~ 'P-value'),
                          pCutoff = 0.01,
                          FCcutoff = 0.25,
                          pointSize = 0.01,
                          labSize = 4,
                          selectLab = unique(c(GeneSymbol[[s]],"Bmp4","Bcas1","Fos","Mbp","Plp1","Mag","Mog","Cnp","Nfkb1")),
                          # shapeCustom = keyvals.shape,
                          #colCustom = keyvals.colour,
                          # col=c("grey30", "forestgreen", "royalblue", "red2"),
                          #                colCustom = keyvals,
                          colAlpha = 0.5,
                          legendPosition = 'bottom',
                          drawConnectors = TRUE,
                          widthConnectors = 0.05,
                          colConnectors = 'black',
                          border = 'full',
                          borderWidth = 1.0,
                          # labCol = "black",
                          # legendPosition = 'none',
                          axisLabSize = 10,
                          titleLabSize = 12,
                          subtitleLabSize = 10,
                          captionLabSize = 10,
                          maxoverlapsConnectors = Inf)

  return(res)

})

names(OLs.merge1.VC) <- c('1','5','9','11','13','16','23','24','25')





NFOL.DE <- FindMarkers(OLs.merge1.id,
                          ident.1 = '9',
                          ident.2 = '13', 
                        #   group.by = "genetype",
                          logfc.threshold = 0.00001, 
                          assay = "RNA",
                          test.use = "MAST",
                          only.pos = FALSE)



require(EnhancedVolcano)
EnhancedVolcano(NFOL.DE,
                          lab = rownames(NFOL.DE),
                          x = 'avg_log2FC',
                          y = 'p_val'   ,
                        #   xlim = c(-max(abs(NFOL.DE$avg_log2FC))-0.1, max(abs(NFOL.DE$avg_log2FC))+0.1),
                        #   ylim = c(-1, max(-log10(NFOL.DE$p_val))+0.1),
                          xlim = c(-0.5, 0.5),
                          ylim = c(-1, 10),
                        #   title = paste0(names(table(obj$broad_type)), " : ",comma(length(colnames(obj)),format = "d")," cells - 5xFAD : ", table(obj$genetype)[1], " cells, WT : ", comma(table(obj$genetype)[2],,format = "d"), " cells"),
                        #   subtitle = "5xFAD versus C57BL/6 : differential expression using MAST",
                          xlab = bquote(~Log[2]~ 'fold change'),
                          ylab = bquote(~Log[10]~ 'P-value'),
                          pCutoff = 0.01,
                          FCcutoff = 0.25,
                          pointSize = 0.01,
                          labSize = 4,
                          selectLab = c("Adam10","Adgrg6","Akt1","Arhgef10","Atrx","Bace1","Dag1","Ddx54",
                    "Dicer1","Edn2","Egr2","Ermn","Fam126a","Gh1","Gjb1","Igf1","Itga6","Itgb4",
                    "Kdr","Lingo1", "Mbp", "Med17","Mmp2","Mog","Mtmr2","Mtor","Ncam1","Nrg1",
                    "Olig2","Oprm1","Plp1","Pmp22","Pou3f1","Ptpn6","Pvrl1","Qki","Slc16a2",
                    "Smarca4","Sod1","Tcf7l2","Tyrobp","Vegfa","Yy1"),#c("Bmp4","Bcas1","Fos","Mbp","Plp1","Mag","Mog","Cnp","Nfkb1",intersect(Infla, NFOL.DE %>% dplyr::filter(p_val < 0.01 & abs(avg_log2FC) > 0.25) %>% rownames()),"Cnp","Mag","Mal","Mog","Opalin","Erbb3","Myrf","Plp1","Plxnb3","Enpp6","Gpr62","Sox10","Bcas1","Prkcq","Gjb1","Gsn","Ninj2","Plekhm1","S1pr5","Ugt8a","Pla2g4a","Lingo1","S1pr4","Fa2h"),#unique(c()),
                          # shapeCustom = keyvals.shape,
                          #colCustom = keyvals.colour,
                          # col=c("grey30", "forestgreen", "royalblue", "red2"),
                          #                colCustom = keyvals,
                          colAlpha = 0.5,
                          legendPosition = 'bottom',
                          drawConnectors = TRUE,
                          widthConnectors = 0.05,
                          colConnectors = 'black',
                          border = 'full',
                          borderWidth = 1.0,
                          # labCol = "black",
                          # legendPosition = 'none',
                          axisLabSize = 10,
                          titleLabSize = 12,
                          subtitleLabSize = 10,
                          captionLabSize = 10,
                          maxoverlapsConnectors = Inf)

    selectedGene <- NFOL.DE %>% dplyr::filter(p_val < 0.01  & abs(avg_log2FC) > 0.25)
    selectedGenes <- selectedGene[!grepl("^mt\\.", rownames(selectedGene)),]
    # selectedGenes <- hipp.DE.15m[['16']] %>% filter(p_val_adj < 0.1) %>% filter(abs(avg_log2FC) > 0.25)
    # selectedGene<-selectedGenes[-which(rownames(selectedGenes) %in% out_features(selectedGenes)),]
    
    gene_lists<- selectedGenes %>% dplyr::select(avg_log2FC)
    

    de<-AnnotationDbi::select(Mm,
                              keys = rownames(gene_lists),
                              columns = c("ENTREZID", "SYMBOL"),
                              keytype = "SYMBOL")
    
    geneLST<-data.frame("SYMBOL"=rownames(gene_lists), "avg_log2FC" = gene_lists$avg_log2FC)
    
    gene_df <- de %>% mutate(rank = rank(de$ENTREZID,  ties.method = "random")) # %>%       
    gene_dfs <-merge(gene_df, geneLST, by = "SYMBOL", all=TRUE) %>%   arrange(desc(rank))
    gene_dfs<-gene_dfs[!is.na(gene_dfs$ENTREZID),]
    
    
      
      gene.list<-gene_dfs$'avg_log2FC'
      names(gene.list)<-gene_dfs$ENTREZID
      
      
      gene.list = sort(gene.list, decreasing = TRUE)
      gene.list <- na.omit(gene.list)
      gse.cc <- gseGO(geneList=gene.list, 
                      ont ="CC",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
      NFOL.cc <- setReadable(gse.cc, 'org.Mm.eg.db', 'ENTREZID')

      gse.bp <- gseGO(geneList=gene.list, 
                      ont ="BP",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)

      NFOL.bp <- setReadable(gse.bp, 'org.Mm.eg.db', 'ENTREZID')

      gse.mf <- gseGO(geneList=gene.list, 
                      ont ="MF",
                      keyType = "ENTREZID",
                      #             nPerm = 10000,
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
      
      NFOL.mf <- setReadable(gse.mf, 'org.Mm.eg.db', 'ENTREZID')

    NFOL.reactome <- gsePathway(gene.list, 
                organism = "mouse",
                  pvalueCutoff = 0.1,
                  minGSSize = 3,
                  maxGSSize = 1000,
                  pAdjustMethod = "none", 
                  verbose = TRUE)

    NFOL.reactome <- setReadable(NFOL.reactome, 'org.Mm.eg.db', 'ENTREZID')

    NFOL.cc[grep("myelin|immune|inflammation|NF|kappaB|ensodome", NFOL.cc$Description),] %>% dplyr::filter(pvalue < 0.1 ) %>% view()
    NFOL.bp[grep("myelin|immune|inflammation|NF|kappaB|ensodome|enshea|lipid", NFOL.bp$Description),] %>% dplyr::filter(pvalue < 0.1 ) %>% view()
    NFOL.mf[grep("myelin|immune|inflammation|NF|kappaB|ensodome|ensheath", NFOL.mf$Description),] %>% dplyr::filter(pvalue < 0.1 ) %>% view()
    NFOL.reactome[grep("myelin|immune|inflammation|NF|kappaB|Myelin|Immune|Inflammation|NF|KappaB", NFOL.reactome$Description),]  %>% dplyr::filter(pvalue < 0.1)

    NFOL.cc %>% dplyr::filter(pvalue < 0.1 & NES < 0 ) %>% view()
    NFOL.mf %>% dplyr::filter(pvalue < 0.1 & NES < 0 ) %>% view()
    NFOL.bp %>% dplyr::filter(pvalue < 0.1 & NES < 0 ) %>% view()
    NFOL.reactome %>% dplyr::filter(pvalue < 0.1 & NES < 0 ) %>% view()

# image_name <- paste0(file.path("/home","choelab","working","OL_sorted","final_Rdata"),"/OL_sorted_",format(Sys.time(), "%Y-%m-%d %H-%M-%S"), ".Rdata")
# save.image(image_name)

#load("/home/choelab/working/OL_sorted/final_Rdata/OL_sorted_2023-07-20 10-36-14.Rdata")
