# ==============================================================================
# Title : Oligodendrocyte precusor cells-microglia cross talk via BMP4 signaling drives microglial neuroprotective response and mitigates Alzheimer's diseas progression
# Author: Jaemyung Jang / Korea Brain Research Institute (piloter2@kbri.re.kr)
# Description: This script performs visualization, differential expression analysis, and 
#   clustering annotation for oligodendrocyte lineage cells (OPC, COP, OL)
# Date: 2024-12-07
# ==============================================================================

# ==============================================================================
# 1. Setup & Configuration
# ==============================================================================

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsci)
library(xlsx)
library(DOSE)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)

cor.71585.merge.genetype.sub <- readRDS("/path/to/bmp4ko_workingOBJ_20230703.RDS")

# ==============================================================================
# 2. Annotation & Cell Type Assignment
# ==============================================================================

Idents(cor.71585.merge.genetype.sub) <- "RNA_snn_res.1.4"
cell_type_map <- c(
  "1"  = "Neuron 1",
  "2"  = "Fib",
  "3"  = "MOL 1",
  "4"  = "COP 4",
  "5"  = "EPC 1",
  "6"  = "Microglia",
  "7"  = "EPC 2",
  "8"  = "COP 3",
  "9"  = "DOL",
  "10" = "Neuron 2",
  "11" = "Endo",
  "12" = "ASC",
  "13" = "CP",
  "14" = "COP 1",
  "15" = "OPC",
  "16" = "COP 2",
  "17" = "MOL 2",
  "18" = "MOL 3",
  "19" = "unknown",  
  "20" = "NFOL 1",
  "21" = "unknown 2",
  "22" = "Neuron 4",
  "23" = "unknown",
  "24" = "NFOL 2",
  "25" = "COP 5",
  "26" = "26",
  "27" = "Neuron 3",
  "28" = "DOL 2",
  "29" = "DOL 3",
  "30" = "DOL 4",
  "31" = "MOL 4",  
  "32" = "32",
  "33" = "33",
  "34" = "NFOL 2",
  "35" = "35",
  "36" = "36",
  "37" = "37"
)

# 3. Assign the 'broad_type' metadata column efficiently
# Map the current Idents (cluster IDs) to the cell types defined above
# This avoids using a slow 'for' loop
current_ids <- as.character(Idents(cor.71585.merge.genetype.sub))
cor.71585.merge.genetype.sub$broad_type <- cell_type_map[current_ids]

for(i in names(celltype)){
  cor.71585.merge.genetype.sub@meta.data$broad_type[colnames(cor.71585.merge.genetype.sub) %in% WhichCells(subset(cor.71585.merge.genetype.sub, ident = i))] <- celltype[i]
}

Idents(cor.71585.merge.genetype.sub) <- as.factor(cor.71585.merge.genetype.sub@meta.data$broad_type)


# ==============================================================================
# 3. Visualization (DimPlot & DotPlot)
# ==============================================================================

# Figure2i
p2i <- DimPlot(subset(cor.71585.merge.genetype.sub, idents = names(table(Idents(cor.71585.merge.genetype.sub)))[table(Idents(cor.71585.merge.genetype.sub))>150]),
             group.by = "broad_type", 
             cols = c(pal_npg("nrc")(10)[c(2:6)],pal_npg("nrc")(10)[c(7,9)],pal_aaas("default")(10)[c(1:8)],pal_nejm("default")(8)[c(1:3)],pal_nejm("default")(8)[c(4:6)],rep("grey90",2)), 
             label =T, label.size = 4, pt.size = 0.01,raster=FALSE)+
             theme(legend.position = "none",  axis.title = element_blank(), 
                   plot.title = element_text(size = 12, face = "bold",hjust = 0), axis.line = element_blank(), 
                   axis.text = element_blank(),  axis.ticks = element_blank()) + 
                   ggtitle("") + NoAxes()

cor.71585.merge.genetype.sub$new_ident <- factor(cor.71585.merge.genetype.sub$broad_type, levels = names(table(cor.71585.merge.genetype.sub$broad_type))[c(1:7,13,18:21,27:30,33,8:12,31:32,23:26,14:17,22,34:35)])

# Supplementary Figure6a
sp6a <- DotPlot(cor.71585.merge.genetype.sub, 
            idents = names(table(Idents(cor.71585.merge.genetype.sub)))[table(Idents(cor.71585.merge.genetype.sub))>150],
            features = unique(c(grep("Pdgfa",rownames(cor.71585.merge.genetype.sub), value = TRUE), "Slc1a2", "Enpp2","Flt1","Epyc","Col1a2","Meg3","Snap25","Ptprz1","Bcas1","Apod","Trf","Cnp","Mobp","Mal","Cryab","Ndrg1","Cnp","Car2","Cst3","Ctss","Dcn")), 
            assay="SCT", cols = c("RdBu"), dot.scale = 8) + 
            RotatedAxis() + coord_flip()+
            theme(legend.title = element_text(size = 12, hjust = 0),
                  legend.text = element_text(size = 12),
                  axis.title = element_text(size = 15, face = "bold"),
                  axis.text = element_text(size = 18))+
            xlab("")+ylab("")

# ==============================================================================
# 4. Differential Expression Analysis (DE)
# ==============================================================================

cor.71585.merge.genetype.sub.DE <- list()

for(id in 1:length(names(table(cor.71585.merge.genetype.sub$broad_type)))) {
        
        res <- c()
        obj <- c()
        obj <- subset(cor.71585.merge.genetype.sub, subset = broad_type == names(table(cor.71585.merge.genetype.sub$broad_type))[id])

        if(dim(obj)[2] > 10) {
        tryCatch(
            { 
              res <- FindMarkers(obj,
                          ident.1 = names(table(cor.71585.merge.genetype.sub$genetype))[2],
                          ident.2 = names(table(cor.71585.merge.genetype.sub$genetype))[1], 
                          group.by = "genetype",
                          logfc.threshold = 0.0001, 
                          assay = "SCT",
                          test.use = "MAST",
                          only.pos = FALSE)},
                error = function(id) {
                    res <- FindMarkers(obj,
                          ident.1 = names(table(cor.71585.merge.genetype.sub$genetype))[2],
                          ident.2 = names(table(cor.71585.merge.genetype.sub$genetype))[1], 
                          group.by = "genetype",
                          logfc.threshold = 0.0001, 
                          assay = "SCT",
                          only.pos = FALSE)
                 },
                finally = {
                cor.71585.merge.genetype.sub.DE[[id]] <- res
                 })
                 } else {
                  cor.71585.merge.genetype.sub.DE[[id]] <- c()
                  }
            # write.csv(res, paste0("/path/to/results/cluster_DE_",names(table(cor.71585.merge.genetype.sub$broad_type))[id],"_",format(Sys.time(), "%Y-%m-%d %H-%M-%S"),".xlsx"))
}

names(cor.71585.merge.genetype.sub.DE) <- names(table(cor.71585.merge.genetype.sub$broad_type))

# ==============================================================================
# 5. Plotting DE Numbers
# ==============================================================================
dfs <- c()
for(i in c(33, 8:12, 14:17, 24:26, 31:32)){
  res <- cor.71585.merge.genetype.sub.DE[[i]] %>% dplyr::filter(p_val < 0.05 & avg_log2FC < -0.1) %>% rownames() %>% length()  %>% as.numeric()
  df<-rbind(data.frame("num"= cor.71585.merge.genetype.sub.DE[[i]] %>% dplyr::filter(p_val < 0.05 & avg_log2FC > 0.1) %>% rownames() %>% length() %>% as.numeric(), "type" = "upregulated", "cell" = names(cor.71585.merge.genetype.sub.DE)[i]), 
      data.frame("num"= 0-res, "type" = "downregulated", "cell" = names(cor.71585.merge.genetype.sub.DE)[i]))
  dfs <- rbind(dfs, df)
}

# 경로 수정: Figure 데이터 저장 (주석)
# write.csv(dfs, paste0("/path/to/results/Figure3c_DE_num_","_",format(Sys.time(), "%Y-%m-%d %H-%M-%S"),".xlsx"))

# Supplementary Figure 6b
sp6b <- dfs %>% 
   mutate(cell = factor(cell, levels = unique(dfs$cell[c(grep("OPC",dfs$cell), grep("COP",dfs$cell), grep("NFOL",dfs$cell), grep("MOL",dfs$cell), grep("DOL",dfs$cell))]))) %>% 
  ggplot(aes(x=num, y=cell, fill=cell)) +
  geom_bar(position="stack", stat="identity", width = 0.5)+
    theme_bw()+
  theme(legend.position = "none", legend.title = element_blank(),
                            legend.text = element_text(size = 20),
                            axis.title = element_text(size = 28, face = "bold"),
                            axis.text = element_text(size = 24))+
  scale_fill_manual(values = c(brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"),brewer.pal(12,"Set3"),brewer.pal(1,"Dark2"))[c(27, 2:6, 25:26, 18:20, 8:12)])+
  ggtitle("Number of DEGs (p<0.05, log2_avg|FC| > 0.1)")+
  xlab("")+
  ylab("")+
  RotatedAxis()

# ==============================================================================
# 6. Gene Set Enrichment Analysis (GSEA) - GO, Reactome, WikiPathways
# ==============================================================================

Mm <- org.Mm.eg.db

bmp4ko_vs_tg6799.cc <- list()
bmp4ko_vs_tg6799.bp <- list()
bmp4ko_vs_tg6799.mf <- list()
bmp4ko_vs_tg6799.reactome <- list()
bmp4ko_vs_tg6799.kegg <- list()
bmp4ko_vs_tg6799.wiki <- list()

genenlist <- list()

for(id in names(cor.71585.merge.genetype.sub.DE)) {
    de <- c()
    selectedGenes <- c()
    gene_list <- c()
    gene_lists <- c()
    res.reactome <- c()
    res.wiki <- c()
    
    if(!is.null(cor.71585.merge.genetype.sub.DE[[id]])) {
    selectedGene <- cor.71585.merge.genetype.sub.DE[[id]] %>% dplyr::filter(p_val < 0.01  & abs(avg_log2FC) > 0.25)
    selectedGenes <- selectedGene[!grepl("^mt\\.", rownames(selectedGene)),]
    
    gene_lists<- selectedGenes %>% dplyr::select(avg_log2FC)
    
    de<-AnnotationDbi::select(Mm,
                              keys = rownames(gene_lists),
                              columns = c("ENTREZID", "SYMBOL"),
                              keytype = "SYMBOL")
    
    geneLST<-data.frame("SYMBOL"=rownames(gene_lists), "avg_log2FC" = gene_lists$avg_log2FC)
    
    gene_df <- de %>% mutate(rank = rank(de$ENTREZID,  ties.method = "random")) 
      
    gene_dfs <-merge(gene_df, geneLST, by = "SYMBOL", all=TRUE) %>%   arrange(desc(rank))
    
    gene_dfs<-gene_dfs[!is.na(gene_dfs$ENTREZID),]
    
    if(length(gene_dfs$ENTREZID)>0) {
      
      gene.list<-gene_dfs$'avg_log2FC'
      names(gene.list)<-gene_dfs$ENTREZID
      
      gene.list = sort(gene.list, decreasing = TRUE)
      gene.list <- na.omit(gene.list)
      genenlist[[id]] <- gene.list
      
      # GO: CC
      gse.cc <- gseGO(geneList=gene.list, 
                      ont ="CC",
                      keyType = "ENTREZID",
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
      bmp4ko_vs_tg6799.cc[[id]] <- setReadable(gse.cc, 'org.Mm.eg.db', 'ENTREZID')

      # GO: BP
      gse.bp <- gseGO(geneList=gene.list, 
                      ont ="BP",
                      keyType = "ENTREZID",
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
      bmp4ko_vs_tg6799.bp[[id]] <- setReadable(gse.bp, 'org.Mm.eg.db', 'ENTREZID')

      # GO: MF
      gse.mf <- gseGO(geneList=gene.list, 
                      ont ="MF",
                      keyType = "ENTREZID",
                      minGSSize = 3,
                      maxGSSize = 1000,
                      pvalueCutoff = 0.1,
                      pAdjustMethod = "none",
                      verbose = TRUE,
                      OrgDb = Mm)
      bmp4ko_vs_tg6799.mf[[id]] <- setReadable(gse.mf, 'org.Mm.eg.db', 'ENTREZID')

      # Reactome
      require(ReactomePA)
      tryCatch( 
        {
       res.reactome <- gsePathway(gene.list, 
                organism = "mouse",
                  pvalueCutoff = 0.1,
                  minGSSize = 3,
                  maxGSSize = 1000,
                  pAdjustMethod = "none", 
                  verbose = TRUE)
        },
        error = function(id) {
                res.reactome <- c()
       }, 
      finally = {
              bmp4ko_vs_tg6799.reactome[[id]]<- res.reactome
      } )

      # WikiPathways
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
               bmp4ko_vs_tg6799.wiki[[id]] <- res.wiki
      } )
    }
    }
    print(id)
}

# ==============================================================================
# 7. Save Results to Excel
# ==============================================================================

wb <- createWorkbook()
# Note: 'cor.71585.merge.genetype.sub.marker2' must be defined prior to this step
lapply(names(table(cor.71585.merge.genetype.sub$broad_type)), function (id) {
  sheet <- createSheet(wb, id)
  addDataFrame(cor.71585.merge.genetype.sub.marker2 %>% dplyr::filter(cluster == id), sheet=sheet, startColumn=1, row.names=FALSE)
})
saveWorkbook(wb, "/path/to/results/Supplementary_file_1.xlsx")


wb <- createWorkbook()
for(id in names(table(cor.71585.merge.genetype.sub$broad_type))) {
  sheet <- createSheet(wb, id)
  addDataFrame(cor.71585.merge.genetype.sub.DE[[id]], sheet=sheet, startColumn=1, row.names=TRUE)
}
saveWorkbook(wb, "/path/to/results/Supplementary_file_2.xlsx")


wb <- createWorkbook()
lapply(names(table(cor.71585.merge.genetype.sub$broad_type)), function (id) {
  sheet <- createSheet(wb, id)
  if(dim(data.frame(bmp4ko_vs_tg6799.cc[[id]]))[1] != 0 & dim(data.frame(bmp4ko_vs_tg6799.cc[[id]]))[2] != 0 ){
  addDataFrame(bmp4ko_vs_tg6799.cc[[id]] %>% data.frame(), sheet=sheet, startColumn=1, row.names=FALSE)
  } 
})
saveWorkbook(wb, "/path/to/results/Supplementary_file_3_GO-CC.xlsx")

wb <- createWorkbook()
lapply(names(table(cor.71585.merge.genetype.sub$broad_type)), function (id) {
  sheet <- createSheet(wb, id)
  if(dim(data.frame(bmp4ko_vs_tg6799.bp[[id]]))[1] != 0 & dim(data.frame(bmp4ko_vs_tg6799.bp[[id]]))[2] != 0 ){
  addDataFrame(bmp4ko_vs_tg6799.bp[[id]] %>% data.frame(), sheet=sheet, startColumn=1, row.names=FALSE)
  }
})
saveWorkbook(wb, "/path/to/results/Supplementary_file_3_GO-BP.xlsx")

wb <- createWorkbook()
lapply(names(table(cor.71585.merge.genetype.sub$broad_type)), function (id) {
  sheet <- createSheet(wb, id)
  if(dim(data.frame(bmp4ko_vs_tg6799.mf[[id]]))[1] != 0 & dim(data.frame(bmp4ko_vs_tg6799.mf[[id]]))[2] != 0 ){
  addDataFrame(bmp4ko_vs_tg6799.mf[[id]] %>% data.frame(), sheet=sheet, startColumn=1, row.names=FALSE)
  }
})
saveWorkbook(wb, "/path/to/results/Supplementary_file_3_GO-MF.xlsx")

wb <- createWorkbook()
lapply(names(table(cor.71585.merge.genetype.sub$broad_type)), function (id) {
  sheet <- createSheet(wb, id)
  if(dim(data.frame(bmp4ko_vs_tg6799.reactome[[id]]))[1] != 0 & dim(data.frame(bmp4ko_vs_tg6799.reactome[[id]]))[2] != 0 ){
  addDataFrame(bmp4ko_vs_tg6799.reactome[[id]] %>% setReadable( 'org.Mm.eg.db', 'ENTREZID')%>% data.frame(), sheet=sheet, startColumn=1, row.names=FALSE)
  }
})
saveWorkbook(wb, "/path/to/results/Supplementary_file_3_reactome.xlsx")


wb <- createWorkbook()
lapply(names(table(cor.71585.merge.genetype.sub$broad_type)), function (id) {
  sheet <- createSheet(wb, id)
  if(dim(data.frame(bmp4ko_vs_tg6799.wiki[[id]]))[1] != 0 & dim(data.frame(bmp4ko_vs_tg6799.wiki[[id]]))[2] != 0 ){
  addDataFrame(bmp4ko_vs_tg6799.wiki[[id]] %>% data.frame(), sheet=sheet, startColumn=1, row.names=FALSE)
  }
})
saveWorkbook(wb, "/path/to/results/Supplementary_file_3_wiki.xlsx")
