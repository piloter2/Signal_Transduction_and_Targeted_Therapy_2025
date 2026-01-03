# ==============================================================================
# Title : Oligodendrocyte precusor cells-microglia cross talk via BMP4 signaling drives microglial neuroprotective response and mitigates Alzheimer's diseas progression
# Author: Jaemyung Jang / Korea Brain Research Institute (piloter2@kbri.re.kr)
# Description: This script performs visualization, differential expression analysis, and 
#   clustering annotation for oligodendrocyte lineage cells (OPC, COP, OL)
# Date: 2025-12-08
# ==============================================================================

# ==============================================================================
# 1. Setup & Configuration
# ==============================================================================

# --- Libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggsci)
  library(tibble)
  library(pheatmap)
  library(ComplexHeatmap)
  library(EnhancedVolcano)
  library(openxlsx)
  library(Seurat)
  library(tidyverse)
  library(cowplot); library(ggpubr)
  library(RColorBrewer)
})

# --- Parameters & Paths ---
# NOTE: Users should adjust these paths to match their local environment.
BASE_DIR    <- "/path/to/"
RESULTS_DIR <- file.path(BASE_DIR, "results")
DATE_SUFFIX <- "20251208"

# Create output directory if it doesn't exist
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# Set seed for reproducibility
set.seed(1234)

# ==============================================================================
# 2. Visualization Functions
# ==============================================================================

#' Wrapper for Seurat DimPlot with customized theme
#' @param obj Seurat object
#' @param group_by Metadata column to group by
#' @param reduction Reduction method (e.g., 'umap')
#' @param label Boolean, whether to label clusters
#' @param title_prefix Prefix for the plot title
#' @param subset_ids Optional: Vector of IDs to subset the object before plotting
plot_custom_dimplot <- function(obj, group_by, reduction = "umap", label = TRUE, title_prefix = "", subset_ids = NULL) {
  
  plot_obj <- obj
  if (!is.null(subset_ids)) {
    plot_obj <- subset(obj, subset = id %in% subset_ids)
  }
  
  n_cells <- comma(length(colnames(plot_obj)), format = "d")
  plot_title <- paste0(title_prefix, " - cells : ", n_cells, " cells")
  
  p <- DimPlot(plot_obj, 
               cols = pal_npg("nrc")(10),
               group.by = group_by, 
               reduction = reduction, 
               label = label, repel = label, 
               label.size = 6, pt.size = 0.01) + 
    theme(axis.title = element_blank(), 
          plot.title = element_text(size = 14, face = "bold", hjust = 0), 
          axis.line = element_blank(), 
          axis.text = element_blank(),  
          axis.ticks = element_blank()) + 
    ggtitle(plot_title)
  
  if (!label) {
    p <- p + NoLegend()
  }
  
  return(p)
}

ref.merge.oligo <- readRDS(paste0(BASE_DIR,"/data/ref.merge.oligo.rds"))

# Generate Plots : Supplementary Figure 2a
sp2a <- plot_custom_dimplot(ref.merge.oligo, group_by = "ref.celltype", label = FALSE)

# ==============================================================================
# 3. Differential Expression (COP Markers)
# ==============================================================================

# Subset for PS2APP analysis
PS2APP <- subset(ref.merge.oligo, subset = id %in% target_subset_ids)
PS2APP <- PrepSCTFindMarkers(PS2APP)

# Find Markers
COP.markers <- FindMarkers(PS2APP,
                           group.by = "celltype", 
                           test.use = "MAST",
                           ident.1 = "COP")

# Save Marker Table
write.xlsx(COP.markers, file.path(RESULTS_DIR, paste0("SuppFigure2a_", DATE_SUFFIX, ".xlsx")), rowNames = TRUE)

# Volcano Plot
# Filter significant genes for labeling
selGenes <- rownames(COP.markers %>% dplyr::filter(p_val_adj < 1e-300 & avg_log2FC > 2))
clean_labels <- setdiff(selGenes, grep("^Gm|Rik$", selGenes, value = TRUE))
specific_labels <- c("Bcas1", "Gpr17", "Enpp6", "Bmp4")

# Generate Plots : Supplementary Figure 2b
EnhancedVolcano(COP.markers,
                lab = rownames(COP.markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "COP vs OPC/Oligodendrocytes",
                subtitle = paste0("COP: ", comma(sum(PS2APP$celltype == "COP")), 
                                  " cells - Others: ", comma(sum(PS2APP$celltype != "COP")), " cells"),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~Log[10]~ 'Adjust P-value'),
                pCutoff = 1e-300, FCcutoff = 2,
                pointSize = 0.25, labSize = 3.2,
                selectLab = specific_labels,
                colAlpha = 0.5, legendPosition = 'bottom',
                drawConnectors = TRUE, widthConnectors = 0.05, colConnectors = 'black',
                border = 'full', borderWidth = 0.25,
                axisLabSize = 10, titleLabSize = 12, subtitleLabSize = 10, captionLabSize = 10,
                maxoverlapsConnectors = 50)

# ==============================================================================
# 4. Comparison Analysis (Heatmap across Datasets)
# ==============================================================================

# Helper function for DE
run_de <- function(obj, ident1, ident2, group_col = "genetype") {
  FindMarkers(obj, ident.1 = ident1, ident.2 = ident2, group.by = group_col, 
              test.use = "MAST", recorrect_umi = FALSE)
}

# Run Comparisons
# Note: Using sorted unique genotypes to ensure correct index matching
# (Adjust indices [12], [10] etc. based on your specific sorted genotype list)
u_genotypes <- sort(unique(OL.COP$genetype))

OL.COP.DEP.DM_with_CupRap_0wks <- run_de(OL.COP, u_genotypes[12], u_genotypes[10]) # Verify indices!
OL.COP.DEP.DM_with_CupRap_3wks <- run_de(OL.COP, u_genotypes[13], u_genotypes[11])
OL.COP.DEP.DM_with_CPZ3        <- run_de(OL.COP, u_genotypes[1],  u_genotypes[9])
OL.COP.DEP.PS2_APP             <- run_de(OL.COP, u_genotypes[5],  u_genotypes[8])
OL.COP.DEP.tau                 <- run_de(OL.COP, u_genotypes[6],  u_genotypes[8])

# Save DE Results
write.xlsx(OL.COP.DEP.DM_with_CPZ3, file.path(RESULTS_DIR, paste0("SuppFigure2b_cpz3_", DATE_SUFFIX, ".xlsx")), rowNames = TRUE)
write.xlsx(OL.COP.DEP.PS2_APP, file.path(RESULTS_DIR, paste0("SuppFigure2b_ps2app_", DATE_SUFFIX, ".xlsx")), rowNames = TRUE)
write.xlsx(OL.COP.DEP.tau, file.path(RESULTS_DIR, paste0("SuppFigure2b_tau_", DATE_SUFFIX, ".xlsx")), rowNames = TRUE)

# Prepare Heatmap Data
target_genes <- c("Bmp4", "Enpp6", "Bcas1", "Gpr17", "Cnp", "Cyfip2", 
                  "Lims2", "Wasf1", "Phyhipl", "Nfasc", "Frmd4a", "Olig2", "Pdgfa")

extract_logFC <- function(res, col_name) {
  res %>% rownames_to_column("Gene") %>% 
    dplyr::select(Gene, avg_log2FC) %>% 
    dplyr::rename(!!col_name := avg_log2FC)
}

df_cpz <- extract_logFC(OL.COP.DEP.DM_with_CPZ3, "CPZ")
df_tau <- extract_logFC(OL.COP.DEP.tau, "P301L")
df_ps2 <- extract_logFC(OL.COP.DEP.PS2_APP, "PS2APP")

heatmap_data <- data.frame(Gene = target_genes) %>%
  left_join(df_cpz, by = "Gene") %>%
  left_join(df_tau, by = "Gene") %>%
  left_join(df_ps2, by = "Gene")

rownames(heatmap_data) <- heatmap_data$Gene
heatmap_data <- heatmap_data[, -1]
heatmap_data[is.na(heatmap_data)] <- 0

# Draw Heatmap : Supplementary Figure 2c
pheatmap(heatmap_data,
         name = "Average Fold Change",
         cellwidth = 25, cellheight = 25,
         cluster_cols = FALSE, cluster_rows = FALSE,
         border_color = "white",
         main = "Average expression")

# Final Check UMAPs
p_check_cop <- DimPlot(OL.COP, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.01) + 
  theme(legend.position = "none") + NoAxes() + ggtitle("")
p_check_merged <- DimPlot(OL.merged, reduction = "umap", label = FALSE, pt.size = 0.01) + 
  theme(legend.position = "none") + NoAxes() + 
  ggtitle(paste0(comma(ncol(OL.merged), format = "d"), " cells"))

# print(p_check_cop)
# print(p_check_merged)
# ==============================================================================
# 5. Visualization Functions
# ==============================================================================
ref.merge <- readRDS(file=paste0("/path/to","/data/ref.merge.rds"))
ref.merge <- UpdateSeuratObject(ref.merge)

if (!"SCT_snn_res.0.3" %in% colnames(ref.merge@meta.data)) {
  ref.merge <- FindClusters(ref.merge, resolution = 0.3, algorithm = 2)
}

celltype_map <- c(
  "0" = "OLG", "1" = "MG", "2" = "MG", "3" = "OPC",
  "4" = "OLG", "5" = "OLG", "6" = "ASC", "7" = "BMEC",
  "8" = "pericytes", "9" = "EN", "10" = "Choroid", "11" = "Reelin",
  "12" = "OLG", "13" = "OLG", "14" = "TCell", "15" = "VSMC",
  "16" = "immune_a", "17" = "Ependymal", "18" = "microOPC", "19" = "OLG",
  "20" = "microOLG", "21" = "NPC", "22" = "C15", "23" = "COP",
  "24" = "microAstro", "25" = "VLMC", "26" = "immune_b", "27" = "OLG",
  "28" = "Choroid"
)

Idents(ref.merge) <- "SCT_snn_res.0.3"
ref.merge <- RenameIdents(ref.merge, celltype_map)
ref.merge$ref.celltype <- Idents(ref.merge)

# Define Order for Plots
plot_order <- c(
  "ASC", "OPC", "COP", "OLG", # Glia / Oligo lineage
  "NPC", "EN", "Reelin",      # Neuronal
  "BMEC", "pericytes", "VLMC", "VSMC", # Vascular
  "MG", "TCell", "immune_a", "immune_b", # Immune
  "Choroid", "Ependymal",     # Others
  "microAstro", "microOPC", "microOLG", "C15" # Micro-clusters
)

cat("Generating UMAP...\n")
# Combine palettes from ggsci
custom_pal <- c(
  pal_aaas("default", alpha = 0.6)(10)[1:8], 
  pal_npg("nrc", alpha = 0.6)(10), 
  pal_nejm("default", alpha = 0.6)(8)
)

# Supplementary Figure 2i
require(ggsci)
p_umap <- DimPlot(ref.merge, label = TRUE, label.size = 6, pt.size = 0.01, raster = FALSE, cols = custom_pal) + 
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        plot.title = element_text(size = 12, face = "bold", hjust = 0),
        axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()) + 
  NoAxes()


### 3.2 Supplementary Figure 2j: Heatmap of Marker Genes
cat("Generating Heatmap...\n")
marker_features <- c(
  "Tmem212", "Mrc1", "Cd3e", "Aif1", "Tagln", "Slc6a13", "Vtn", "Cldn5", 
  "Snap25", "Pax6", "Plp1", "Bmp4", "Pdgfra", "Slc1a3"
)

# Calculate Average Expression
avg_list <- AverageExpression(ref.merge, features = marker_features, assays = "SCT", slot = "data")
avg_mat <- avg_list$SCT

# Reorder columns (clusters) and rows (genes)
valid_cols <- plot_order[plot_order %in% colnames(avg_mat)]
avg_mat <- avg_mat[marker_features, valid_cols]

# Draw Heatmap
pheatmap(
  avg_mat,
  scale = "row",
  border_color = "white",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
)

### 3.3 Supplementary Figure 2k: Violin Plot (Bmp4)
cat("Generating Violin Plot...\n")
p_vln <- VlnPlot(ref.merge, features = "Bmp4", assay = "SCT", pt.size = 0.01)

# Reorder factor levels for x-axis
p_vln$data$ident <- factor(p_vln$data$ident, levels = plot_order)

p_vln_final <- p_vln +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(p_vln_final)
cat("Finding All Markers (Method: MAST)... This may take time.\n")

# Prepare object for DE analysis
ref.merge <- PrepSCTFindMarkers(ref.merge, assay = "SCT")

# Find Markers
ref.merge.markers <- FindAllMarkers(
  ref.merge, 
  logfc.threshold = 0.1,
  only.pos = TRUE,
  assay = "SCT", 
  test.use = "MAST", 
  recorrect_umi = FALSE
)

# Export to Excel
cat("Saving markers to Excel...\n")
df_list <- split(ref.merge.markers, ref.merge.markers$cluster)
write_xlsx(df_list, path = MARKER_OUT_PATH)

