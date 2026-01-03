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

# ==============================================================================
# 3. Cell Type Annotation Update
# ==============================================================================

# Define Cluster to Cell Type Mapping
# Using a named vector for efficient mapping instead of repeated subset calls
cluster_mapping <- c(
  "4" = "OPC", "6" = "OPC", "12" = "OPC", "14" = "OPC", "18" = "OPC", "19" = "OPC", "23" = "OPC",
  "11" = "COP", "21" = "COP",
  "0" = "Oligodendrocytes", "1" = "Oligodendrocytes", "2" = "Oligodendrocytes", "3" = "Oligodendrocytes",
  "5" = "Oligodendrocytes", "7" = "Oligodendrocytes", "8" = "Oligodendrocytes", "9" = "Oligodendrocytes",
  "10" = "Oligodendrocytes", "13" = "Oligodendrocytes", "15" = "Oligodendrocytes", "16" = "Oligodendrocytes",
  "17" = "Oligodendrocytes", "20" = "Oligodendrocytes"
)

# Apply mapping
# Note: Ensure 'idents' corresponds to the clustering resolution used for mapping (e.g., SCT_snn_res.0.4)
# Assuming OL.mergedd active ident is the cluster ID
current_clusters <- as.character(Idents(OL.mergedd))
new_labels <- cluster_mapping[current_clusters]
# Assign original cluster ID if no mapping found (optional safety check)
# new_labels[is.na(new_labels)] <- current_clusters[is.na(new_labels)]
OL.mergedd$new_cell_type <- new_labels

# ==============================================================================
# 4. Genotype Annotation Update
# ==============================================================================

# Subset for OL.COP specific analysis
OL.COP <- subset(OL.merged, idents = c(11, 21))

# Define Genotype Mapping based on ID
# Assuming 'id' corresponds to specific experimental conditions or samples
# Refactored from repetitive 'which' statements to a loop or vectorized approach
genotype_map <- c(
  "1" = "PS2/APP/P301L", "2" = "PS2/APP/P301L/TREM2KO", "3" = "PS2/APP/P301L",
  "4" = "PS2/APP/P301L/TREM2KO", "5" = "PS2/APP/P301L", "6" = "P301L",
  "7" = "WT", "8" = "WT", "9" = "WT", "10" = "PS2/APP/P301L/TREM2KO",
  "11" = "P301L", "12" = "P301L", "13" = "WT", "14" = "PS2/APP",
  "15" = "WT", "16" = "PS2/APP", "17" = "WT", "18" = "PS2/APP",
  "19" = "WT", "20" = "WT", "21" = "WT", "22" = "PS2/APP",
  "23" = "WT", "24" = "PS2/APP", 
  "25" = "WT/Ctrl_0wks", "26" = "WT/Ctrl_3wks", 
  "27" = "WT/CupRap_3wks", "28" = "WT/CupRap_3wks",
  "29" = "NesCre/Ctrl_3wks", "30" = "NesCre/Ctrl_3wks",
  "31" = "NesCre/CupRap_3wks", "32" = "NesCre/CupRap_3wks",
  "33" = "WT/CupRap_0wk", "34" = "WT/CupRap_0wk",
  "35" = "WT/CTRL", "36" = "WT/CPZ3", "37" = "WT/XPro1595_CPZ3"
)

# Apply Genotype Mapping
# Note: This logic assumes 'id' in metadata matches the index order of table(OL.COP$id)
# Ideally, map directly from ID string if possible. Assuming 'id' is a factor or character matching the map keys.
# Since the original code used index-based mapping from table(), we replicate that logic safely.
unique_ids <- names(table(OL.COP$id))
for(i in seq_along(unique_ids)){
  target_id <- unique_ids[i]
  if(as.character(i) %in% names(genotype_map)){
     OL.COP$genetype[OL.COP$id == target_id] <- genotype_map[as.character(i)]
  }
}

# Handle Special Case for GSE278199
OL.COP$genetype[OL.COP$basic == "GSE278199" & is.na(OL.COP$genetype)] <- "WT/CPZ"


# ==============================================================================
# 5. UMAP Plots Generation
# ==============================================================================

# Define specific sample IDs for subset plotting (e.g., PS2APP cohort)
target_subset_ids <- names(table(OL.mergedd$id))[c(7:9, 13:24)]

# Generate Plots : Supplementary Figure 2a
sp2a <- plot_custom_dimplot(OL.mergedd, group_by = "new_cell_type", label = FALSE)

# ==============================================================================
# 6. Differential Expression (COP Markers)
# ==============================================================================

# Subset for PS2APP analysis
PS2APP <- subset(OL.mergedd, subset = id %in% target_subset_ids)
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
# 7. Comparison Analysis (Heatmap across Datasets)
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
