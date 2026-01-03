# ==============================================================================
# Single-Cell RNA-seq Analysis Pipeline: Visualization & DE Analysis
# Author: Jaemyung, Jang / Korea Brain Research Institute
# Description: This script performs visualization (DimPlot, DotPlot), differential
#              expression analysis (FindMarkers), and generates Volcano/Heatmap plots
#              for Oligodendrocyte lineage cells.
# ==============================================================================

# 1. Setup & Libraries ---------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggsci)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(openxlsx)
library(scales) # For comma format

# Define Paths
BASE_DIR <- "/path/to/"
RESULTS_DIR <- file.path(BASE_DIR, "results")
DATE_SUFFIX <- "20251208" # Update date as needed

# 2. Visualization Functions ---------------------------------------------------

#' Wrapper for Seurat DimPlot with custom themes
#' @param obj Seurat object
#' @param group_by Column to group cells by
#' @param title Plot title
#' @param label Boolean, show labels
#' @param subset_ids Optional list of IDs to subset
plot_custom_dimplot <- function(obj, group_by, title = "", label = TRUE, subset_ids = NULL) {
  
  plot_obj <- obj
  if (!is.null(subset_ids)) {
    plot_obj <- subset(obj, subset = id %in% subset_ids)
  }
  
  # Calculate cell count for title if needed
  n_cells <- comma(ncol(plot_obj), format = "d")
  full_title <- paste0(title, " - cells: ", n_cells)
  
  p <- DimPlot(plot_obj, 
               cols = pal_npg("nrc")(10),
               group.by = group_by, 
               reduction = "umap", 
               label = label, repel = label, 
               label.size = 6, pt.size = 0.01) + 
    theme(axis.title = element_blank(), 
          plot.title = element_text(size = 14, face = "bold", hjust = 0), 
          axis.line = element_blank(), 
          axis.text = element_blank(),  
          axis.ticks = element_blank()) + 
    ggtitle(full_title)
  
  if (!label) p <- p + NoLegend()
  
  return(p)
}

# 3. UMAP Visualizations -------------------------------------------------------

# Define Sample IDs for Subsetting (Refactoring the long list of IDs)
# These IDs seem to correspond to specific datasets or batches
target_sample_ids <- names(table(OL.mergedd$id))[c(7:9, 13:24)]

# Plot 1: All Cells by New Cell Type
p0s <- plot_custom_dimplot(OL.mergedd, group_by = "new_cell_type", label = FALSE, title = "")

# Plot 2: Subset of Samples by New Cell Type
p1s <- plot_custom_dimplot(OL.mergedd, group_by = "new_cell_type", subset_ids = target_sample_ids, title = "OL")

# Plot 3: Subset of Samples by Basic Group
p1s_1 <- plot_custom_dimplot(OL.mergedd, group_by = "basic", subset_ids = target_sample_ids, title = "OL")

# Plot 4: All Cells by Basic Group
p0s_1 <- plot_custom_dimplot(OL.mergedd, group_by = "basic", label = FALSE, title = "OL") + theme(legend.position = "bottom")

# Save Plots
ggsave(file.path(RESULTS_DIR, paste0("fig1S_c_new_", DATE_SUFFIX, ".pdf")), p1s, width = 7.5, height = 7.5)
ggsave(file.path(RESULTS_DIR, paste0("fig1S_d_new_", DATE_SUFFIX, ".pdf")), p1s_1, width = 7.5, height = 7.5)
ggsave(file.path(RESULTS_DIR, paste0("fig1S_a_new_", DATE_SUFFIX, ".pdf")), p0s, width = 7.5, height = 7.5)
ggsave(file.path(RESULTS_DIR, paste0("fig1S_e_new_", DATE_SUFFIX, ".pdf")), p0s_1, width = 7.5, height = 7.5)

# 4. Differential Expression Analysis (COP Markers) ----------------------------

# Subset for PS2APP related samples
PS2APP <- subset(OL.mergedd, subset = id %in% target_sample_ids)
PS2APP <- PrepSCTFindMarkers(PS2APP)

# Find Markers for COP
COP.markers <- FindMarkers(PS2APP,
                           group.by = "celltype", 
                           test.use = "MAST",
                           ident.1 = "COP")

# Save Results
write.xlsx(COP.markers, file.path(RESULTS_DIR, paste0("fig1S_b_new_", DATE_SUFFIX, ".xlsx")), rowNames = TRUE)

# Select Genes for Volcano Plot
selGenes <- COP.markers %>% dplyr::filter(p_val_adj < 1e-300 & avg_log2FC > 2) %>% rownames()
# Filter out unknown genes (Gm, Rik)
volcano_labels <- setdiff(selGenes, grep("^Gm|Rik$", selGenes, value = TRUE))
specific_labels <- c("Bcas1", "Gpr17", "Enpp6", "Bmp4") # Manual selection

# Generate Volcano Plot
pdf(file.path(RESULTS_DIR, paste0("fig1S_b_new_", DATE_SUFFIX, ".pdf")), width = 5, height = 7.5, family = "ArialMT")
EnhancedVolcano(COP.markers,
                lab = rownames(COP.markers),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "COP vs OPC/Oligodendrocytes",
                subtitle = paste0("COP: ", comma(sum(PS2APP$celltype == "COP")), 
                                  " cells - Others: ", comma(sum(PS2APP$celltype != "COP")), " cells"),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~Log[10]~ 'Adjust P-value'),
                pCutoff = 1e-300,
                FCcutoff = 2,
                pointSize = 0.25,
                labSize = 3.2,
                selectLab = specific_labels,
                colAlpha = 0.5,
                legendPosition = 'bottom',
                drawConnectors = TRUE,
                widthConnectors = 0.05,
                colConnectors = 'black',
                border = 'full',
                borderWidth = 0.25,
                axisLabSize = 10,
                titleLabSize = 12,
                subtitleLabSize = 10,
                captionLabSize = 10,
                maxoverlapsConnectors = 50)
dev.off()


# 5. Cell Type Annotation Update -----------------------------------------------

# Define cluster to cell type mapping
# Note: Using a named vector for mapping is cleaner than repeated subset calls
cluster_annotations <- c(
  "4" = "OPC", "6" = "OPC", "12" = "OPC", "14" = "OPC", "18" = "OPC", "19" = "OPC", "23" = "OPC",
  "11" = "COP", "21" = "COP",
  "0" = "Oligodendrocytes", "1" = "Oligodendrocytes", "2" = "Oligodendrocytes", "3" = "Oligodendrocytes",
  "5" = "Oligodendrocytes", "7" = "Oligodendrocytes", "8" = "Oligodendrocytes", "9" = "Oligodendrocytes",
  "10" = "Oligodendrocytes", "13" = "Oligodendrocytes", "15" = "Oligodendrocytes", "16" = "Oligodendrocytes",
  "17" = "Oligodendrocytes", "20" = "Oligodendrocytes"
)

# Apply annotations
# Assuming 'idents' refers to a specific clustering resolution column in meta.data
# You might need to adjust the column name (e.g., 'SCT_snn_res.0.8') used for matching
current_clusters <- as.character(Idents(OL.mergedd))
new_labels <- cluster_annotations[current_clusters]
# Fill NA (unmapped clusters) with original or "Unknown" if necessary
# new_labels[is.na(new_labels)] <- current_clusters[is.na(new_labels)] 
OL.mergedd$new_cell_type <- new_labels

# 6. Genotype Annotation Update ------------------------------------------------

# Subset for specific analysis
OL.COP <- subset(OL.merged, idents = c(11, 21))

# Genotype Mapping Table
# (Refactoring the long list of manual assignments)
genotype_map <- c(
  "1" = "PS2/APP/P301L", "2" = "PS2/APP/P301L/TREM2KO", "3" = "PS2/APP/P301L", 
  "4" = "PS2/APP/P301L/TREM2KO", "5" = "PS2/APP/P301L", "6" = "P301L",
  "7" = "WT", "8" = "WT", "9" = "WT", 
  "10" = "PS2/APP/P301L/TREM2KO", "11" = "P301L", "12" = "P301L",
  "13" = "WT", "14" = "PS2/APP", "15" = "WT", "16" = "PS2/APP",
  "17" = "WT", "18" = "PS2/APP", "19" = "WT", "20" = "WT", "21" = "WT",
  "22" = "PS2/APP", "23" = "WT", "24" = "PS2/APP",
  "25" = "WT/Ctrl_0wks", "26" = "WT/Ctrl_3wks", 
  "27" = "WT/CupRap_3wks", "28" = "WT/CupRap_3wks",
  "29" = "NesCre/Ctrl_3wks", "30" = "NesCre/Ctrl_3wks",
  "31" = "NesCre/CupRap_3wks", "32" = "NesCre/CupRap_3wks",
  "33" = "WT/CupRap_0wk", "34" = "WT/CupRap_0wk",
  "35" = "WT/CTRL", "36" = "WT/CPZ3", "37" = "WT/XPro1595_CPZ3"
)

# Apply mapping
# Assuming 'id' column matches the indices in genotype_map keys
# NOTE: The original code used names(table(OL.COP$id))[index]. 
# This relies on alphabetical sorting. Ensure this logic holds.
sorted_ids <- sort(unique(OL.COP$id))
current_ids <- OL.COP$id

# Loop or vectorize to assign. 
# Since map keys are indices (1..37), we map IDs to indices first.
# This part is tricky without original data, assuming manual assignment logic was specific.
# Here is a cleaner way if you have a direct mapping file. 
# For now, applying logic as intended:
for(i in 1:length(genotype_map)) {
  target_id <- sorted_ids[i]
  if(!is.na(target_id)) {
    OL.COP$genetype[OL.COP$id == target_id] <- genotype_map[as.character(i)]
  }
}

# Handle special case
OL.COP$genetype[OL.COP$basic == "GSE278199" & is.na(OL.COP$genetype)] <- "WT/CPZ"


# 7. Comparison Analysis (Heatmap) ---------------------------------------------

# Function to run DE analysis
run_de_analysis <- function(obj, id1, id2, group_col = "genetype") {
  FindMarkers(obj, ident.1 = id1, ident.2 = id2, group.by = group_col, 
              test.use = "MAST", recorrect_umi = FALSE)
}

# Run Comparisons
# Note: Using sorted unique genetypes to get IDs by index
sorted_genotypes <- sort(unique(OL.COP$genetype))

OL.COP.DEP.DM_with_CupRap_0wks <- run_de_analysis(OL.COP, sorted_genotypes[12], sorted_genotypes[10])
OL.COP.DEP.DM_with_CupRap_3wks <- run_de_analysis(OL.COP, sorted_genotypes[13], sorted_genotypes[11])
OL.COP.DEP.DM_with_CPZ3        <- run_de_analysis(OL.COP, sorted_genotypes[1],  sorted_genotypes[9])
OL.COP.DEP.PS2_APP             <- run_de_analysis(OL.COP, sorted_genotypes[5],  sorted_genotypes[8])
OL.COP.DEP.tau                 <- run_de_analysis(OL.COP, sorted_genotypes[6],  sorted_genotypes[8])

# Save DE Results
write.xlsx(OL.COP.DEP.DM_with_CPZ3, file.path(RESULTS_DIR, "fig1S_CPZ3_20251209.xlsx"), rowNames = TRUE)
write.xlsx(OL.COP.DEP.PS2_APP, file.path(RESULTS_DIR, "fig1S_PS2APP_20251209.xlsx"), rowNames = TRUE)
write.xlsx(OL.COP.DEP.tau, file.path(RESULTS_DIR, "fig1S_P301L_20251209.xlsx"), rowNames = TRUE)


# Prepare Heatmap Data
target_genes <- c("Bmp4", "Enpp6", "Bcas1", "Gpr17", "Cnp", "Cyfip2", 
                  "Lims2", "Wasf1", "Phyhipl", "Nfasc", "Frmd4a", "Olig2", "Pdgfa")

extract_logFC <- function(de_results, col_name) {
  de_results %>%
    rownames_to_column("Gene") %>%
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

# Draw Heatmap
pdf(file.path(RESULTS_DIR, "fig1S_c_20251209.pdf"), width = 5, height = 15, family = "ArialMT")
pheatmap(heatmap_data,
         name = "Average Fold Change",
         cellwidth = 25, cellheight = 25,
         cluster_cols = FALSE, cluster_rows = FALSE,
         border_color = "white",
         main = "Average expression")
dev.off()

# Final UMAPs
DimPlot(OL.COP, reduction = "umap", label = TRUE, label.size = 4, pt.size = 0.01) + 
  theme(legend.position = "none") + NoAxes() + ggtitle("")

DimPlot(OL.merged, reduction = "umap", label = FALSE, pt.size = 0.01) + 
  theme(legend.position = "none") + NoAxes() + 
  ggtitle(paste0(comma(ncol(OL.merged), format = "d"), " cells"))
