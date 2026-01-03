################################################################################
# Helper Functions for Mouse Single Cell Analysis
# Target Species: Mus musculus
# Updated: 2023-06-01
################################################################################

# ==============================================================================
# 1. Environment & Libraries
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(moments) # For skewness in scStress
  library(org.Mm.eg.db) # Mouse Annotation
})

# ==============================================================================
# 2. QC & Preprocessing (Mouse Optimized)
# ==============================================================================

#' Preprocessing Wrapper (Process A) - Mouse Optimized
#' Includes dynamic thresholding for outliers
#' @param run Seurat Object
#' @return Filtered Seurat Object
processA <- function(run) {
  # 1. Calculate QC Metrics using MOUSE gene patterns
  # Note: Mouse mito genes usually start with lower case 'mt-' or 'mt.'
  run[["percent.mt"]] <- PercentageFeatureSet(run, pattern = "^mt\\.|^mt-")
  
  # Mouse Ribosomal (Rpl/Rps)
  run[["percent.ribosomal"]] <- PercentageFeatureSet(run, pattern = "^Rpl|^Rps")
  
  # Mouse Hemoglobin (Hba/Hbb)
  run[["percent.Hb"]] <- PercentageFeatureSet(run, pattern = "^Hb[ab]")
  
  # Complexity Metric
  run$log10GenesPerUMI <- log10(run$nFeature_RNA) / log10(run$nCount_RNA)
  
  # 2. Dynamic Threshold Calculation (Median + MAD)
  # Using robust statistics to identify outliers
  stats <- list(
    nCount = c(median(run$nCount_RNA), mad(run$nCount_RNA)),
    nFeature = c(median(run$nFeature_RNA), mad(run$nFeature_RNA)),
    pctMT = c(median(run$percent.mt), mad(run$percent.mt))
  )
  
  thresh_nCount <- stats$nCount[1] + 5 * stats$nCount[2]
  thresh_nFeature <- stats$nFeature[1] + 5 * stats$nFeature[2]
  thresh_mt <- stats$pctMT[1] + 2 * stats$pctMT[2]
  
  # Ensure MT threshold isn't too strict (minimum 5%)
  thresh_mt <- max(thresh_mt, 5)
  
  # 3. Filtering
  run <- subset(run, subset = (nFeature_RNA >= 500) & 
                  (percent.mt <= thresh_mt) & 
                  (log10GenesPerUMI > 0.80) & 
                  (nCount_RNA <= thresh_nCount) & 
                  (nFeature_RNA <= thresh_nFeature))
  
  return(run)
}

# ==============================================================================
# 3. Data Loading Helpers
# ==============================================================================

#' Custom Read10X for GEO (Handles .gz files)
#' @param data.dir Directory containing the sample files
#' @param sample.names File prefix (e.g. "GSM123_")
Read10X_GEO <- function(data.dir, sample.names = NULL, gene.column = 2) {
  if (!dir.exists(data.dir)) stop("Directory not found: ", data.dir)
  
  # Define file paths
  files <- list(
    barcode = file.path(data.dir, paste0(sample.names, 'barcodes.tsv.gz')),
    matrix = file.path(data.dir, paste0(sample.names, 'matrix.mtx.gz')),
    features = file.path(data.dir, paste0(sample.names, 'features.tsv.gz')),
    genes = file.path(data.dir, paste0(sample.names, 'genes.tsv.gz')) # Legacy support
  )
  
  if (!file.exists(files$barcode)) stop("Barcode file missing")
  if (!file.exists(files$matrix)) stop("Matrix file missing")
  
  # Determine feature file to use
  feat_file <- if (file.exists(files$features)) files$features else files$genes
  if (!file.exists(feat_file)) stop("Features/Genes file missing")
  
  # Load Data
  data <- Matrix::readMM(file = files$matrix)
  cell.names <- readLines(files$barcode)
  colnames(data) <- cell.names
  
  feature.names <- read.delim(feat_file, header = FALSE, stringsAsFactors = FALSE)
  
  # Handle Row names (Gene Symbols)
  # Force unique names to avoid duplicates error
  gene_ids <- if (ncol(feature.names) >= gene.column) feature.names[, gene.column] else feature.names[, 1]
  rownames(data) <- make.unique(gene_ids)
  
  return(list(data))
}

# ==============================================================================
# 4. Stress Analysis (Mouse)
# ==============================================================================

#' Calculate Stress Score and Filter Cells
#' @param obj Seurat Object
#' @param stress_genes Character vector of mouse stress genes
#' @param cut Quantile cutoff (default 0.95)
scStress <- function(obj, stress_genes, cut = 0.95){
  
  # 1. Intersect genes
  valid_genes <- intersect(rownames(obj), stress_genes)
  if (length(valid_genes) < 5) {
    warning("Too few stress genes found in data. Skipping stress filtering.")
    return(obj)
  }
  
  # 2. Extract Data for PCA
  # Using counts slot is standard for this specific stress method (as per original code)
  stress_mat <- GetAssayData(obj, slot = "counts")[valid_genes, ]
  stress_mat <- as.matrix(stress_mat)
  stress_mat <- t(stress_mat[apply(stress_mat, 1, var) != 0, ]) # Remove zero var genes
  
  # 3. PCA
  pca <- prcomp(stress_mat, center = TRUE, scale. = TRUE)
  
  # 4. Handle Skewness
  # Stress score should be positively correlated with stress gene expression
  pc1 <- pca$x[, 1]
  if (moments::skewness(pc1) < 0) {
    pc1 <- -pc1
  }
  
  # 5. Define Cutoff and Filter
  cutoff_val <- quantile(pc1, probs = cut)
  stress_cells <- names(pc1[pc1 > cutoff_val])
  
  # Optional: Store score in metadata before filtering (for checking)
  obj <- AddMetaData(obj, metadata = pc1, col.name = "Stress_Score_PC1")
  
  message(paste0("Filtering ", length(stress_cells), " stress cells (>", cut*100, "th percentile)."))
  
  # Return subset
  return(subset(obj, cells = setdiff(colnames(obj), stress_cells)))
}

# ==============================================================================
# 5. Utilities
# ==============================================================================

#' Convert Human Gene List to Mouse Symbols
#' Requires biomaRt
convert_human_to_mouse <- function(human_genes){
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 <- getLDS(
    attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
    values = human_genes, mart = human, 
    attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = TRUE
  )
  
  return(unique(genesV2[, 2]))
}

#' Identify Feature Types to Exclude
#' Helper for finding MT, HB, Ribosomal genes
out_features <- function(x){
  genes <- rownames(x)
  out <- c(
    grep("^mt\\.|^mt-", genes, value = TRUE),
    grep("^Hb[ab]", genes, value = TRUE),
    grep("^Rpl|^Rps", genes, value = TRUE)
  )
  return(unique(out))
}
