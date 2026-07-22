#!/usr/bin/env Rscript
# Test script for updated pseudobulk PCA function

library(Seurat)
library(Matrix)

# Create a minimal test Seurat object
set.seed(123)
n_genes <- 2000
n_cells <- 1000

# Create random count matrix
counts <- matrix(
  rpois(n_genes * n_cells, lambda = 5),
  nrow = n_genes,
  ncol = n_cells
)
rownames(counts) <- paste0("Gene", 1:n_genes)
colnames(counts) <- paste0("Cell", 1:n_cells)

# Create Seurat object
pbmc3k <- CreateSeuratObject(counts = counts, project = "test")

# Add fake donor information (simulating biological replicates)
pbmc3k$donor <- sample(x = c('A', 'B', 'C', 'D'), size = ncol(pbmc3k), replace = TRUE)

# Add fake condition information
pbmc3k$condition <- sample(x = c('Control', 'Treatment'),
                           size = ncol(pbmc3k), replace = TRUE)

# Add fake cluster information
pbmc3k$seurat_clusters <- sample(x = c('0', '1', '2', '3'),
                                  size = ncol(pbmc3k), replace = TRUE)
Idents(pbmc3k) <- pbmc3k$seurat_clusters

cat("Testing new run_pseudobulk_pca() function...\n\n")

# Test function definition
run_pseudobulk_pca <- function(obj, group_col) {
  # Use AggregateExpression to create pseudobulk Seurat object
  pb_seurat <- AggregateExpression(
    obj,
    assays = "RNA",
    return.seurat = TRUE,
    group.by = group_col
  )

  # Apply Seurat normalization workflow
  pb_seurat <- NormalizeData(pb_seurat, verbose = FALSE)
  pb_seurat <- FindVariableFeatures(pb_seurat, verbose = FALSE)
  pb_seurat <- ScaleData(pb_seurat, verbose = FALSE)
  pb_seurat <- RunPCA(pb_seurat, npcs = min(10, ncol(pb_seurat) - 1), verbose = FALSE)

  # Extract PCA embeddings
  pca_embeddings <- Embeddings(pb_seurat, reduction = "pca")

  # Extract variance explained
  stdev <- Stdev(pb_seurat, reduction = "pca")
  var_explained <- stdev^2 / sum(stdev^2)
  var_pct <- round(var_explained * 100, 1)

  # Determine number of valid dimensions
  n_dims <- min(2, ncol(pca_embeddings))

  # Create dataframe for plotting
  pca_df <- data.frame(
    PC1 = pca_embeddings[, 1],
    PC2 = if (n_dims >= 2) pca_embeddings[, 2] else 0,
    SampleID = rownames(pca_embeddings)
  )

  list(pca_df = pca_df, var_pct = var_pct, n_dims = n_dims)
}

# Test 1: Run PCA on all cells by donor
cat("Test 1: Running pseudobulk PCA by donor...\n")
tryCatch({
  pca_res <- run_pseudobulk_pca(pbmc3k, "donor")
  cat("✓ PCA completed successfully\n")
  cat("  - Number of samples:", nrow(pca_res$pca_df), "\n")
  cat("  - PC1 variance explained:", pca_res$var_pct[1], "%\n")
  cat("  - PC2 variance explained:", pca_res$var_pct[2], "%\n")
  cat("  - Sample IDs:", paste(pca_res$pca_df$SampleID, collapse = ", "), "\n\n")
}, error = function(e) {
  cat("✗ Test 1 failed:", e$message, "\n\n")
})

# Test 2: Run PCA on subset (simulate per-cluster analysis)
cat("Test 2: Running pseudobulk PCA on subset (cluster 0)...\n")
tryCatch({
  Idents(pbmc3k) <- pbmc3k$seurat_clusters
  cells_cl0 <- WhichCells(pbmc3k, idents = "0")
  pbmc3k_cl0 <- subset(pbmc3k, cells = cells_cl0)

  pca_res_cl0 <- run_pseudobulk_pca(pbmc3k_cl0, "donor")
  cat("✓ Subset PCA completed successfully\n")
  cat("  - Number of samples:", nrow(pca_res_cl0$pca_df), "\n")
  cat("  - PC1 variance explained:", pca_res_cl0$var_pct[1], "%\n\n")
}, error = function(e) {
  cat("✗ Test 2 failed:", e$message, "\n\n")
})

# Test 3: Verify normalization is applied (check that values are log-normalized)
cat("Test 3: Verifying normalization workflow...\n")
tryCatch({
  pb_seurat <- AggregateExpression(
    pbmc3k,
    assays = "RNA",
    return.seurat = TRUE,
    group.by = "donor"
  )

  cat("✓ Raw pseudobulk created\n")

  # Check raw counts
  raw_counts_sample <- GetAssayData(pb_seurat, assay = "RNA", layer = "counts")[1:5, 1]
  cat("  - Sample raw counts (first 5 genes):", paste(round(raw_counts_sample, 2), collapse = ", "), "\n")

  # Apply normalization
  pb_seurat <- NormalizeData(pb_seurat, verbose = FALSE)
  norm_data_sample <- GetAssayData(pb_seurat, assay = "RNA", layer = "data")[1:5, 1]
  cat("  - Sample normalized data (first 5 genes):", paste(round(norm_data_sample, 2), collapse = ", "), "\n")

  # Verify normalization happened (log-normalized values should be different from raw counts)
  if (all(norm_data_sample != raw_counts_sample)) {
    cat("✓ Normalization confirmed: values differ from raw counts\n\n")
  } else {
    cat("⚠ Warning: Normalized values identical to raw counts\n\n")
  }
}, error = function(e) {
  cat("✗ Test 3 failed:", e$message, "\n\n")
})

cat("All tests completed!\n")
