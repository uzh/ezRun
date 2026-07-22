#!/usr/bin/env Rscript
# Generate test HTML report for ScSeuratCompare with updated PCA

library(Seurat)
library(qs2)
library(writexl)

# Create output directory
output_dir <- "/home/pgueguen/git/ezRun/test_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

cat("Creating test dataset...\n")

# Create a more realistic test Seurat object
set.seed(42)
n_genes <- 2000
n_cells <- 3000

# Create count matrix with some structure
counts <- matrix(
  rpois(n_genes * n_cells, lambda = 3),
  nrow = n_genes,
  ncol = n_cells
)
rownames(counts) <- paste0("Gene_", 1:n_genes)
colnames(counts) <- paste0("Cell_", 1:n_cells)

# Add some differential expression between groups
control_cells <- 1:1500
treatment_cells <- 1501:3000
counts[1:100, treatment_cells] <- counts[1:100, treatment_cells] * 3

# Create Seurat object
scData <- CreateSeuratObject(counts = counts, project = "TestCompare")

# Add metadata
scData$Sample <- rep(c("Sample_A", "Sample_B", "Sample_C", "Sample_D",
                       "Sample_E", "Sample_F"), each = 500)
scData$condition <- rep(c("Control", "Treatment"), each = 1500)
scData$celltype <- sample(c("CellType_1", "CellType_2", "CellType_3"),
                          n_cells, replace = TRUE)
scData$seurat_clusters <- as.character(sample(0:4, n_cells, replace = TRUE))

# Process the object
scData <- NormalizeData(scData, verbose = FALSE)
scData <- FindVariableFeatures(scData, verbose = FALSE)
scData <- ScaleData(scData, verbose = FALSE)
scData <- RunPCA(scData, npcs = 20, verbose = FALSE)
scData <- RunUMAP(scData, dims = 1:10, verbose = FALSE)

# SCTransform
scData <- SCTransform(scData, verbose = FALSE)

cat("Saving test Seurat object...\n")
qs_save(scData, "scData.qs2", nthreads = 4)

# Create param object
param <- list(
  name = "Test ScSeuratCompare Report",
  grouping = "condition",
  sampleGroup = "Treatment",
  refGroup = "Control",
  CellIdentity = "celltype",
  replicateGrouping = "Sample",
  DE.method = "wilcox",
  DE.regress = "Batch",
  pseudoBulkMode = FALSE,
  sccomp.variability = FALSE,
  refBuild = "Homo_sapiens/Ensembl/GRCh38/Annotation/Release_110-2023-10-30",
  cores = 4
)

cat("Saving param.rds...\n")
saveRDS(param, "param.rds")

# Get actual genes from the Seurat object (genes have dashes, not underscores)
actual_genes <- rownames(scData)

# Create mock conserved markers using real gene names
consMarkers <- data.frame(
  cluster = rep(c("CellType_1", "CellType_2", "CellType_3"), each = 10),
  gene = actual_genes[1:30],
  avg_avg_fc = runif(30, 1.5, 4),
  p_val = runif(30, 0, 0.01),
  p_val_adj = runif(30, 0, 0.05)
)

cat("Saving consMarkers.xlsx...\n")
write_xlsx(consMarkers, "consMarkers.xlsx")

# Create mock differential genes using real gene names
diffGenes <- data.frame(
  cluster = rep(c("CellType_1", "CellType_2", "CellType_3"), each = 20),
  gene = actual_genes[31:90],
  avg_log2FC = c(runif(30, 0.5, 2), runif(30, -2, -0.5)),
  p_val = runif(60, 0, 0.01),
  p_val_adj = runif(60, 0, 0.05),
  pct.1 = runif(60, 0.3, 0.8),
  pct.2 = runif(60, 0.2, 0.7),
  diff_pct = runif(60, 0.1, 0.3)
)

cat("Saving diffGenes.xlsx...\n")
write_xlsx(diffGenes, "diffGenes.xlsx")

# Copy the updated Rmd template
cat("Copying Rmd template...\n")
file.copy(
  "/home/pgueguen/git/ezRun/inst/templates/ScSeuratCompare.Rmd",
  "ScSeuratCompare.Rmd",
  overwrite = TRUE
)

# Render the report
cat("\nRendering HTML report...\n")
rmarkdown::render(
  "ScSeuratCompare.Rmd",
  output_file = "ScSeuratCompare_Test_Report.html",
  envir = new.env()
)

cat("\n✓ Report generated successfully!\n")
cat("\nReport location:\n")
cat(file.path(output_dir, "ScSeuratCompare_Test_Report.html"), "\n")
