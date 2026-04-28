context("multiOmics utils")

test_that("detectModalities flags only RNA for plain CellRanger output", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "filtered_feature_bc_matrix.h5"))
  m <- detectModalities(d)
  expect_true(m$hasRNA)
  expect_false(m$hasADT)
  expect_false(m$hasVDJ_T)
  expect_false(m$hasVDJ_B)
  expect_false(m$hasATAC)
})

test_that("detectModalities flags VDJ-T and VDJ-B siblings (CellRanger Multi)", {
  root <- tempfile(); dir.create(file.path(root, "count"), recursive = TRUE)
  file.create(file.path(root, "count", "sample_filtered_feature_bc_matrix.h5"))
  dir.create(file.path(root, "vdj_t")); file.create(file.path(root, "vdj_t", "filtered_contig_annotations.csv"))
  dir.create(file.path(root, "vdj_b")); file.create(file.path(root, "vdj_b", "filtered_contig_annotations.csv"))
  m <- detectModalities(file.path(root, "count"))
  expect_true(m$hasVDJ_T)
  expect_true(m$hasVDJ_B)
})

test_that("detectModalities flags ATAC siblings (CellRanger ARC)", {
  d <- tempfile(); dir.create(d)
  file.create(file.path(d, "filtered_feature_bc_matrix.h5"))
  file.create(file.path(d, "atac_fragments.tsv.gz"))
  file.create(file.path(d, "atac_peaks.bed"))
  m <- detectModalities(d)
  expect_true(m$hasATAC)
})

test_that("h5HasAntibodyCapture returns TRUE when feature_type contains Antibody Capture", {
  skip_if_not_installed("hdf5r")
  h5 <- tempfile(fileext = ".h5")
  on.exit(unlink(h5))
  hf <- hdf5r::H5File$new(h5, mode = "w")
  g <- hf$create_group("matrix")
  fg <- g$create_group("features")
  fg[["feature_type"]] <- c("Gene Expression", "Antibody Capture")
  hf$close_all()
  expect_true(h5HasAntibodyCapture(h5))
})

test_that("h5HasAntibodyCapture returns FALSE when only Gene Expression features", {
  skip_if_not_installed("hdf5r")
  h5 <- tempfile(fileext = ".h5")
  on.exit(unlink(h5))
  hf <- hdf5r::H5File$new(h5, mode = "w")
  g <- hf$create_group("matrix")
  fg <- g$create_group("features")
  fg[["feature_type"]] <- c("Gene Expression", "Gene Expression")
  hf$close_all()
  expect_false(h5HasAntibodyCapture(h5))
})

test_that("detectModalities flags ADT when H5 has Antibody Capture features", {
  skip_if_not_installed("hdf5r")
  d <- tempfile(); dir.create(d)
  h5 <- file.path(d, "filtered_feature_bc_matrix.h5")
  hf <- hdf5r::H5File$new(h5, mode = "w")
  g <- hf$create_group("matrix")
  fg <- g$create_group("features")
  fg[["feature_type"]] <- c("Gene Expression", "Antibody Capture")
  hf$close_all()
  m <- detectModalities(d)
  expect_true(m$hasRNA)
  expect_true(m$hasADT)
})

test_that("detectModalities handles FGCZ CellRanger Multi layout (mtx dir as CountMatrix)", {
  # FGCZ writes CountMatrix as the mtx directory; H5 lives one level up; VDJ two levels up.
  root <- tempfile()
  cm  <- file.path(root, "per_sample_outs", "S-cellRanger", "count", "sample_filtered_feature_bc_matrix")
  dir.create(cm, recursive = TRUE)
  file.create(file.path(dirname(cm), "sample_filtered_feature_bc_matrix.h5"))
  vdj_t <- file.path(dirname(dirname(cm)), "vdj_t"); dir.create(vdj_t)
  file.create(file.path(vdj_t, "filtered_contig_annotations.csv"))
  m <- detectModalities(cm)
  expect_true(m$hasRNA)
  expect_true(m$hasVDJ_T)
})

test_that("readADTCounts returns only Antibody Capture features from a multi-modal H5", {
  skip_on_cran(); skip_if_not_installed("Seurat"); skip_if_not_installed("hdf5r")
  # Build a minimal valid CellRanger H5 with 3 GEX + 2 ADT features and 4 cells.
  h5 <- tempfile(fileext = ".h5")
  on.exit(unlink(h5))
  n_feat <- 5; n_cells <- 4
  ids <- c(paste0("ENSG", 1:3), paste0("ADT", 1:2))
  names <- c(paste0("Gene", 1:3), paste0("ADT", 1:2))
  ftypes <- c(rep("Gene Expression", 3), rep("Antibody Capture", 2))
  bcs <- paste0("BC", 1:n_cells)
  # Dense counts matrix: 5 features x 4 cells
  M <- matrix(as.integer(c(1,2,3,4,5, 6,7,8,9,10, 11,12,13,14,15, 16,17,18,19,20)),
              nrow = n_feat, ncol = n_cells)
  # CSC sparse storage as 10x writes it
  sp <- Matrix::Matrix(M, sparse = TRUE)
  hf <- hdf5r::H5File$new(h5, mode = "w")
  g <- hf$create_group("matrix")
  g[["data"]] <- as.integer(sp@x)
  g[["indices"]] <- as.integer(sp@i)
  g[["indptr"]] <- as.integer(sp@p)
  g[["shape"]] <- as.integer(c(n_feat, n_cells))
  g[["barcodes"]] <- bcs
  fg <- g$create_group("features")
  fg[["id"]] <- ids
  fg[["name"]] <- names
  fg[["feature_type"]] <- ftypes
  hf$close_all()

  adt <- readADTCounts(h5, sampleName = "S1")
  expect_equal(nrow(adt), 2L)
  expect_setequal(rownames(adt), c("ADT1", "ADT2"))
  expect_equal(ncol(adt), n_cells)
})

test_that("processADT attaches ADT assay and adt.umap reduction to a Seurat object", {
  skip_on_cran(); skip_if_not_installed("Seurat")
  set.seed(1)
  rna <- Seurat::CreateSeuratObject(counts = matrix(rpois(2000, 5), nrow = 100, ncol = 20))
  adt_counts <- matrix(rpois(200, 50), nrow = 10, ncol = 20)
  rownames(adt_counts) <- paste0("ADT", 1:10)
  colnames(adt_counts) <- colnames(rna)
  obj <- processADT(rna, adt_counts, normMethod = "CLR")  # CLR for fast test
  expect_true("ADT" %in% Seurat::Assays(obj))
  expect_true("adt.umap" %in% names(obj@reductions))
})
