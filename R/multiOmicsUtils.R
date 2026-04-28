###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Detect modalities present in a single-cell input directory.
##' @description Inspects a CellRanger / CellRanger Multi / CellRanger ARC output
##'   directory and reports which modalities are available for downstream
##'   processing.
##' @param countMatrixPath Path to a directory containing CellRanger count outputs
##'   (the value of the dataset.tsv `CountMatrix` column).
##' @return Named list with logical flags: hasRNA, hasADT, hasVDJ_T, hasVDJ_B, hasATAC.
##' @export
detectModalities <- function(countMatrixPath) {
  stopifnot(dir.exists(countMatrixPath))

  h5_filtered <- list.files(countMatrixPath, pattern = "filtered_feature_bc_matrix\\.h5$",
                            full.names = TRUE, recursive = FALSE)
  hasRNA <- length(h5_filtered) > 0

  hasADT <- hasRNA && h5HasAntibodyCapture(h5_filtered[1])

  parent <- dirname(countMatrixPath)  # e.g. per_sample_outs/<sample>/
  hasVDJ_T <- file.exists(file.path(parent, "vdj_t", "filtered_contig_annotations.csv"))
  hasVDJ_B <- file.exists(file.path(parent, "vdj_b", "filtered_contig_annotations.csv"))

  hasATAC <- file.exists(file.path(countMatrixPath, "atac_fragments.tsv.gz")) &&
             file.exists(file.path(countMatrixPath, "atac_peaks.bed"))

  list(hasRNA = hasRNA, hasADT = hasADT,
       hasVDJ_T = hasVDJ_T, hasVDJ_B = hasVDJ_B,
       hasATAC = hasATAC)
}

##' @title Add an ADT assay to a Seurat object with normalization, PCA and UMAP.
##' @description Attaches a Seurat v5 ADT assay built from raw antibody counts,
##'   normalizes (ADTnorm or CLR), runs PCA on the ADT features and a UMAP
##'   embedding stored as `adt.umap`. Cells absent from `adtCounts` are dropped.
##' @param obj Seurat object (already has RNA assay + clusters).
##' @param adtCounts Sparse matrix or matrix of ADT raw counts (features x cells).
##'   Cell barcodes must match `colnames(obj)`.
##' @param normMethod "ADTnorm" (default) or "CLR".
##' @param npcs Number of ADT PCs (clipped to nrow(adt) - 1).
##' @return Seurat object with `assays$ADT`, `reductions$adt.pca`, `reductions$adt.umap`.
##'   `DefaultAssay` is restored to "RNA" before return.
##' @export
processADT <- function(obj, adtCounts, normMethod = c("ADTnorm", "CLR"), npcs = 18) {
  normMethod <- match.arg(normMethod)
  shared <- intersect(colnames(adtCounts), colnames(obj))
  if (length(shared) < ncol(obj) * 0.5) {
    warning(sprintf("Only %d/%d cells overlap ADT barcodes - check sample prefixing.",
                    length(shared), ncol(obj)))
  }
  obj <- obj[, shared]
  adt_aligned <- adtCounts[, shared, drop = FALSE]

  obj[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = adt_aligned)
  Seurat::DefaultAssay(obj) <- "ADT"

  if (normMethod == "ADTnorm" && requireNamespace("ADTnorm", quietly = TRUE)) {
    norm_mat <- ADTnorm::ADTnorm(
      cell_x_adt = t(as.matrix(adt_aligned)),
      cell_x_feature = data.frame(sample = rep("all", ncol(adt_aligned))),
      marker_to_process = rownames(adt_aligned),
      save_outpath = tempdir(),
      study_name = "ezrun",
      save_intermediate_fig = FALSE
    )
    obj <- Seurat::SetAssayData(obj, assay = "ADT", layer = "data",
                                new.data = t(norm_mat))
  } else {
    obj <- Seurat::NormalizeData(obj, assay = "ADT",
                                 normalization.method = "CLR", margin = 2,
                                 verbose = FALSE)
  }

  npcs <- min(npcs, nrow(adt_aligned) - 1)
  obj <- Seurat::ScaleData(obj, assay = "ADT", verbose = FALSE)
  obj <- Seurat::RunPCA(obj, assay = "ADT", reduction.name = "adt.pca",
                        npcs = npcs, features = rownames(obj[["ADT"]]),
                        approx = FALSE, verbose = FALSE)
  n_neighbors <- min(30L, max(2L, ncol(obj) - 1L))
  obj <- Seurat::RunUMAP(obj, assay = "ADT", reduction = "adt.pca",
                         dims = seq_len(npcs), reduction.name = "adt.umap",
                         n.neighbors = n_neighbors, verbose = FALSE)

  Seurat::DefaultAssay(obj) <- "RNA"
  obj
}

##' @title Check whether a CellRanger HDF5 has any Antibody Capture features.
##' @param h5path Path to a CellRanger filtered_feature_bc_matrix.h5 file.
##' @return TRUE if any feature has type "Antibody Capture", FALSE otherwise
##'   (or if the file is missing / empty / unreadable).
##' @keywords internal
h5HasAntibodyCapture <- function(h5path) {
  if (!file.exists(h5path)) return(FALSE)
  if (file.size(h5path) == 0) return(FALSE)
  tryCatch({
    hf <- hdf5r::H5File$new(h5path, mode = "r")
    on.exit(hf$close_all(), add = TRUE)
    if (!hf$exists("matrix/features/feature_type")) return(FALSE)
    ft <- hf[["matrix/features/feature_type"]]$read()
    any(ft == "Antibody Capture")
  }, error = function(e) FALSE)
}
