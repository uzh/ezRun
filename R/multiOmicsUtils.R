###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Locate the CellRanger filtered feature H5 from a CountMatrix path.
##' @description The dataset.tsv `CountMatrix` column may point either at the matrix
##'   mtx directory (`*/sample_filtered_feature_bc_matrix/`, Read10X layout) or at
##'   a directory that already contains the H5 file. This helper returns the first
##'   matching `*filtered_feature_bc_matrix.h5` found in either the path itself or
##'   its parent directory.
##' @param countMatrixPath Path to the CountMatrix value.
##' @return Character path to the H5 file, or character(0) if none found.
##' @keywords internal
findFilteredH5 <- function(countMatrixPath) {
  candidates <- c(countMatrixPath, dirname(countMatrixPath))
  for (d in candidates) {
    if (!dir.exists(d)) next
    h5 <- list.files(d, pattern = "filtered_feature_bc_matrix\\.h5$",
                     full.names = TRUE, recursive = FALSE)
    if (length(h5) > 0) return(h5[1])
  }
  character(0)
}

##' @title Detect modalities present in a single-cell input directory.
##' @description Inspects a CellRanger / CellRanger Multi / CellRanger ARC output
##'   directory and reports which modalities are available for downstream
##'   processing. Robust to the CountMatrix path being either the matrix mtx
##'   directory (Read10X layout) or a count parent directory.
##' @param countMatrixPath Path to a directory containing or adjacent to
##'   CellRanger count outputs (the value of the dataset.tsv `CountMatrix` column).
##' @return Named list with logical flags: hasRNA, hasADT, hasVDJ_T, hasVDJ_B, hasATAC.
##' @export
detectModalities <- function(countMatrixPath) {
  stopifnot(dir.exists(countMatrixPath))

  h5 <- findFilteredH5(countMatrixPath)
  hasRNA <- length(h5) > 0

  hasADT <- hasRNA && h5HasAntibodyCapture(h5)

  # VDJ siblings can live one or two levels up from the CountMatrix.
  # Plain CellRanger:        <sample>/<count_or_outs>/vdj_t/...
  # CellRanger Multi (FGCZ): per_sample_outs/<sample>-cellRanger/{count,vdj_t}/...
  vdj_candidates <- unique(c(dirname(countMatrixPath),
                             dirname(dirname(countMatrixPath))))
  hasVDJ_T <- any(vapply(vdj_candidates, function(p) {
    file.exists(file.path(p, "vdj_t", "filtered_contig_annotations.csv"))
  }, logical(1)))
  hasVDJ_B <- any(vapply(vdj_candidates, function(p) {
    file.exists(file.path(p, "vdj_b", "filtered_contig_annotations.csv"))
  }, logical(1)))

  # ATAC siblings live alongside or just above the count matrix.
  atac_candidates <- unique(c(countMatrixPath, dirname(countMatrixPath)))
  hasATAC <- any(vapply(atac_candidates, function(p) {
    file.exists(file.path(p, "atac_fragments.tsv.gz")) &&
      file.exists(file.path(p, "atac_peaks.bed"))
  }, logical(1)))

  list(hasRNA = hasRNA, hasADT = hasADT,
       hasVDJ_T = hasVDJ_T, hasVDJ_B = hasVDJ_B,
       hasATAC = hasATAC)
}

##' @title Extract Antibody Capture (ADT) counts from a CellRanger H5.
##' @description Reads a multi-modal CellRanger filtered_feature_bc_matrix.h5 and
##'   returns the Antibody Capture submatrix with `<sample>_<barcode>` colnames so
##'   it aligns with a Seurat object built upstream (which prefixes barcodes by
##'   sample name). When the H5 has only a single feature type, that matrix is
##'   returned unchanged when it is Antibody Capture, otherwise an empty matrix.
##' @param h5path Path to the CellRanger HDF5.
##' @param sampleName Sample identifier used to prefix barcodes (matches
##'   `obj$Sample` / `obj$cellBarcode`).
##' @return Sparse matrix with ADT features as rows and prefixed cell barcodes
##'   as columns, or NULL if no ADT features are present.
##' @export
readADTCounts <- function(h5path, sampleName) {
  stopifnot(file.exists(h5path))
  cts <- Seurat::Read10X_h5(h5path, use.names = TRUE, unique.features = TRUE)
  if (is.list(cts)) {
    if (!"Antibody Capture" %in% names(cts)) return(NULL)
    adt <- cts[["Antibody Capture"]]
  } else {
    # Single feature type - inspect the H5 directly to confirm it's ADT.
    if (!h5HasAntibodyCapture(h5path)) return(NULL)
    adt <- cts
  }
  if (!is.null(adt) && ncol(adt) > 0) {
    colnames(adt) <- paste0(sampleName, "_", colnames(adt))
  }
  adt
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
