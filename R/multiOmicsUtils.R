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
