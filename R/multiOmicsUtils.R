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

  # ADT detection requires reading the H5; populated in detectModalities once
  # h5HasAntibodyCapture() is wired in (see Task 1.2).
  hasADT <- FALSE

  parent <- dirname(countMatrixPath)  # e.g. per_sample_outs/<sample>/
  hasVDJ_T <- file.exists(file.path(parent, "vdj_t", "filtered_contig_annotations.csv"))
  hasVDJ_B <- file.exists(file.path(parent, "vdj_b", "filtered_contig_annotations.csv"))

  hasATAC <- file.exists(file.path(countMatrixPath, "atac_fragments.tsv.gz")) &&
             file.exists(file.path(countMatrixPath, "atac_peaks.bed"))

  list(hasRNA = hasRNA, hasADT = hasADT,
       hasVDJ_T = hasVDJ_T, hasVDJ_B = hasVDJ_B,
       hasATAC = hasATAC)
}
