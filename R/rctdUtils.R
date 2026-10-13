###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## RCTD (spacexr) helpers for spot-based data. XeniumSeurat and VisiumHDSeurat
## keep their own inline code paths; these functions follow the same rules
## (reference coercion, gene-overlap check, UMI_min_sigma fallback, fixed seed).

RCTD_REFERENCE_DIR <- "/srv/GT/databases/RCTD_References"

##' @title Resolve the RCTD reference path from the app parameters
##' @description \code{rctdFile} (full path) wins over \code{rctdReference},
##'   which is relative to \code{RCTD_REFERENCE_DIR} and may carry a
##'   " (description)" suffix as used in the SUSHI dropdowns.
##' @return a path or NULL when no reference was requested.
resolve_rctd_reference_path <- function(param) {
  if (ezIsSpecified(param$rctdFile) && param$rctdFile != "") {
    return(param$rctdFile)
  }
  if (ezIsSpecified(param$rctdReference) && param$rctdReference != "None") {
    refRelative <- sub(" \\([^)]+\\)$", "", param$rctdReference)
    return(file.path(RCTD_REFERENCE_DIR, refRelative))
  }
  return(NULL)
}

##' @title Load an RCTD reference
##' @description Accepts a spacexr \code{Reference}, a list holding one in
##'   \code{$reference}, or a Seurat object with raw counts and a cell-type
##'   column.
##' @param refPath path to an .rds or .qs2 file.
##' @param cellTypeColumn cell-type column used for Seurat references.
##' @return a spacexr \code{Reference}.
load_rctd_reference <- function(refPath, cellTypeColumn = NULL) {
  if (!file.exists(refPath)) {
    stop("RCTD reference not found: ", refPath)
  }
  refObj <- ezLoadRobj(refPath)
  if (inherits(refObj, "Reference")) {
    return(refObj)
  }
  if (is.list(refObj) && "reference" %in% names(refObj)) {
    return(refObj$reference)
  }
  if (inherits(refObj, "Seurat")) {
    candidates <- c(
      cellTypeColumn,
      "author_cell_type",
      "cell_type",
      "celltype",
      "CellType",
      "class",
      "subclass"
    )
    ctCol <- intersect(candidates, colnames(refObj@meta.data))[1]
    if (is.na(ctCol)) {
      stop("No cell-type column found in the Seurat reference ", refPath)
    }
    cellTypes <- as.factor(refObj@meta.data[[ctCol]])
    names(cellTypes) <- colnames(refObj)
    counts <- Seurat::GetAssayData(refObj, layer = "counts")
    return(spacexr::Reference(counts, cellTypes, Matrix::colSums(counts)))
  }
  stop("Not a spacexr Reference, list(reference=) or Seurat object: ", refPath)
}

##' @title Name the query genes like the reference
##' @description References built from atlases often use Ensembl IDs (e.g. the
##'   Allen brain references) while the Seurat object carries (uniquified)
##'   symbols. Renames the query rows to gene IDs when that gives the larger
##'   overlap; all genes are kept, so per-spot UMI totals do not change.
##' @param counts gene x spot count matrix with symbol row names.
##' @param geneIds gene IDs named by the row names of \code{counts}.
##' @param refGenes genes of the reference.
##' @return the count matrix, with gene IDs as row names if those match better.
rename_rctd_query_genes <- function(counts, geneIds, refGenes) {
  bySymbol <- sum(rownames(counts) %in% refGenes)
  ids <- unname(geneIds[rownames(counts)])
  byId <- sum(!is.na(ids) & ids %in% refGenes)
  if (byId > bySymbol) {
    keep <- !is.na(ids) & !duplicated(ids)
    counts <- counts[keep, , drop = FALSE]
    rownames(counts) <- ids[keep]
  }
  nShared <- max(byId, bySymbol)
  if (nShared < 100) {
    stop(sprintf(
      "RCTD reference shares only %d genes with the data (wrong species or ID type?)",
      nShared
    ))
  }
  ezLog(sprintf(
    "RCTD: %d genes shared with the reference (matched by %s)",
    nShared,
    if (byId > bySymbol) "gene ID" else "symbol"
  ))
  return(counts)
}

##' @title Match query gene names to the reference
##' @description \code{rename_rctd_query_genes} followed by restriction to the
##'   shared genes.
##' @inheritParams rename_rctd_query_genes
##' @return the count matrix restricted to genes shared with the reference.
match_rctd_genes <- function(counts, geneIds, refGenes) {
  counts <- rename_rctd_query_genes(counts, geneIds, refGenes)
  shared <- intersect(rownames(counts), refGenes)
  return(counts[shared, , drop = FALSE])
}

##' @title Main cell type per bin from a doublet-mode RCTD run
##' @description The cell type with the highest weight, except for bins that
##'   RCTD classified as \code{"reject"} (it could not fit them confidently),
##'   which get NA, as in XeniumSeurat.
##' @param weights bin x cell type weight matrix.
##' @param resultsDf RCTD results_df with a \code{spot_class} column, row
##'   names = bins; NULL keeps every bin.
##' @return named character vector, one entry per row of \code{weights}.
rctd_main_type <- function(weights, resultsDf = NULL) {
  mainType <- colnames(weights)[max.col(weights, ties.method = "first")]
  names(mainType) <- rownames(weights)
  if (!is.null(resultsDf) && "spot_class" %in% colnames(resultsDf)) {
    spotClass <- as.character(resultsDf[rownames(weights), "spot_class"])
    mainType[spotClass %in% "reject"] <- NA
  }
  return(mainType)
}

##' @title Run RCTD on spots
##' @param counts gene x spot raw counts (genes already matched).
##' @param coords data.frame with columns x and y, row names = spots.
##' @param reference spacexr Reference.
##' @param cores number of cores.
##' @param umiMin minimal UMIs per spot.
##' @param doubletMode "full" for multi-cell spots (Visium SD).
##' @param seed random seed (spacexr samples spots to fit sigma).
##' @return list with \code{weights} (spot x cell type, rows sum to 1) and
##'   \code{resultsDf} (NULL in full mode).
run_rctd_spots <- function(
  counts,
  coords,
  reference,
  cores = 1,
  umiMin = 100,
  doubletMode = "full",
  seed = 42
) {
  coords <- coords[colnames(counts), c("x", "y")]
  nUmi <- Matrix::colSums(counts)
  puck <- spacexr::SpatialRNA(coords, counts, nUmi)
  umiMinSigma <- 300
  if (sum(nUmi > umiMinSigma) < 100) {
    umiMinSigma <- max(1, floor(stats::quantile(nUmi, 0.75, names = FALSE)))
    ezLog("RCTD: few spots above 300 UMIs; UMI_min_sigma set to ", umiMinSigma)
  }
  ## spacexr starts fresh worker processes that do not inherit the BLAS
  ## thread pinning of this session: without these variables every worker
  ## uses all cores of the node (observed 8 workers x 17 threads, 4x slower)
  threadVars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
  oldVars <- Sys.getenv(threadVars, unset = NA)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(threadVars)), threadVars)))
  on.exit(
    {
      for (v in threadVars) {
        if (is.na(oldVars[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, stats::setNames(list(oldVars[[v]]), v))
      }
    },
    add = TRUE
  )
  set.seed(seed)
  rctd <- spacexr::create.RCTD(
    puck,
    reference,
    max_cores = cores,
    UMI_min = umiMin,
    UMI_min_sigma = umiMinSigma
  )
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = doubletMode)
  weights <- as.matrix(spacexr::normalize_weights(rctd@results$weights))
  return(list(weights = weights, resultsDf = rctd@results$results_df))
}
