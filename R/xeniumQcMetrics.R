###################################################################
# Xenium per-cell QC metrics and MECR
#
# Helpers for the XeniumSeurat app. Kept as pure functions operating on
# a data.frame / matrix so they can be unit-tested without Xenium data.
#
# - computeXeniumQcMetrics(): derives per-cell QC metrics (signal density,
#   nucleus:cell area ratio, transcripts per nucleus, negative-control
#   fraction) and MAD-based outlier flags. Outlier flags use the "outlier."
#   prefix so they are reported but NOT folded into the "qc."/discard
#   removal logic (cells are flagged, not removed).
# - addXeniumCellQc(): thin Seurat wrapper around computeXeniumQcMetrics().
# - computeXeniumMECR(): Mutually Exclusive Co-expression Rate over a set of
#   marker genes (Jaccard co-detection per pair); a contamination / doublet
#   QC metric in the spirit of SpatialQM.
###################################################################

## Convert an isOutlier() result to a plain logical with NA -> FALSE.
.flagFalseNA <- function(x) {
  x <- as.logical(x)
  x[is.na(x)] <- FALSE
  x
}

##' @title Compute Xenium per-cell QC metrics and outlier flags
##' @description Derives per-cell QC metrics not already present and flags
##'   MAD-based outliers. Optional source columns are handled gracefully:
##'   a metric is only added when its inputs exist.
##' @param df data.frame of per-cell metadata. Recognised columns:
##'   a count column matching "^nCount_Xenium" (or supplied via \code{count_col}),
##'   \code{cell_area}, \code{nucleus_area}, \code{nucleus_count} and any of
##'   \code{nCount_ControlProbe}, \code{nCount_ControlCodeword},
##'   \code{nCount_BlankCodeword}.
##' @param nmads numeric, MAD multiple for outlier detection (default 3).
##' @param count_col optional explicit transcript-count column name.
##' @return The input data.frame with added metric and \code{outlier.*} columns.
##' @export
computeXeniumQcMetrics <- function(df, nmads = 3, count_col = NULL) {
  if (is.null(count_col)) {
    hit <- grep("^nCount_Xenium", colnames(df), value = TRUE)
    count_col <- if (length(hit) > 0) hit[1] else NULL
  }
  has_counts <- !is.null(count_col) && count_col %in% colnames(df)
  counts <- if (has_counts) df[[count_col]] else NULL

  ## --- derived per-cell metrics ---
  if (has_counts && "cell_area" %in% colnames(df)) {
    sig <- counts / df$cell_area
    sig[!is.finite(sig)] <- NA
    df$signal_density <- sig
  }
  if (all(c("nucleus_area", "cell_area") %in% colnames(df))) {
    ratio <- df$nucleus_area / df$cell_area
    ratio[!is.finite(ratio)] <- 0 # NA / no nucleus -> 0
    df$nucleus_cell_ratio <- ratio
  }
  if (has_counts && "nucleus_count" %in% colnames(df)) {
    tpn <- counts / df$nucleus_count
    tpn[!is.finite(tpn)] <- NA # nucleus_count == 0 -> NA
    df$transcripts_per_nucleus <- tpn
  }
  ctrl_cols <- intersect(
    c("nCount_ControlProbe", "nCount_ControlCodeword", "nCount_BlankCodeword"),
    colnames(df)
  )
  if (has_counts && length(ctrl_cols) > 0) {
    ctrl <- rowSums(as.matrix(df[, ctrl_cols, drop = FALSE]), na.rm = TRUE)
    frac <- ctrl / (counts + ctrl)
    frac[!is.finite(frac)] <- 0
    df$neg_control_frac <- frac
  }

  ## --- outlier flags (flagged & reported, NOT removed) ---
  if ("cell_area" %in% colnames(df)) {
    df$outlier.cell_area <- .flagFalseNA(scater::isOutlier(
      df$cell_area,
      log = TRUE, type = "both", nmads = nmads
    ))
  }
  if ("signal_density" %in% colnames(df)) {
    sig <- df$signal_density
    sig[!is.finite(sig) | sig <= 0] <- NA
    df$outlier.signal_density <- .flagFalseNA(scater::isOutlier(
      sig,
      log = TRUE, type = "lower", nmads = nmads
    ))
  }
  if ("neg_control_frac" %in% colnames(df)) {
    df$outlier.neg_control <- .flagFalseNA(scater::isOutlier(
      df$neg_control_frac,
      type = "higher", nmads = nmads
    ))
  }
  if ("nucleus_count" %in% colnames(df)) {
    # anucleate (0) or multinucleate (>1) cells are segmentation concerns
    df$outlier.nucleus_count <- df$nucleus_count == 0 | df$nucleus_count > 1
  }
  df
}

##' @title Add Xenium per-cell QC metrics to a Seurat object
##' @param object a Seurat object whose meta.data carries the Xenium per-cell
##'   columns (see \code{computeXeniumQcMetrics}).
##' @param nmads numeric, MAD multiple for outlier detection (default 3).
##' @return The Seurat object with the new metric / outlier columns added.
##' @export
addXeniumCellQc <- function(object, nmads = 3) {
  md <- computeXeniumQcMetrics(object@meta.data, nmads = nmads)
  new_cols <- setdiff(colnames(md), colnames(object@meta.data))
  for (cl in new_cols) {
    object[[cl]] <- md[[cl]]
  }
  object
}

##' @title Mutually Exclusive Co-expression Rate (MECR)
##' @description For each pair of marker genes, the Jaccard co-detection rate
##'   = (cells expressing both) / (cells expressing either). Markers that label
##'   distinct populations are expected to be mutually exclusive, so a high MECR
##'   indicates transcript spill-over / segmentation contamination.
##' @param counts gene x cell count matrix (dense or sparse).
##' @param marker_genes character vector of marker genes; those absent from
##'   \code{counts} are dropped.
##' @param detection_threshold counts strictly above this are "detected" (0).
##' @return list with \code{pairwise} (marker x marker rate matrix),
##'   \code{mecr} (mean over unique pairs), \code{n_pairs}, \code{markers_used}.
##' @export
computeXeniumMECR <- function(counts, marker_genes, detection_threshold = 0) {
  markers_used <- intersect(unique(marker_genes), rownames(counts))
  n <- length(markers_used)
  if (n < 2) {
    return(list(
      pairwise = NULL, mecr = NA_real_, n_pairs = 0,
      markers_used = markers_used
    ))
  }
  detected <- as.matrix(counts[markers_used, , drop = FALSE]) > detection_threshold
  pairwise <- matrix(
    NA_real_, n, n,
    dimnames = list(markers_used, markers_used)
  )
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      both <- sum(detected[i, ] & detected[j, ])
      either <- sum(detected[i, ] | detected[j, ])
      pairwise[i, j] <- if (either > 0) both / either else 0
    }
  }
  upper <- pairwise[upper.tri(pairwise)]
  list(
    pairwise = pairwise,
    mecr = mean(upper, na.rm = TRUE),
    n_pairs = sum(upper.tri(pairwise)),
    markers_used = markers_used
  )
}
