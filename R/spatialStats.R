###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## Spot-level spatial statistics for spot-based (Visium SD) data. All functions
## work on a sparse spot x spot adjacency matrix and are independent of Seurat
## so that they can be unit-tested on synthetic grids.

##' @title Hexagonal spot neighbour graph
##' @description Builds the spot adjacency of a Visium SD capture area from the
##'   array coordinates. On the Visium array a spot's direct neighbours sit in
##'   the same row two columns away or in an adjacent row one column away; with
##'   \code{x = array_col} and \code{y = sqrt(3) * array_row} all six of them
##'   are exactly 2 units away. \code{rings = 2} adds the next ring (18
##'   neighbours in total).
##' @param arrayRow,arrayCol integer array coordinates.
##' @param spotNames names of the spots, used as dimnames.
##' @param rings number of hexagonal rings to include.
##' @return a symmetric binary \code{dgCMatrix} without self-loops.
spatial_hex_graph <- function(arrayRow, arrayCol, spotNames, rings = 1) {
  stopifnot(length(arrayRow) == length(arrayCol))
  stopifnot(length(arrayRow) == length(spotNames))
  x <- as.numeric(arrayCol)
  y <- sqrt(3) * as.numeric(arrayRow)
  maxDist <- 2 * rings + 0.1
  ## rings = 2 has at most 18 neighbours; 3 rings at most 36
  k <- min(length(x), 3 * rings * (rings + 1) + 1)
  nn <- RANN::nn2(cbind(x, y), k = k, searchtype = "radius", radius = maxDist)
  idx <- nn$nn.idx
  from <- rep(seq_len(nrow(idx)), ncol(idx))
  to <- as.vector(idx)
  keep <- to > 0 & to != from
  adj <- Matrix::sparseMatrix(
    i = from[keep],
    j = to[keep],
    x = 1,
    dims = c(length(x), length(x)),
    dimnames = list(spotNames, spotNames)
  )
  adj <- (adj + Matrix::t(adj)) > 0
  return(methods::as(methods::as(adj, "dMatrix"), "CsparseMatrix"))
}

##' @title Row-standardise an adjacency matrix
##' @param adj sparse adjacency.
##' @param includeSelf add self-loops before standardising (used for smoothing).
##' @return sparse matrix whose rows sum to 1 (rows without neighbours stay 0).
spatial_row_standardise <- function(adj, includeSelf = FALSE) {
  if (includeSelf) {
    adj <- adj + Matrix::Diagonal(nrow(adj))
  }
  rs <- Matrix::rowSums(adj)
  rs[rs == 0] <- 1
  return(Matrix::Diagonal(x = 1 / rs) %*% adj)
}

##' @title Moran's I for many features
##' @description Global Moran's I of every row of \code{mat} with the
##'   analytical expectation and variance under the normality assumption
##'   (Cliff & Ord). Computed in gene chunks with sparse algebra, so all
##'   expressed genes can be tested in seconds.
##' @param mat feature x spot matrix (dense or sparse), e.g. SCT data.
##' @param adj spot x spot adjacency in the same spot order.
##' @param chunkSize number of genes per chunk.
##' @return data.frame with moransI, expected, zScore, pValue (one-sided,
##'   positive autocorrelation) and fdr, ordered by decreasing moransI.
spatial_morans_i <- function(mat, adj, chunkSize = 1000) {
  stopifnot(ncol(mat) == nrow(adj))
  w <- spatial_row_standardise(adj)
  n <- ncol(mat)
  s0 <- sum(w)
  wSym <- w + Matrix::t(w)
  s1 <- 0.5 * sum(wSym^2)
  s2 <- sum((Matrix::rowSums(w) + Matrix::colSums(w))^2)
  expected <- -1 / (n - 1)
  varI <- (n^2 * s1 - n * s2 + 3 * s0^2) / ((n^2 - 1) * s0^2) - expected^2
  chunks <- split(seq_len(nrow(mat)), ceiling(seq_len(nrow(mat)) / chunkSize))
  moransI <- unlist(lapply(chunks, function(rows) {
    z <- as.matrix(mat[rows, , drop = FALSE])
    z <- z - rowMeans(z)
    num <- rowSums(z * as.matrix(z %*% Matrix::t(w)))
    den <- rowSums(z^2)
    res <- (n / s0) * num / den
    res[den == 0] <- NA
    res
  }), use.names = FALSE)
  zScore <- (moransI - expected) / sqrt(varI)
  pValue <- stats::pnorm(zScore, lower.tail = FALSE)
  res <- data.frame(
    gene = rownames(mat),
    moransI = moransI,
    expected = expected,
    zScore = zScore,
    pValue = pValue,
    fdr = stats::p.adjust(pValue, method = "BH"),
    row.names = rownames(mat)
  )
  return(res[order(res$moransI, decreasing = TRUE, na.last = TRUE), ])
}

##' @title Neighbourhood smoothing of features
##' @param mat feature x spot matrix.
##' @param adj spot x spot adjacency.
##' @return dense feature x spot matrix with each spot replaced by the mean of
##'   itself and its neighbours.
spatial_smooth <- function(mat, adj) {
  w <- spatial_row_standardise(adj, includeSelf = TRUE)
  res <- as.matrix(mat %*% Matrix::t(w))
  dimnames(res) <- dimnames(mat)
  return(res)
}

##' @title Domain adjacency enrichment
##' @description Counts spot-spot edges between every pair of labels and
##'   compares them with label permutations. A positive z-score means two
##'   domains touch more often than expected from their sizes.
##' @param labels factor or character, one per spot (same order as adj).
##' @param adj spot x spot adjacency.
##' @param nPerm number of permutations.
##' @param seed random seed.
##' @return list with observed edge counts, expected counts and z-scores
##'   (label x label matrices).
spatial_adjacency_enrichment <- function(labels, adj, nPerm = 500, seed = 38) {
  labels <- droplevels(as.factor(labels))
  stopifnot(length(labels) == nrow(adj))
  lev <- levels(labels)
  count_edges <- function(lab) {
    ind <- Matrix::sparseMatrix(
      i = seq_along(lab),
      j = as.integer(lab),
      x = 1,
      dims = c(length(lab), length(lev))
    )
    m <- as.matrix(Matrix::t(ind) %*% adj %*% ind)
    dimnames(m) <- list(lev, lev)
    m
  }
  observed <- count_edges(labels)
  set.seed(seed)
  perms <- vapply(
    seq_len(nPerm),
    function(i) count_edges(factor(sample(labels), levels = lev)),
    observed
  )
  expected <- apply(perms, c(1, 2), mean)
  sdPerm <- apply(perms, c(1, 2), stats::sd)
  zScore <- (observed - expected) / sdPerm
  zScore[!is.finite(zScore)] <- 0
  return(list(observed = observed, expected = expected, zScore = zScore))
}

##' @title Spatial co-localisation of continuous scores
##' @description Plain per-spot correlation and the spatially lagged
##'   correlation (score of A in a spot vs the neighbourhood mean of B) for
##'   every pair of columns, e.g. RCTD cell-type weights.
##' @param scores spot x feature matrix.
##' @param adj spot x spot adjacency.
##' @return list with matrices \code{spot} and \code{lagged}.
spatial_colocalisation <- function(scores, adj) {
  scores <- as.matrix(scores)
  keep <- apply(scores, 2, stats::sd) > 0
  scores <- scores[, keep, drop = FALSE]
  lagged <- as.matrix(spatial_row_standardise(adj) %*% scores)
  spotCor <- stats::cor(scores)
  lagCor <- stats::cor(scores, lagged)
  lagCor <- (lagCor + t(lagCor)) / 2
  dimnames(lagCor) <- dimnames(spotCor)
  return(list(spot = spotCor, lagged = lagCor))
}

##' @title Spatial label smoothness
##' @param labels one label per spot.
##' @param adj spot x spot adjacency.
##' @return per-spot fraction of neighbours with the same label (NA for spots
##'   without neighbours).
spatial_label_smoothness <- function(labels, adj) {
  labels <- as.character(labels)
  a <- methods::as(adj, "TsparseMatrix")
  same <- labels[a@i + 1] == labels[a@j + 1]
  nSame <- tapply(same, factor(a@i + 1, levels = seq_along(labels)), sum)
  nAll <- Matrix::rowSums(adj)
  res <- as.numeric(nSame) / nAll
  res[nAll == 0] <- NA
  names(res) <- rownames(adj)
  return(res)
}

##' @title Adjusted Rand index
##' @param a,b two labelings of the same items.
##' @return numeric(1) in [-1, 1]; 1 means identical partitions.
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  comb2 <- function(x) x * (x - 1) / 2
  sumComb <- sum(comb2(tab))
  sumA <- sum(comb2(rowSums(tab)))
  sumB <- sum(comb2(colSums(tab)))
  expectedIndex <- sumA * sumB / comb2(sum(tab))
  maxIndex <- (sumA + sumB) / 2
  if (maxIndex == expectedIndex) {
    return(1)
  }
  return((sumComb - expectedIndex) / (maxIndex - expectedIndex))
}

##' @title Assign region labels to domains
##' @description Turns a domain x region score matrix (e.g. mean z-scored
##'   UCell scores) into one label per domain: the best region if its score is
##'   positive and exceeds the runner-up by \code{minMargin}, otherwise "mixed"
##'   followed by the two best regions.
##' @param scoreMat domain x region numeric matrix.
##' @param minMargin minimal difference between best and second best score.
##' @return data.frame with domain, label, best, second, margin.
assign_domain_labels <- function(scoreMat, minMargin = 0.25) {
  stopifnot(!is.null(rownames(scoreMat)), !is.null(colnames(scoreMat)))
  res <- lapply(rownames(scoreMat), function(d) {
    s <- sort(scoreMat[d, ], decreasing = TRUE)
    best <- names(s)[1]
    second <- if (length(s) > 1) names(s)[2] else NA_character_
    margin <- if (length(s) > 1) s[[1]] - s[[2]] else Inf
    label <- if (s[[1]] > 0 && margin >= minMargin) {
      best
    } else {
      paste0("mixed (", best, "/", second, ")")
    }
    data.frame(
      domain = d,
      label = label,
      best = best,
      second = second,
      margin = margin
    )
  })
  return(do.call(rbind, res))
}
