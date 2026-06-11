###################################################################
# Xenium cell-type spatial co-occurrence
#
# Permutation-based neighbourhood enrichment between cell types on a
# fixed-radius spatial graph. From a SINGLE label-permutation null it
# returns BOTH:
#   - log2enrich: Giotto cellProximityEnrichment-style effect size
#       log2((observed + 1) / (expected + 1)), abundance-controlled by
#       the permutation null (a rare type is not "enriched" just for
#       being rare).
#   - zscore: squidpy nhood_enrichment-style (observed - mean_perm)/sd_perm.
# Reporting both lets the two methods be cross-checked for agreement.
#
# Pure function operating on a coordinate matrix + a cell-type vector so
# it can be unit-tested without Xenium data.
###################################################################

##' @title Cell-type spatial co-occurrence (neighbourhood enrichment)
##' @description Builds a fixed-radius spatial neighbour graph and tests,
##'   per cell-type pair, whether the two types are neighbours more (enriched)
##'   or less (avoidant) than under random label assignment. The label
##'   permutation keeps the graph fixed, so it controls for cell-type
##'   abundance and overall density.
##' @param coords numeric matrix / data.frame of cell coordinates; the first
##'   two columns are used as x, y (microns for Xenium).
##' @param celltypes vector (factor/character) of cell-type labels aligned to
##'   the rows of \code{coords}. \code{NA} labels are dropped.
##' @param radius neighbour radius in coordinate units (default 30 = 30 um).
##' @param n_perm number of label permutations for the null (default 1000).
##' @param min_cells drop cell types with fewer than this many cells
##'   (default 20).
##' @param max_cells if more cells than this, randomly subsample to keep the
##'   permutation tractable; graph + observed + null are all computed on the
##'   subsample (default 1e5).
##' @param seed RNG seed for subsampling + permutation (default 42).
##' @return list with K x K matrices \code{observed}, \code{expected},
##'   \code{sd}, \code{log2enrich}, \code{zscore}, \code{pval} (two-sided
##'   empirical), plus \code{levels}, \code{n_cells}, \code{n_total},
##'   \code{n_edges}, \code{mean_degree}, \code{radius}, \code{n_perm},
##'   \code{dropped_types}. Matrices are NULL when fewer than two cell types
##'   remain.
##' @export
computeCelltypeCooccurrence <- function(coords, celltypes, radius = 30,
                                        n_perm = 1000, min_cells = 20,
                                        max_cells = 1e5, seed = 42) {
  coords <- as.matrix(coords[, 1:2, drop = FALSE])
  storage.mode(coords) <- "double"
  ct <- as.character(celltypes)

  ## --- drop cells with missing label or coordinates ---
  keep <- !is.na(ct) & stats::complete.cases(coords)
  coords <- coords[keep, , drop = FALSE]
  ct <- ct[keep]

  ## --- drop rare cell types ---
  tab <- table(ct)
  rare <- names(tab)[tab < min_cells]
  dropped_types <- rare
  if (length(rare) > 0) {
    keep2 <- !(ct %in% rare)
    coords <- coords[keep2, , drop = FALSE]
    ct <- ct[keep2]
  }
  n_total <- length(ct)

  empty_result <- function(levels) {
    list(
      observed = NULL, expected = NULL, sd = NULL, log2enrich = NULL,
      zscore = NULL, pval = NULL, levels = levels, n_cells = length(ct),
      n_total = n_total, n_edges = 0, mean_degree = 0, radius = radius,
      n_perm = n_perm, dropped_types = dropped_types
    )
  }
  if (length(unique(ct)) < 2) {
    return(empty_result(unique(ct)))
  }

  ## --- subsample for tractability (graph + null on the same cells) ---
  set.seed(seed)
  if (n_total > max_cells) {
    idx <- sample.int(n_total, max_cells)
    coords <- coords[idx, , drop = FALSE]
    ct <- ct[idx]
  }
  n <- length(ct)
  ct <- factor(ct)
  levels_ct <- levels(ct)
  K <- length(levels_ct)

  ## --- fixed-radius spatial neighbour graph (kd-tree) ---
  nn <- dbscan::frNN(coords, eps = radius)
  deg <- lengths(nn$id)
  from <- rep.int(seq_len(n), deg)
  to <- unlist(nn$id, use.names = FALSE)
  n_edges <- length(to) # directed (frNN is symmetric -> each undirected edge twice)
  if (n_edges == 0) {
    res <- empty_result(levels_ct)
    res$n_cells <- n
    return(res)
  }

  A <- Matrix::sparseMatrix(i = from, j = to, x = 1, dims = c(n, n))
  ## cell x celltype indicator (n x K)
  X <- Matrix::t(Matrix::fac2sparse(ct))

  pair_counts <- function(ind) {
    as.matrix(Matrix::crossprod(ind, A %*% ind)) # t(ind) %*% (A %*% ind)
  }
  symm <- function(m) (m + t(m)) / 2

  observed <- symm(pair_counts(X))

  ## --- permutation null (streaming mean / sd / tail counts) ---
  s <- matrix(0, K, K)
  ss <- matrix(0, K, K)
  ge <- matrix(0L, K, K)
  le <- matrix(0L, K, K)
  for (p in seq_len(n_perm)) {
    Xp <- X[sample.int(n), , drop = FALSE]
    Cp <- symm(pair_counts(Xp))
    s <- s + Cp
    ss <- ss + Cp * Cp
    ge <- ge + (Cp >= observed)
    le <- le + (Cp <= observed)
  }
  expected <- s / n_perm
  v <- (ss - n_perm * expected^2) / (n_perm - 1)
  v[v < 0] <- 0
  sdm <- sqrt(v)

  log2enrich <- log2((observed + 1) / (expected + 1))
  zscore <- (observed - expected) / sdm
  zscore[!is.finite(zscore)] <- 0 # sd == 0 -> no variation -> z = 0
  p_gt <- (ge + 1) / (n_perm + 1)
  p_lt <- (le + 1) / (n_perm + 1)
  ## pmin() drops the matrix dim -> rebuild as a K x K matrix
  pval <- matrix(pmin(1, 2 * pmin(p_gt, p_lt)), nrow = K, ncol = K)

  dn <- list(levels_ct, levels_ct)
  for (m in c("observed", "expected", "sdm", "log2enrich", "zscore", "pval")) {
    tmp <- get(m)
    dimnames(tmp) <- dn
    assign(m, tmp)
  }

  list(
    observed = observed, expected = expected, sd = sdm,
    log2enrich = log2enrich, zscore = zscore, pval = pval,
    levels = levels_ct, n_cells = n, n_total = n_total,
    n_edges = n_edges, mean_degree = mean(deg), radius = radius,
    n_perm = n_perm, dropped_types = dropped_types
  )
}
