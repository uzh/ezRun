make_counts <- function(genes, nSpots = 5) {
  m <- Matrix::Matrix(1, nrow = length(genes), ncol = nSpots, sparse = TRUE)
  dimnames(m) <- list(genes, paste0("s", seq_len(nSpots)))
  m
}

test_that("query genes are renamed to IDs only when IDs match better", {
  symbols <- paste0("Gene", 1:200)
  ids <- stats::setNames(sprintf("ENSMUSG%011d", 1:200), symbols)
  counts <- make_counts(symbols)

  byId <- rename_rctd_query_genes(counts, ids, unname(ids[1:150]))
  expect_equal(rownames(byId), unname(ids))
  ## all genes are kept, so UMI totals per spot stay the same
  expect_equal(Matrix::colSums(byId), Matrix::colSums(counts))

  bySymbol <- rename_rctd_query_genes(counts, ids, symbols[1:150])
  expect_equal(rownames(bySymbol), symbols)

  expect_error(rename_rctd_query_genes(counts, ids, c("x", "y")), "shares only")
})

test_that("rejected bins get no main cell type", {
  w <- rbind(
    b1 = c(A = 0.7, B = 0.3),
    b2 = c(A = 0.2, B = 0.8),
    b3 = c(A = 0.6, B = 0.4)
  )
  rdf <- data.frame(
    spot_class = factor(c("singlet", "doublet_certain", "reject")),
    row.names = c("b1", "b2", "b3")
  )
  expect_equal(rctd_main_type(w, rdf), c(b1 = "A", b2 = "B", b3 = NA))
  ## without results_df every bin keeps its highest-weight type
  expect_equal(rctd_main_type(w), c(b1 = "A", b2 = "B", b3 = "A"))
})

test_that("match_rctd_genes keeps only the shared genes", {
  symbols <- paste0("Gene", 1:200)
  ids <- stats::setNames(sprintf("ENSMUSG%011d", 1:200), symbols)
  m <- match_rctd_genes(make_counts(symbols), ids, unname(ids[51:200]))
  expect_equal(rownames(m), unname(ids[51:200]))
})
