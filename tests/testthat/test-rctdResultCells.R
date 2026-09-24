## The VisiumHD app reads rctd-py results with rhdf5 because anndataR 1.2 drops
## the obs table of an anndata >= 0.13 file (job 386669). rctdResultCells is the
## guard that the rows still line up with the query bins.

.appEnv <- function() {
  e <- new.env()
  f <- test_path("../../R/app-VisiumHDSeurat.R")
  for (x in parse(f)) {
    if (!grepl("setRefClass", paste(deparse(x), collapse = ""))) eval(x, e)
  }
  e
}

test_that("rctdResultCells accepts rows in query order", {
  e <- .appEnv()
  q <- c("s_016um_1-1", "s_016um_2-1", "s_016um_3-1")
  expect_identical(e$rctdResultCells(q, q), q)
})

test_that("rctdResultCells refuses a reordered or truncated result", {
  e <- .appEnv()
  q <- c("a", "b", "c")
  expect_error(e$rctdResultCells(c("b", "a", "c"), q), "query order")
  expect_error(e$rctdResultCells(c("a", "b"), q), "2 rows for 3")
})
