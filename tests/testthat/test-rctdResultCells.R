## The VisiumHD and SpatialSeurat apps read rctd-py results with rhdf5 because
## anndataR 1.2 drops the obs table of an anndata >= 0.13 file (job 386669).
## rctdResultCells is the guard that the rows still line up with the query.

.rctdEnv <- function() {
  e <- new.env()
  sys.source(test_path("../../R/rctdUtils.R"), envir = e)
  e
}

test_that("rctdResultCells accepts rows in query order", {
  e <- .rctdEnv()
  q <- c("s_016um_1-1", "s_016um_2-1", "s_016um_3-1")
  expect_identical(e$rctdResultCells(q, q), q)
})

test_that("rctdResultCells refuses a reordered or truncated result", {
  e <- .rctdEnv()
  q <- c("a", "b", "c")
  expect_error(e$rctdResultCells(c("b", "a", "c"), q), "query order")
  expect_error(e$rctdResultCells(c("a", "b"), q), "2 rows for 3")
})

test_that("readRctdPyResult marks the spots full mode left unfitted", {
  skip_if_not_installed("rhdf5")
  e <- .rctdEnv()
  f <- tempfile(fileext = ".h5ad")
  rhdf5::h5createFile(f)
  rhdf5::h5createGroup(f, "obs")
  rhdf5::h5createGroup(f, "obsm")
  rhdf5::h5createGroup(f, "uns")
  rhdf5::h5write(c("s1", "s2", "s3"), f, "obs/_index")
  ## rctd-py writes spots x types row-major; rhdf5 reads it transposed
  w <- rbind(c(0.2, 0.8), c(NaN, NaN), c(0.5, 0.5))
  rhdf5::h5write(t(w), f, "obsm/rctd_weights")
  rhdf5::h5write(c("Astro", "Neuron"), f, "uns/rctd_cell_type_names")
  res <- e$readRctdPyResult(f, mode = "full")
  expect_identical(res$obsNames, c("s1", "s2", "s3"))
  expect_identical(res$kept, c(TRUE, FALSE, TRUE))
  expect_equal(res$weights[1, ], c(0.2, 0.8))
  expect_identical(res$cellTypes, c("Astro", "Neuron"))
})

test_that("rctd-py compiles on a GPU or for large inputs only", {
  e <- .rctdEnv()
  withr::local_envvar(CUDA_VISIBLE_DEVICES = "")
  expect_false(e$rctdPyCompile(2319))
  expect_true(e$rctdPyCompile(99739))
  withr::local_envvar(CUDA_VISIBLE_DEVICES = "0")
  expect_true(e$rctdPyCompile(2319))
  withr::local_envvar(CUDA_VISIBLE_DEVICES = "NoDevFiles")
  expect_false(e$rctdPyCompile(2319))
})
