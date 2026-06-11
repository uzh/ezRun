## Tests for Xenium per-cell QC metric and MECR helpers.
## These functions are pure (operate on data.frame / matrix) so they can be
## tested without loading Xenium data or installing the package.

## When run outside an installed ezRun (e.g. during development), source the
## helper directly so the tests are self-contained.
if (!exists("computeXeniumQcMetrics")) {
  helper <- file.path("..", "..", "R", "xeniumQcMetrics.R")
  if (!file.exists(helper)) {
    helper <- file.path("R", "xeniumQcMetrics.R")
  }
  source(helper)
}

context("Xenium per-cell QC metrics")

test_that("derived per-cell metrics are computed correctly", {
  df <- data.frame(
    nCount_Xenium = c(100, 200, 90, 40),
    cell_area = c(50, 100, 45, 20),
    nucleus_area = c(25, 50, NA, 10),
    nucleus_count = c(1, 2, 0, 1),
    nCount_ControlProbe = c(5, 0, 5, 0),
    nCount_ControlCodeword = c(5, 0, 5, 0)
  )
  out <- computeXeniumQcMetrics(df, nmads = 3)

  # signal density = counts / area
  expect_equal(out$signal_density, c(2, 2, 2, 2))

  # nucleus:cell area ratio; NA nucleus -> 0
  expect_equal(out$nucleus_cell_ratio, c(0.5, 0.5, 0, 0.5))

  # transcripts per nucleus; nucleus_count == 0 -> NA
  expect_equal(out$transcripts_per_nucleus, c(100, 100, NA, 40))

  # negative-control fraction = controls / (counts + controls)
  expect_equal(out$neg_control_frac, c(10 / 110, 0, 10 / 100, 0))
})

test_that("nucleus_count outliers flag anucleate and multinucleate cells", {
  df <- data.frame(
    nCount_Xenium = rep(100, 5),
    cell_area = rep(50, 5),
    nucleus_area = rep(25, 5),
    nucleus_count = c(0, 1, 1, 2, 1)
  )
  out <- computeXeniumQcMetrics(df, nmads = 3)
  expect_equal(out$outlier.nucleus_count, c(TRUE, FALSE, FALSE, TRUE, FALSE))
})

test_that("cell area outliers are flagged via MAD on a clear outlier", {
  areas <- c(seq(45, 55, length.out = 20), 5000)
  df <- data.frame(
    nCount_Xenium = rep(100, 21),
    cell_area = areas,
    nucleus_area = rep(25, 21),
    nucleus_count = rep(1, 21)
  )
  out <- computeXeniumQcMetrics(df, nmads = 3)
  expect_true(out$outlier.cell_area[21])      # 5000 is a gross outlier
  expect_false(any(out$outlier.cell_area[1:20]))
})

test_that("missing optional columns are handled gracefully", {
  df <- data.frame(
    nCount_Xenium = c(100, 200),
    cell_area = c(50, 100)
  )
  out <- computeXeniumQcMetrics(df, nmads = 3)
  expect_true("signal_density" %in% colnames(out))
  expect_false("transcripts_per_nucleus" %in% colnames(out)) # no nucleus_count
  expect_false("neg_control_frac" %in% colnames(out))        # no control cols
})

context("Xenium MECR")

test_that("MECR computes Jaccard co-detection for marker pairs", {
  counts <- matrix(
    c(1, 1, 0, 0,   # gene A detected in cells 1,2
      0, 0, 1, 1,   # gene B detected in cells 3,4 (exclusive with A)
      1, 1, 1, 0),  # gene C detected in cells 1,2,3
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), paste0("c", 1:4))
  )
  res <- computeXeniumMECR(counts, marker_genes = c("A", "B", "C"))

  expect_equal(res$pairwise["A", "B"], 0)            # never co-detected
  expect_equal(res$pairwise["A", "C"], 2 / 3, tolerance = 1e-8)
  expect_equal(res$pairwise["B", "C"], 1 / 4, tolerance = 1e-8)
  expect_equal(res$mecr, mean(c(0, 2 / 3, 1 / 4)), tolerance = 1e-8)
})

test_that("MECR drops marker genes absent from the count matrix", {
  counts <- matrix(
    c(1, 0, 1, 0,
      0, 1, 1, 0),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("A", "B"), paste0("c", 1:4))
  )
  res <- computeXeniumMECR(counts, marker_genes = c("A", "B", "ZZZ_absent"))
  expect_setequal(res$markers_used, c("A", "B"))
  expect_equal(res$n_pairs, 1)
})
