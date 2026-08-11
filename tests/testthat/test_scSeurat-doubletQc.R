context("ScSeurat: zero-count cells must not reach scDblFinder")

## Regression test for the silent loss of doublet filtering seen in production
## (2026-08-11):
##
##   WARN scDblFinder failed (likely too few cells), skipping doublet detection:
##     error in evaluating the argument 'x' in selecting a method for function
##     't': x / y: operation not supported on SparseArray object x when y
##     contains zeros (result wouldn't be sparse)
##
## The message is misleading -- the sample had plenty of cells. scDblFinder
## builds artificial doublets in scDblFinder:::createDoublets and rescales them
## with `x2 %*% diag(ls/colSums(x2))`, falling back to `t(t(x2)/colSums(x2))`.
## Both blow up once a parent cell has zero counts. Empty cells reach that call
## because addCellQcToSeurat only bounds the library size when either a fixed
## nUMI threshold or nmad is set; with both unspecified qc.lib is a flat FALSE.
##
## The knock-on effect is what users notice: doublet detection is skipped, the
## error is downgraded to a warning, and -- since doublets are the only filter
## left in that configuration -- the QC table reports that nothing was removed.

makeScData <- function(nEmpty = 3, nGood = 40, nGenes = 200) {
  set.seed(7)
  counts <- matrix(
    rpois(nGenes * (nGood + nEmpty), 3),
    nrow = nGenes,
    dimnames = list(
      paste0("g", seq_len(nGenes)),
      paste0("c", seq_len(nGood + nEmpty))
    )
  )
  counts[, seq_len(nEmpty)] <- 0L
  Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
}

## no fixed thresholds and no nmad -- the configuration that let empty cells
## through in production
noThresholdParam <- list(
  nUMI = NA,
  ngenes = NA,
  perc_mito = NA,
  perc_riboprot = NA,
  nmad = NA,
  keepDoublets = FALSE
)

test_that("scDblFinder chokes on empty cells (the bug)", {
  skip_if_not_installed("scDblFinder")
  skip_if_not_installed("SparseArray")

  set.seed(7)
  counts <- matrix(
    rpois(20 * 6, 3), nrow = 20,
    dimnames = list(paste0("g", 1:20), paste0("c", 1:6))
  )
  counts[, 5:6] <- 0L

  # the two empty cells are paired into one artificial doublet
  dblIdx <- matrix(c(5L, 6L, 1L, 2L, 3L, 4L), ncol = 2, byrow = TRUE)
  clusters <- factor(c(1, 1, 2, 2, 1, 2))

  expect_error(
    suppressWarnings(scDblFinder:::createDoublets(
      as(counts, "SVT_SparseMatrix"),
      dblIdx,
      clusters = clusters,
      adjustSize = 1,
      halfSize = 0
    )),
    "SparseArray"
  )
})

test_that("addCellQcToSeurat flags zero-count cells even with no thresholds set", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scater")

  scData <- makeScData()
  isEmpty <- Matrix::colSums(
    Seurat::GetAssayData(scData, layer = "counts")
  ) == 0
  expect_equal(sum(isEmpty), 3)

  scData <- suppressWarnings(addCellQcToSeurat(
    scData,
    param = noThresholdParam,
    BPPARAM = BiocParallel::SerialParam()
  ))

  expect_equal(unname(scData$qc.lib), unname(isEmpty))
  expect_false(any(scData$useCell[isEmpty]))
})

test_that("the cells handed to scDblFinder all have counts", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scater")

  scData <- makeScData()
  scData <- suppressWarnings(addCellQcToSeurat(
    scData,
    param = noThresholdParam,
    BPPARAM = BiocParallel::SerialParam()
  ))

  # reproduce the subset addCellQcToSeurat passes to scDblFinder: useCell as it
  # stands before the doublet class narrows it further
  passedQc <- !(scData$qc.lib |
    scData$qc.nexprs |
    scData$qc.mito |
    scData$qc.riboprot)
  dblCounts <- Seurat::GetAssayData(scData, layer = "counts")[, passedQc]

  expect_gt(min(Matrix::colSums(dblCounts)), 0)
})

test_that("a skipped doublet call leaves the other cells in place", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scater")

  # too few cells for scDblFinder -- the genuine "too few cells" path
  scData <- makeScData(nEmpty = 1, nGood = 5, nGenes = 30)
  scData <- suppressWarnings(addCellQcToSeurat(
    scData,
    param = noThresholdParam,
    BPPARAM = BiocParallel::SerialParam()
  ))

  expect_true(all(scData$doubletClass == "undetermined"))
  expect_false(any(scData$qc.doublet))
  # the fallback must not discard everything
  expect_equal(sum(scData$useCell), 5)
})
