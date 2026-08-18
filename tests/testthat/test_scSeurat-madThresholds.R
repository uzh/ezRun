context("ScSeurat: a MAD threshold on a bounded metric can be unreachable")

library(Seurat)

## percent_mito is a percentage, so it cannot exceed 100. MAD-based filtering
## uses median + nmads * MAD, which on a badly degraded sample lands ABOVE 100 --
## at which point the filter is mathematically incapable of flagging anything,
## and says so nowhere.
##
## Production case (o41361 / p31662, 2026-08-17). CellBender over-called cells
## 2.2-4.2x on every sample of the order, so on the degraded ones the ambient
## droplets are the majority and set the median. Measured from each sample's own
## allCellsMeta.qs2:
##
##   sample  median mito  MAD threshold  cells flagged
##   hIr51        62.4 %          229.6              0
##   hIr53        29.6 %          157.0              0
##   hCB54        74.8 %          187.0              0
##   pUM32        99.8 %          100.6              0
##
## Four of twelve samples had ZERO cells removed for mitochondrial content while
## the other eight lost 20-30 % of their barcodes to it, so no two samples in the
## batch were filtered comparably. In the report a threshold of 187 % is
## indistinguishable from a threshold that found nothing to flag, which is why
## this went unnoticed.

## Build an object with an exact per-cell mitochondrial fraction: every cell gets
## the same 1000 counts, split between MT- genes and the rest.
makeMitoScData <- function(mitoFrac, nGenes = 200, total = 1000) {
  nCells <- length(mitoFrac)
  nMito <- 20
  counts <- matrix(
    0L,
    nrow = nGenes,
    ncol = nCells,
    dimnames = list(
      c(paste0("MT-", seq_len(nMito)), paste0("g", seq_len(nGenes - nMito))),
      paste0("c", seq_len(nCells))
    )
  )
  for (j in seq_len(nCells)) {
    mt <- round(total * mitoFrac[j])
    counts[seq_len(nMito), j] <- as.integer(diff(round(seq(0, mt, length.out = nMito + 1))))
    rest <- total - mt
    idx <- (nMito + 1):nGenes
    counts[idx, j] <- as.integer(diff(round(seq(0, rest, length.out = length(idx) + 1))))
  }
  Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
}

madParam <- list(
  nUMI = NA,
  ngenes = NA,
  perc_mito = NA,
  perc_riboprot = NA,
  nmad = 3,
  keepDoublets = FALSE
)

## doublet detection is irrelevant here, so keep it out of the way and fast
withFakeDoublets <- function(env = parent.frame()) {
  localFakeScDblFinder(function(x, ...) fakeDoubletTable(colnames(x)), env = env)
}

test_that("an unreachable percent_mito MAD threshold is reported, not silent", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("withr")
  withFakeDoublets()

  ## Mimics hCB54's shape: a low-mito quarter, a mid band, and a high-mito half.
  ## That gives median 87.5 % and a MAD threshold of 143 %.
  ## Note a flat "most cells at 100 %" does NOT reproduce it -- the median goes to
  ## 100 and the MAD collapses to 0, so the threshold is exactly 100 and thus
  ## still reachable. The spread is what pushes it out of range.
  frac <- c(rep(0.05, 50), seq(0.5, 0.75, length.out = 50), rep(1, 100))
  scData <- makeMitoScData(frac)

  # addCellQcToSeurat computes percent_mito itself; do it here too so the
  # threshold can be checked before the call
  scData <- Seurat::PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
  th <- attr(
    scuttle::isOutlier(scData$percent_mito, nmads = 3, type = "higher"),
    "thresholds"
  )[["higher"]]
  expect_gt(th, 100) # the construction really is pathological

  logs <- captureEzLog(
    scData <- suppressWarnings(addCellQcToSeurat(
      scData,
      param = madParam,
      BPPARAM = BiocParallel::SerialParam()
    ))
  )

  # the pathology itself: nothing can be flagged
  expect_false(any(scData$qc.mito))
  # and it must be said out loud
  expect_true(any(grepl("percent_mito", logs, fixed = TRUE)))
  expect_true(any(grepl("cannot flag", logs, fixed = TRUE)))
})

test_that("a reachable percent_mito MAD threshold flags cells and warns nothing", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("withr")
  withFakeDoublets()

  ## a healthy sample: most cells low mito, a handful genuinely high
  set.seed(11)
  frac <- c(runif(190, 0.01, 0.08), runif(10, 0.5, 0.9))
  scData <- makeMitoScData(frac)

  # addCellQcToSeurat computes percent_mito itself; do it here too so the
  # threshold can be checked before the call
  scData <- Seurat::PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
  th <- attr(
    scuttle::isOutlier(scData$percent_mito, nmads = 3, type = "higher"),
    "thresholds"
  )[["higher"]]
  expect_lt(th, 100)

  logs <- captureEzLog(
    scData <- suppressWarnings(addCellQcToSeurat(
      scData,
      param = madParam,
      BPPARAM = BiocParallel::SerialParam()
    ))
  )

  # negative control: the warning must not fire when the threshold is usable
  expect_true(any(scData$qc.mito))
  expect_false(any(grepl("cannot flag", logs, fixed = TRUE)))
})
