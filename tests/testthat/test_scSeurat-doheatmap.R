context("ScSeurat: DoHeatmap must not kill the report on empty features")

## Regression test for the ScSeurat report crash seen in production
## (p41818/o41876 2026-06-16, p31662/o36310 2026-07-12):
##
##   Error in `names(x) <- value`:
##     'names' attribute [3] must be the same length as the vector [1]
##   1. Seurat::DoHeatmap(scData, features = unique(genesToPlot))
##   2.   Seurat::SingleRasterMap(...)
##   3.     base::`colnames<-`(...)
##
## Cause: DoHeatmap only errors when *some* requested features are missing --
## `any(!features %in% possible.features)` is FALSE for character(0) -- so an
## empty feature set sails through, SingleRasterMap gets a 0-feature matrix,
## Melt() collapses it to a single column, and colnames<- (which assigns 3
## names) blows up. genesToPlot is empty whenever no markers were found, e.g. a
## low-cell sample that collapses to a single cluster.
##
## ScSeurat.Rmd wraps DoHeatmap in tryCatch so the heatmap is skipped and the
## rest of the report still renders.

test_that("empty feature set makes DoHeatmap fatal (the bug)", {
  skip_if_not_installed("Seurat")
  suppressMessages(library(Seurat))

  set.seed(1)
  counts <- matrix(
    rpois(200, 5), nrow = 20,
    dimnames = list(paste0("g", 1:20), paste0("c", 1:10))
  )
  scData <- CreateSeuratObject(counts)
  scData <- NormalizeData(scData, verbose = FALSE)
  scData <- suppressWarnings(ScaleData(scData, verbose = FALSE))
  # low-cell reality: everything lands in one cluster -> no markers found
  Idents(scData) <- factor(rep("0", 10))
  genesToPlot <- character(0)

  expect_error(
    DoHeatmap(scData, features = unique(genesToPlot)),
    "must be the same length as the vector"
  )
})

test_that("the tryCatch guard skips the heatmap and the report continues", {
  skip_if_not_installed("Seurat")
  suppressMessages(library(Seurat))

  set.seed(1)
  counts <- matrix(
    rpois(200, 5), nrow = 20,
    dimnames = list(paste0("g", 1:20), paste0("c", 1:10))
  )
  scData <- CreateSeuratObject(counts)
  scData <- NormalizeData(scData, verbose = FALSE)
  scData <- suppressWarnings(ScaleData(scData, verbose = FALSE))
  Idents(scData) <- factor(rep("0", 10))
  genesToPlot <- character(0)

  # mirrors the guard in inst/templates/ScSeurat.Rmd
  caught <- FALSE
  suppressMessages(
    tryCatch(
      print(DoHeatmap(scData, features = unique(genesToPlot))),
      error = function(e) caught <<- TRUE
    )
  )
  expect_true(caught)

  # the rest of the report still renders
  expect_s3_class(DotPlot(scData, features = "g1"), "ggplot")
})
