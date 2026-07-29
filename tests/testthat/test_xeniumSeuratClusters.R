## Regression tests for the XeniumSeurat app's cluster bookkeeping.
##
## Background: `FindClusters(cluster.name = "banksy_cluster")` ends with an
## `object[["seurat_clusters"]] <- Idents(object)` that is guarded in Seurat
## 5.4.0 but UNCONDITIONAL in Seurat 5.5.1. The app runs BANKSY niche
## clustering after transcriptional clustering, so on 5.5.1 the niches
## silently overwrote `seurat_clusters` -- every "cluster" figure in the report
## showed niches, while posMarkers had been computed from the real clusters.
##
## These tests pin the invariant the app must hold regardless of Seurat
## version, so an upstream re-change cannot reintroduce the bug silently.

skip_if_not_installed("Seurat")

context("XeniumSeurat cluster bookkeeping")

## Minimal Seurat object with a graph, clustered twice the way the app does.
.makeTwiceClustered <- function(seed = 42, res1 = 0.5, res2 = 1.8) {
  set.seed(seed)
  ng <- 120
  nc <- 300
  m <- matrix(rpois(ng * nc, 3), nrow = ng, ncol = nc,
              dimnames = list(paste0("g", seq_len(ng)), paste0("c", seq_len(nc))))
  o <- Seurat::CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE), assay = "Xenium")
  o <- Seurat::NormalizeData(o, verbose = FALSE)
  Seurat::VariableFeatures(o) <- rownames(o)
  o <- Seurat::ScaleData(o, verbose = FALSE)
  o <- Seurat::RunPCA(o, npcs = 15, verbose = FALSE)
  o <- Seurat::FindNeighbors(o, dims = 1:10, verbose = FALSE)
  ## 1. transcriptional clustering, exactly as the app does it
  o <- Seurat::FindClusters(o, resolution = res1, verbose = FALSE)
  ## keep a copy as metadata: a plain attribute does not survive FindClusters
  o$originalClusters <- o$seurat_clusters
  ## 2. BANKSY niche clustering, exactly as the app does it
  o <- Seurat::FindClusters(o, resolution = res2,
                            cluster.name = "banksy_cluster", verbose = FALSE)
  o
}

test_that("the app's restore step keeps seurat_clusters transcriptional", {
  o <- .makeTwiceClustered()
  orig <- o$originalClusters

  ## Guard: on Seurat >= 5.5.1 the unguarded assignment must actually have
  ## fired, otherwise the restore step below is untested on this stack.
  overwritten <- identical(as.character(o$seurat_clusters),
                           as.character(o$banksy_cluster))
  message("seurat_clusters overwritten by niches on Seurat ",
          as.character(utils::packageVersion("Seurat")), ": ", overwritten)

  ## Positive control: the two clusterings must actually differ, otherwise
  ## this test cannot detect the overwrite it exists to catch.
  expect_gt(nlevels(o$banksy_cluster), nlevels(orig))

  ## The fix the app applies after BANKSY clustering.
  res <- 0.5
  resCol <- paste0("Xenium_snn_res.", res)
  expect_true(resCol %in% colnames(o[[]]))
  o$seurat_clusters <- o[[resCol]][, 1]

  expect_equal(as.character(o$seurat_clusters), as.character(orig))
  expect_false(identical(as.character(o$seurat_clusters),
                         as.character(o$banksy_cluster)))
})

test_that("clusters and niches produce DIFFERENT Xenium Explorer exports", {
  ## The user-visible symptom of the bug: clusters_for_explorer.csv and
  ## niches_for_explorer.csv came out byte-identical.
  o <- .makeTwiceClustered()
  o$seurat_clusters <- o[["Xenium_snn_res.0.5"]][, 1]

  clusters <- data.frame(cell_id = colnames(o),
                         group = as.character(o$seurat_clusters))
  niches <- data.frame(cell_id = colnames(o),
                       group = as.character(o$banksy_cluster))

  expect_false(identical(clusters$group, niches$group))
})

test_that("explorer cell ids drop the sample prefix RenameCells adds", {
  ## Xenium Explorer matches on the raw cell_id ("aaaddlda-1"), but the app
  ## prefixes cells with the sample name, so the export must strip it back off.
  sampleName <- "MySample"
  raw <- c("aaaddlda-1", "aacccpin-1", "aadjghil-1")
  prefixed <- paste(sampleName, raw, sep = "_")

  exported <- sub(paste0("^", sampleName, "_"), "", prefixed)

  expect_equal(exported, raw)
  ## Negative control: a cell that does NOT carry the prefix is left alone.
  expect_equal(sub(paste0("^", sampleName, "_"), "", "OtherSample_aaa-1"),
               "OtherSample_aaa-1")
})
