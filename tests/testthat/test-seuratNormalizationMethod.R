## Tests for the normalizationMethod knob (SCTransform | LogNormalize) added
## 2026-09-03 for ScSeuratCombine.
##
## The failure this guards against is silent: the log-normalization path has to
## produce scale.data itself (SCTransform makes its own), or RunPCA and
## DoHeatmap get nothing to work with. It also has to leave the clustering
## columns named RNA_snn_res.<r> rather than SCT_snn_res.<r>, which is what the
## report's clustree prefix is derived from.

suppressPackageStartupMessages(library(Seurat))

makeSampleList <- function(nCells = 150L, nGenes = 400L) {
  set.seed(38)
  makeOne <- function(tag) {
    counts <- matrix(
      rnbinom(nGenes * nCells, size = 2, mu = 5),
      nrow = nGenes,
      dimnames = list(
        paste0("Gene", seq_len(nGenes)),
        paste0(tag, "_Cell", seq_len(nCells))
      )
    )
    obj <- CreateSeuratObject(as(counts, "dgCMatrix"))
    obj$Sample <- tag
    obj$Condition <- tag
    obj
  }
  list(A = makeOne("A"), B = makeOne("B"))
}

testParam <- function(normalizationMethod) {
  list(
    name = "ScSeuratCombine",
    nfeatures = 200,
    npcs = 10,
    resolution = 0.6,
    normalizationMethod = normalizationMethod,
    SCT.regress.CellCycle = FALSE
  )
}

## Exactly the expression inst/templates/ScSeuratCombine.Rmd uses to pick the
## clustree prefix. If this stops finding a column the report stops rendering.
rmdClustreePrefix <- function(scData) {
  resCols <- grep("_snn_res\\.", colnames(scData@meta.data), value = TRUE)
  if (length(resCols) == 0) {
    return(NA_character_)
  }
  paste0(sub("_snn_res\\..*$", "", resCols[1]), "_snn_res.")
}

scaleDataCols <- function(scData, assay) {
  layer <- SeuratObject::LayerData(scData, assay = assay, layer = "scale.data")
  if (is.null(layer)) 0L else ncol(layer)
}

test_that("the knob defaults to SCTransform when the parameter is absent", {
  expect_equal(seuratNormalizationMethod(list()), "SCTransform")
  expect_equal(seuratNormalizationMethod(list(normalizationMethod = "")), "SCTransform")
  expect_equal(
    seuratNormalizationMethod(list(normalizationMethod = "LogNormalize")),
    "LogNormalize"
  )
  expect_equal(seuratAnalysisAssay(list()), "SCT")
  expect_equal(
    seuratAnalysisAssay(list(normalizationMethod = "LogNormalize")),
    "RNA"
  )
  expect_equal(
    seuratAnalysisAssay(
      list(normalizationMethod = "LogNormalize"),
      rawAssay = "Spatial"
    ),
    "Spatial"
  )
})

test_that("log-normalizing per sample equals log-normalizing the merged object", {
  ## The claim that makes JoinLayers lossless in seuratScaleMergedLogNorm: the
  ## transform is per cell, so the sample split cannot change a value.
  scDataList <- makeSampleList(nCells = 40L, nGenes = 60L)

  perSample <- suppressWarnings(
    seuratNormalizeSampleList(scDataList, testParam("LogNormalize"))
  )
  perSample <- merge(perSample[[1]], perSample[-1], merge.data = TRUE)
  perSample <- JoinLayers(perSample, assay = "RNA")

  merged <- merge(scDataList[[1]], scDataList[-1])
  merged <- JoinLayers(merged, assay = "RNA")
  merged <- NormalizeData(
    merged,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  a <- SeuratObject::LayerData(perSample, assay = "RNA", layer = "data")
  b <- SeuratObject::LayerData(merged, assay = "RNA", layer = "data")
  expect_equal(dim(a), dim(b))
  expect_equal(max(abs(a - b[rownames(a), colnames(a)])), 0)
})

test_that("LogNormalize clusters in the RNA assay and still produces scale.data", {
  scData <- suppressWarnings(
    cellClustNoCorrection(makeSampleList(), testParam("LogNormalize"))
  )

  ## no SCT assay at all, so PrepSCTFindMarkers must be skipped downstream
  expect_false("SCT" %in% Seurat::Assays(scData))
  expect_equal(DefaultAssay(scData), "RNA")

  ## the thing that silently breaks RunPCA and DoHeatmap if it is forgotten
  expect_gt(scaleDataCols(scData, "RNA"), 0L)
  expect_equal(scaleDataCols(scData, "RNA"), ncol(scData))

  ## layers joined, not left split per sample
  expect_equal(
    sum(grepl("^data", SeuratObject::Layers(scData[["RNA"]]))),
    1L
  )

  expect_true(all(c("pca", "umap", "tsne") %in% Reductions(scData)))
  expect_equal(rmdClustreePrefix(scData), "RNA_snn_res.")
  expect_gt(nlevels(Idents(scData)), 0L)
})

test_that("SCTransform stays the historical behaviour", {
  ## Positive control for the assertions above: the same checks distinguish the
  ## two paths, so they are measuring the normalization and not something else.
  scData <- suppressWarnings(
    cellClustNoCorrection(makeSampleList(), testParam("SCTransform"))
  )

  expect_true("SCT" %in% Seurat::Assays(scData))
  expect_equal(DefaultAssay(scData), "SCT")
  expect_gt(scaleDataCols(scData, "SCT"), 0L)
  expect_equal(rmdClustreePrefix(scData), "SCT_snn_res.")
})
