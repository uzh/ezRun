context("Test cell-type annotation column selection: multiOmicsUtils.R")

## pickCellTypeColumn() only ever touches obj@meta.data, so a two-slot stand-in
## exercises the same code path without adding Seurat to Suggests (it is not
## declared there, so a real fixture would silently skip in CI).
setClass("FakeSeuratForPick", representation(meta.data = "data.frame"))

fakeObj <- function(...) {
  md <- as.data.frame(list(...), stringsAsFactors = FALSE, check.names = FALSE)
  new("FakeSeuratForPick", meta.data = md)
}

## The eight columns AzimuthAPI::CloudAzimuth lands on scData. Verified against
## delivered objects under /srv/gstore/projects/p31662/oNA_ScSeurat_2026-02-05
## and /srv/gstore/projects/p36614/o5495_o5444_ScSeurat_2026-03-12, and against
## a live CloudAzimuth call. Identical in AzimuthAPI 0.1.0 and 1.0.0.
panHumanMeta <- function() {
  fakeObj(
    nCount_RNA                = 1000,
    seurat_clusters           = "0",
    full_hierarchical_labels  = "immune|lymphoid|T cell",
    final_level_labels        = "CD4 T cell",
    final_level_confidence    = 0.91,
    full_consistent_hierarchy = TRUE,
    azimuth_broad             = "immune",
    azimuth_medium            = "T cell",
    azimuth_fine              = "CD4 T cell",
    azimuth_label             = "CD4 T cell"
  )
}

test_that("a Pan-Human-only object yields a Pan-Human label column, not NULL", {
  hit <- pickCellTypeColumn(panHumanMeta())
  expect_false(is.null(hit))
  expect_true(hit %in% c("azimuth_label", "azimuth_medium", "azimuth_broad",
                         "azimuth_fine", "final_level_labels"))
})

test_that("Pan-Human outranks scType", {
  obj <- panHumanMeta()
  obj@meta.data$sctype_classification <- "T cells"
  hit <- pickCellTypeColumn(obj)
  expect_false(identical(hit, "sctype_classification"))
  expect_identical(hit, "azimuth_label")
})

test_that("Pan-Human outranks the Azimuth tissue reference and SingleR", {
  obj <- panHumanMeta()
  obj@meta.data$Azimuth.celltype.l2 <- "CD4 Naive"
  obj@meta.data$MonacoImmuneData_cluster <- "T cells"
  expect_identical(pickCellTypeColumn(obj), "azimuth_label")
})

test_that("non-label Pan-Human columns are never selected", {
  ## final_level_confidence is numeric, full_consistent_hierarchy is logical and
  ## full_hierarchical_labels is a pipe-delimited path with near-per-cell
  ## cardinality. CloudAzimuth writes all three, none works as a group.by.
  obj <- fakeObj(
    nCount_RNA                = 1000,
    final_level_confidence    = 0.91,
    full_consistent_hierarchy = TRUE,
    full_hierarchical_labels  = "immune|lymphoid|T cell"
  )
  expect_null(pickCellTypeColumn(obj))
})

test_that("cellxgene label transfer is found under its real column name", {
  obj <- fakeObj(
    nCount_RNA = 1000,
    predicted.id.cellxgene.authorlabel = "T cell",
    prediction.score.max = 0.8
  )
  expect_identical(pickCellTypeColumn(obj), "predicted.id.cellxgene.authorlabel")
})

test_that("cellxgene outranks Pan-Human (documented FGCZ priority)", {
  obj <- panHumanMeta()
  obj@meta.data$predicted.id.cellxgene.authorlabel <- "T cell"
  expect_identical(pickCellTypeColumn(obj), "predicted.id.cellxgene.authorlabel")
})

test_that("CyteTypeR still wins over everything", {
  obj <- panHumanMeta()
  obj@meta.data$predicted.id.cellxgene.authorlabel <- "T cell"
  obj@meta.data$sctype_classification <- "T cells"
  obj@meta.data$CyteTypeR_annotation <- "CD4 T"
  expect_identical(pickCellTypeColumn(obj), "CyteTypeR_annotation")
})

test_that("an object with no annotation column returns NULL", {
  expect_null(pickCellTypeColumn(fakeObj(nCount_RNA = 1000,
                                         nFeature_RNA = 500,
                                         seurat_clusters = "0")))
})

test_that("non-Pan-Human writers stay reachable (regression guards)", {
  ## Azimuth::RunAzimuth run manually upstream, per README.md.
  expect_identical(
    pickCellTypeColumn(fakeObj(predicted.celltype.l1 = "T",
                               predicted.celltype.l2 = "CD4 T")),
    "predicted.celltype.l2")
  ## ezRun's own manual cluster-labelling apps. intersect() is case-sensitive,
  ## so the generic "celltype"/"CellType" spellings do not cover these.
  expect_identical(pickCellTypeColumn(fakeObj(cellType = "CD4 T")), "cellType")
  expect_identical(pickCellTypeColumn(fakeObj(cellTypeIntegrated = "CD4 T")),
                   "cellTypeIntegrated")
  ## BD Rhapsody objects loaded by loadBDRhapsody(), which readRDS's BD's own
  ## *_Seurat.rds and bypasses attachUpstreamAnnotations entirely.
  expect_identical(pickCellTypeColumn(fakeObj(Cell_Type_Experimental = "T cell")),
                   "Cell_Type_Experimental")
})
