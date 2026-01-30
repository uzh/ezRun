context("Functions from datasets.r")

test_that("Tests makeMimimal...EndReadDataset()", {
  fqDir = system.file("extdata/yeast_10k", package = "ezRun", mustWork = TRUE)
  species = "Example"
  minSingleDs = makeMinimalSingleEndReadDataset(fqDir, species)
  expect_is(minSingleDs, "data.frame")
})

test_that("Tests ezDesignFromDataset() and similar functions", {
  file = system.file(
    "extdata/yeast_10k/dataset.tsv",
    package = "ezRun",
    mustWork = TRUE
  )
  ds = EzDataset$new(
    file = file,
    dataRoot = system.file(".", package = "ezRun", mustWork = TRUE)
  )
  design = ezDesignFromDataset(ds$meta)
  expect_is(design, "data.frame")
  cond1 = ezConditionsFromDesign(design)
  expect_is(cond1, "character")
  cond2 = ezConditionsFromDataset(ds$meta)
  expect_identical(cond1, cond2)
  replicates = addReplicate(apply(design, 1, paste, collapse = "_"))
  expect_is(replicates, "character")
  expect_true(all(grepl("_", replicates)))
})
