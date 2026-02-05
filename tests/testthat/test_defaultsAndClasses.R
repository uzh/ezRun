context("Tests defaults and classes: 00defaults.R; 01classes.R; 02references.R")

test_that("Tests ezParam()", {
  globalParams = ezParam()
  expect_is(globalParams, "list")
  modifiedParam = ezParam(list(ram = 7))
  expect_false(globalParams$ram == modifiedParam$ram)
})

test_that("Tests ezIsSpecified()", {
  is = ezIsSpecified("a")
  expect_is(is, "logical")
  isnt = ezIsSpecified("")
  expect_false(isnt)
  alsoisnt = ezIsSpecified(NULL)
  expect_false(alsoisnt)
  isnteither = ezIsSpecified(NA)
  expect_false(isnteither)
})

test_that("Tests ezFrame()", {
  frame = ezFrame(
    first = 1:3,
    second = 5,
    "with space" = "text",
    row.names = letters[1:3]
  )
  expect_is(frame, "data.frame")
})

test_that("Tests EzDataset", {
  param = list(dataRoot = system.file(package = "ezRun", mustWork = TRUE))
  ds = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  path = ds$getFullPaths("Read1")
  expect_is(ds, "EzDataset")
  expect_equal(ds$getNames(), c("wt_1", "wt_2", "mut_1", "mut_2"))
  expect_is(ds$file, "character")
  expect_equal(ds$getLength(), 4)
  expect_equal(
    ds$getColumn("Read1"),
    c(
      wt_1 = "extdata/yeast_10k/wt_1_R1.fastq.gz",
      wt_2 = "extdata/yeast_10k/wt_2_R1.fastq.gz",
      mut_1 = "extdata/yeast_10k/mut_1_R1.fastq.gz",
      mut_2 = "extdata/yeast_10k/mut_2_R1.fastq.gz"
    )
  )
  expect_is(ds$meta, "data.frame")
  expect_is(ds$columnHasTag("File"), "logical")
  expect_true(ds$columnHasTag("File")[2])
  expect_is(path, "character")
  expect_true(all(grepl("^/", path)))
  expect_equal(ds$subset("wt_1")$getLength(), 1)
})

test_that("Tests copying an EzDataset and the method setColumn", {
  ds = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    )
  )
  expect_is(ds, "EzDataset")
  ds2 = ds$copy()
  expect_is(ds2, "EzDataset")
  expect_equal(ds, ds2)
  ds$setColumn("Read1", 1)
  expect_is(all.equal(ds, ds2, check.attributes = F), "character")
  expect_is(ds$meta[1, "Read1 [File]"], "numeric")
})

test_that("Tests EzApp", {
  ds = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    )
  )
  app = EzApp$new(runMethod = function(input, output, param) {}, name = "Test")
  expect_is(app, "EzApp")
  app2 = app$copy()
  expect_equal(app, app2)
  expect_is(app$runMethod, "function")
  expect_is(app$name, "character")
  expect_is(app$appDefaults, "data.frame")
  expect_null(app$run(
    input = ds,
    output = ds,
    param = list(process_mode = "DATASET")
  ))
})

test_that("Tests EzRef constructor", {
  refBuild = "Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18"
  ref = EzRef(param = ezParam(list(refBuild = refBuild)))
  expect_is(ref, "EzRef")
  expect_identical(ref["refBuild"], refBuild)
  expect_is(getOrganism(ref), "character")
})
