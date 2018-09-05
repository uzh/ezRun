context("Tests the functions in bamio.r; ngsio.r and rawData.R")

bamFile = system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
param = ezParam()
param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
file = system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE)
input = EzDataset$new(file=file, dataRoot=param$dataRoot)
rawData = loadCountDataset(input, param)

test_that("Tests ezBamSeqNames() and ezBamSeqLengths() from bamio.r", {
  lengths = ezBamSeqLengths(bamFile)
  expect_is(lengths, "integer")
  names = ezBamSeqNames(bamFile)
  expect_is(names, "character")
})

test_that("Tests ezScanBam() from bamio.r", {
  scanned = ezScanBam(bamFile)
  expect_is(scanned, "list")
  expect_is(scanned$seq, "DNAStringSet")
  expect_is(scanned$qual, "PhredQuality")
})

test_that("Tests ezReadGappedAlignments() and ezReadPairedAlignments() from bamio.r", {
  gapped = ezReadGappedAlignments(bamFile)
  paired = ezReadPairedAlignments(bamFile)
  expect_is(gapped, "GAlignments")
  expect_is(paired, "GAlignments")
  expect_gt(length(gapped), length(paired)*1.8)
})

test_that("Tests getBamMultiMatching() from bamio.r", {
  multi = getBamMultiMatching(param, bamFile, nReads=10000)
  expect_is(multi, "integer")
  expect_identical(names(multi), as.character(0:(length(multi)-1)))
})

test_that("Tests loadCountDataset() from ngsio.r", {
  expect_is(rawData, "list")
  expect_identical(rawData$dataset, input$meta)
  expect_is(rawData$counts, "matrix")
  expect_is(rawData$isLog, "logical")
  expect_is(rawData$presentFlag, "matrix")
  expect_is(rawData$seqAnno, "data.frame")
  expect_is(rawData$featureLevel, "character")
  expect_is(rawData$type, "character")
  expect_is(rawData$countName, "character")
})

test_that("Tests countReadsInFastq() from ngsio.r", {
  file2 = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
  input2 = EzDataset$new(file=file2, dataRoot=param$dataRoot)
  fqFiles = input2$getFullPaths("Read1")
  names(fqFiles) = NULL
  counted = countReadsInFastq(fqFiles)
  expect_is(counted, "numeric")
  expect_identical(length(counted), nrow(input2$meta))
  expect_identical(names(counted), fqFiles)
})

test_that("Tests getSignal() and getLog2Signal() from rawData.R", {
  ## TODOP: get rawData with signal to test these.
})

test_that("Tests getRpkm() and getTpm() from rawData.R", {
  ## TODOP: get rawData with proper seqAnno to test these.
})

test_that("Tests selectFeatures() from rawData.R", {
  keep = "YFR014C"
  selected = selectFeatures(rawData, keep)
  expect_equal(length(selected$counts), 4)
})
