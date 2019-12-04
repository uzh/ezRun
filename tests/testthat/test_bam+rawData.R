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


test_that("Tests countReadsInFastq() from ngsio.r", {
  file2 = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
  input2 = EzDataset$new(file=file2, dataRoot=param$dataRoot)
  fqFiles = input2$getFullPaths("Read1")
  names(fqFiles) = NULL
  counted = countReadsInFastq(fqFiles)
  expect_is(counted, "integer")
  expect_identical(length(counted), nrow(input2$meta))
  expect_identical(names(counted), fqFiles)
})

