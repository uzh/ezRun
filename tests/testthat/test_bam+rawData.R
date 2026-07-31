context("Tests the functions in bamio.r; ngsio.r and rawData.R")

bamFile = system.file(
  "extdata",
  "ex1.bam",
  package = "Rsamtools",
  mustWork = TRUE
)
param = ezParam()
param$dataRoot = system.file(package = "ezRun", mustWork = TRUE)
file = system.file(
  "extdata/yeast_10k_STAR_counts/dataset.tsv",
  package = "ezRun",
  mustWork = TRUE
)
input = EzDataset$new(file = file, dataRoot = param$dataRoot)
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
  expect_gt(length(gapped), length(paired) * 1.8)
})

test_that("Tests getBamMultiMatching() from bamio.r", {
  multi = getBamMultiMatching(param, bamFile, nReads = 10000)
  expect_is(multi, "integer")
  expect_identical(names(multi), as.character(0:(length(multi) - 1)))
})

test_that("bamHasMarkedDuplicates() reads the marking off the bam header", {
  starBam = system.file(
    "extdata/yeast_10k_STAR/wt_1.bam",
    package = "ezRun",
    mustWork = TRUE
  )
  ## this fixture was aligned with markDuplicates=TRUE, its header carries the program
  ## record of picard and 88 of its reads have the duplicate flag set
  expect_true(bamHasMarkedDuplicates(starBam))
  ## ex1.bam was never marked
  expect_false(bamHasMarkedDuplicates(bamFile))
  expect_false(bamHasMarkedDuplicates("does_not_exist.bam"))
  expect_false(bamHasMarkedDuplicates(NA))

  skip_if(Sys.which("samtools") == "", "samtools not available")
  ## the decision has to come from the program record and from nothing else: dropping it
  ## from the header of the very same file must turn the answer around. A false positive
  ## would make getDupRateFromBam() skip picard and report the duplication rates of an
  ## unmarked bam.
  header = system2("samtools", c("view", "-H", starBam), stdout = TRUE)
  headerFile = tempfile(fileext = ".sam")
  unmarkedBam = tempfile(fileext = ".bam")
  on.exit(file.remove(headerFile, unmarkedBam))
  ## every record has to go, not just the one picard wrote: the ones added afterwards link
  ## back to it with PP:MarkDuplicates and that counts as marked as well
  writeLines(
    grep(
      "MarkDuplicates|markdup",
      header,
      invert = TRUE,
      value = TRUE,
      ignore.case = TRUE
    ),
    headerFile
  )
  system2(
    "samtools",
    c("reheader", shQuote(headerFile), shQuote(starBam)),
    stdout = unmarkedBam
  )
  expect_false(bamHasMarkedDuplicates(unmarkedBam))
})

test_that("Tests ezBamFlagKeepOnly() and ezBamFlagSkip() from bamio.r", {
  ## the three states of a scanBamFlag() argument: NA does not filter, TRUE keeps only the
  ## records with the flag, FALSE only those without it. A logical parameter must never be
  ## passed straight through, its FALSE would select what it is meant to leave out.
  expect_identical(ezBamFlagKeepOnly(TRUE), TRUE)
  expect_identical(ezBamFlagKeepOnly(FALSE), NA)
  expect_identical(ezBamFlagKeepOnly(NA), NA)
  expect_identical(ezBamFlagKeepOnly(NULL), NA)

  expect_identical(ezBamFlagSkip(TRUE), FALSE)
  expect_identical(ezBamFlagSkip(FALSE), NA)
  expect_identical(ezBamFlagSkip(NA), NA)
  expect_identical(ezBamFlagSkip(NULL), NA)

  ## NA must let both kinds of record through, which is what makes it the "off" state
  noFilter = scanBamFlag(isDuplicate = NA)
  expect_true(bitwAnd(noFilter[[1]], 0x400) != 0)
  expect_true(bitwAnd(noFilter[[2]], 0x400) != 0)
})

test_that("getBamMultiMatching() uses the NH tag and falls back without it", {
  starBam = system.file(
    "extdata/yeast_10k_STAR/wt_1.bam",
    package = "ezRun",
    mustWork = TRUE
  )
  ## the aligned bam files carry NH, so both implementations must agree. param$paired has
  ## to describe the library, both implementations count fragments only if it does.
  pairedParam = ezParam(list(paired = TRUE))
  fromNH = getBamMultiMatchingFromNH(pairedParam, starBam)
  fromQnames = getBamMultiMatchingFromQnames(pairedParam, starBam)
  expect_false(is.null(fromNH))
  expect_identical(sum(fromNH), sum(fromQnames))
  expect_identical(fromNH[["1"]], fromQnames[["1"]])

  ## ex1.bam has no NH tag, so the read name based counting has to take over
  expect_null(getBamMultiMatchingFromNH(param, bamFile))
  expect_identical(
    getBamMultiMatching(param, bamFile),
    getBamMultiMatchingFromQnames(param, bamFile)
  )
})

test_that("getBamMultiMatching() does not drop the proper pairs when they are not required", {
  starBam = system.file(
    "extdata/yeast_10k_STAR/wt_1.bam",
    package = "ezRun",
    mustWork = TRUE
  )
  pairedParam = ezParam(list(paired = TRUE))
  pairedParam$keepProperPairsOnly = TRUE
  onlyProper = getBamMultiMatching(pairedParam, starBam)
  ## keepProperPairsOnly=FALSE must switch the filter off rather than invert it, i.e. it
  ## can only ever add fragments to the count
  pairedParam$keepProperPairsOnly = FALSE
  anyPair = getBamMultiMatching(pairedParam, starBam)
  expect_gte(sum(anyPair), sum(onlyProper))
})


test_that("Tests countReadsInFastq() from ngsio.r", {
  file2 = system.file(
    "extdata/yeast_10k/dataset.tsv",
    package = "ezRun",
    mustWork = TRUE
  )
  input2 = EzDataset$new(file = file2, dataRoot = param$dataRoot)
  fqFiles = input2$getFullPaths("Read1")
  names(fqFiles) = NULL
  counted = countReadsInFastq(fqFiles)
  expect_is(counted, "integer")
  expect_identical(length(counted), nrow(input2$meta))
  expect_identical(names(counted), fqFiles)
})
