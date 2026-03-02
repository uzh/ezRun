context("Tests the functions in util.r and util-genome.R")

m1 = matrix(1:20, 5)
rownames(m1) = 1:5
l1 = list(a = 1:3, b = c(2, 5), c = 4:8)
numbers = runif(100)

test_that("Tests lastVal()", {
  x = 1:5
  expect_equal(lastVal(x), rev(x)[1])
})

test_that("Tests vennFromSets()", {
  expect_error(vennFromSets(list(a = 1)))
  expect_error(vennFromSets(list(a = 1, b = 2, c = 3, d = 4)))
  expect_error(vennFromSets(list(1)))
  expect_error(vennFromSets(list(1, 2, 3)))
})

test_that("Tests tableFromSets()", {
  expect_error(tableFromSets(list(a = 1)))
  expect_error(tableFromSets(list(a = 1, b = 2, c = 3, d = 4)))
  expect_error(tableFromSets(list(1)))
  expect_error(tableFromSets(list(1, 2, 3)))
})

test_that("Tests shrinkToRange()", {
  low = runif(1, min = 0.4, max = 0.6)
  high = runif(1, min = 0.6)
  shrunkRange = shrinkToRange(runif(100), c(low, high))
  expect_false(all(shrunkRange < low))
  expect_false(all(shrunkRange > high))
  expect_gt(length(unique(shrunkRange)), 1)
})

test_that("Tests interleaveMatricesByColumn()", {
  matrix1 = matrix(1:10, 2)
  matrix2 = matrix(11:20, 2)
  colnames(matrix1) = as.character(1:5)
  rownames(matrix1) = c("a", "b")
  colnames(matrix2) = as.character(1:5)
  mNew = interleaveMatricesByColumn(matrix1, matrix2)
  expect_is(mNew, "matrix")
  expect_identical(
    length(as.vector(mNew)),
    length(as.vector(matrix1)) + length(as.vector(matrix2))
  )
})

test_that("Tests ezNorm()", {
  numbers = matrix(runif(1000), 50)
  none = ezNorm(numbers)
  expect_identical(numbers, none)
  quantile = ezNorm(numbers, method = "quantile")
  expect_lt(sd(quantile), sd(none))
  logMean = ezNorm(numbers, method = "logMean")
  expect_gt(max(logMean), max(none))
  median = ezNorm(numbers, method = "median")
  expect_gt(max(median), max(none))
  vsn = ezNorm(numbers, method = "vsn")
  expect_gt(mean(vsn), max(c(none, quantile, logMean, median)))
})

test_that("Tests ezCut()", {
  cut = ezCut(1:10, breaks = c(2, 5, 7), prefix = letters[1:4])
  expect_is(cut, "factor")
})

test_that("Tests isError()", {
  isnt = isError("error")
  is = isError(list(a = 3:5, error = 3))
  expect_false(isnt)
  expect_true(is)
})

test_that("Tests ezGrepl()", {
  grepl1 = ezGrepl(c(2, 4), 1:100)
  grepl2 = ezGrepl(c(2, 4), 1:100, combine = "and")
  expect_is(grepl1, "logical")
  expect_gt(length(which(grepl1)), length(which(grepl2)))
})

test_that("Tests ezSplit()", {
  split1 = ezSplit(letters[1:5], "b")
  split2 = ezSplit(rep("abcde", 4), letters[1:4])
  expect_is(split1, "matrix")
})

test_that("Tests trimWhiteSpace()", {
  toTrim = "    bla    "
  trimmed = trimWhiteSpace(toTrim)
  expect_lt(nchar(trimmed), nchar(toTrim))
})

test_that("Tests hasFilesafeCharacters()", {
  expect_true(hasFilesafeCharacters("a"))
  expect_false(hasFilesafeCharacters("a\n"))
  expect_true(all(hasFilesafeCharacters(c("1", "2"))))
  expect_false(all(hasFilesafeCharacters(list("1", "2 x"))))
})

test_that("Tests ezMatrix()", {
  em1 = ezMatrix(1, rows = 1:4, cols = 1:3)
  em2 = ezMatrix(3:6, dim = c(4, 6))
  expect_is(em1, "matrix")
  expect_is(em2, "matrix")
})

test_that("Tests ezScaleColumns()", {
  scaledMatrix = ezScaleColumns(m1, 1:4)
  expect_identical(scaledMatrix[, 1], m1[, 1])
  expect_equal(scaledMatrix[, 2], m1[, 2] * 2)
})

test_that("Tests ezGeomean()", {
  one = ezGeomean(numbers)
  two = exp(mean(log(numbers)))
  expect_identical(one, two)
})


test_that("Tests inverseMapping() and makeMultiMapping()", {
  invMapped = inverseMapping(l1)
  expect_is(invMapped, "list")
  multiMapped = makeMultiMapping(l1)
  expect_is(multiMapped, "data.frame")
})

test_that("Tests ezMclapply()", {
  applied = ezMclapply(l1, sum)
  expect_is(applied, "list")
})

test_that("Tests ezDuplicated() and ezMultiplicated()", {
  v1 = c(1, 2, 3, 4, 5, 4, 3, 2, 3, 4, 5, 6, 7, 8, 7)
  isDup1 = ezDuplicated(v1)
  isDup2 = ezDuplicated(v1, "all")
  expect_is(isDup1, "logical")
  expect_gt(length(which(isDup2)), length(which(isDup1)))
  isMul1 = ezMultiplicated(v1)
  isMul2 = ezMultiplicated(v1, 3)
  isMul3 = ezMultiplicated(v1, 2, "all")
  expect_is(isMul1, "logical")
  expect_gt(length(which(isMul3)), length(which(isMul2)))
})

test_that("Tests ezReplicateNumber()", {
  x = c("a", "c", "a", "b", "c")
  repli = ezReplicateNumber(x)
  expect_is(repli, "integer")
})

test_that("Tests ezCollapse()", {
  list1 = list(a = c(1, "", 6), c = c("rsrg", "yjrt", NA, 6))
  collapse1 = ezCollapse(list1, sep = "_")
  collapse2 = ezCollapse(list1, na.rm = T, empty.rm = T, uniqueOnly = T)
  expect_is(collapse1, "character")
  expect_equal(length(collapse2), 1)
})

test_that("Tests ezSplitLongLabels()", {
  nSplit = 22
  a = paste(letters[1:nSplit], collapse = "")
  b = paste(letters[1:(nSplit + 1)], collapse = "")
  charVec = c(a, b)
  splittedLabels = ezSplitLongLabels(charVec, nSplit)
  expect_false(grepl("\n", splittedLabels[1]))
  expect_true(grepl("\n", splittedLabels[2]))
})

test_that("Tests functions in util-genome.R", {
  coord = makeCoordinate("chrm", 3, 150, 300)
  expect_is(coord, "character")
  splitted = splitCoordinate(coord)
  expect_is(splitted, "data.frame")
  region = splitRegion(coord)
  expect_is(region, "list")
})

test_that("Tests ezFindRobj()", {
  # Create temp directory with test files

tmpDir <- tempfile()
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  # Create test files
  testData <- list(a = 1, b = 2)
  rdsPath <- file.path(tmpDir, "test.rds")
  saveRDS(testData, rdsPath)

  # Test with full path including extension
  expect_equal(ezFindRobj(rdsPath), rdsPath)

  # Test with base path (no extension)
  basePath <- file.path(tmpDir, "test")
  expect_equal(ezFindRobj(basePath), rdsPath)

  # Test non-existent file returns NULL
  expect_null(ezFindRobj(file.path(tmpDir, "nonexistent")))
  expect_null(ezFindRobj(file.path(tmpDir, "nonexistent.rds")))

  # Test priority: .qs2 > .qs > .rds
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2Path <- file.path(tmpDir, "test.qs2")
    qs2::qs_save(testData, qs2Path)
    expect_equal(ezFindRobj(basePath), qs2Path)
  }
})

test_that("Tests ezRobjExists()", {
  # Create temp directory with test file
  tmpDir <- tempfile()
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  testData <- list(a = 1, b = 2)
  rdsPath <- file.path(tmpDir, "exists.rds")
  saveRDS(testData, rdsPath)

  # Test existing file
  expect_true(ezRobjExists(rdsPath))
  expect_true(ezRobjExists(file.path(tmpDir, "exists")))

  # Test non-existent file
  expect_false(ezRobjExists(file.path(tmpDir, "nonexistent")))
  expect_false(ezRobjExists(file.path(tmpDir, "nonexistent.rds")))
})

test_that("Tests ezLoadRobj()", {
  # Create temp directory with test files
  tmpDir <- tempfile()
  dir.create(tmpDir)
  on.exit(unlink(tmpDir, recursive = TRUE))

  testData <- list(a = 1:10, b = "test")

  # Test loading .rds file
  rdsPath <- file.path(tmpDir, "data.rds")
  saveRDS(testData, rdsPath)
  loaded <- ezLoadRobj(rdsPath)
  expect_equal(loaded, testData)

  # Test auto-detection with base path
  basePath <- file.path(tmpDir, "data")
  loaded2 <- ezLoadRobj(basePath)
  expect_equal(loaded2, testData)

  # Test error on non-existent file
  expect_error(ezLoadRobj(file.path(tmpDir, "nonexistent")))
  expect_error(ezLoadRobj(file.path(tmpDir, "nonexistent.rds")))

  # Test qs format if available
  if (requireNamespace("qs", quietly = TRUE)) {
    qsPath <- file.path(tmpDir, "qsdata.qs")
    qs::qsave(testData, qsPath)
    loadedQs <- ezLoadRobj(qsPath)
    expect_equal(loadedQs, testData)
  }

  # Test qs2 format if available
  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2Path <- file.path(tmpDir, "qs2data.qs2")
    qs2::qs_save(testData, qs2Path)
    loadedQs2 <- ezLoadRobj(qs2Path)
    expect_equal(loadedQs2, testData)
  }
})

test_that("Tests .ezDefaultThreads()", {
  threads <- ezRun:::.ezDefaultThreads()
  expect_is(threads, "integer")
  expect_gte(threads, 1L)
  expect_lte(threads, 4L)
})
