context("Tests the functions in ngsTools.r and rangesTools.r")

test_that("Tests strandValue() and strandName()", {
  values = c(plus = "+", minus = "-")
  names = c("+" = "plus", "-" = "minus")
  expect_identical(strandValue(1), values[1])
  expect_identical(strandValue(2), values[2])
  expect_identical(strandName(1), names[1])
  expect_identical(strandName(2), names[2])
})

test_that("Tests fixStrand() and flipStrand()", {
  strandValues = c("+", "+", "-", "*")
  expect_identical(fixStrand(strandValues), strandValues)
  expect_identical(
    fixStrand(strandValues, "antisense"),
    flipStrand(strandValues)
  )
  expect_is(fixStrand(strandValues, "both"), "Rle")
})

test_that("Tests getTuxedoLibraryType()", {
  expect_is(getTuxedoLibraryType("sense"), "character")
  expect_is(getTuxedoLibraryType("antisense"), "character")
  expect_is(getTuxedoLibraryType("both"), "character")
})

test_that("Tests isValidCigar()", {
  is = isValidCigar("3M5G")
  isnt = isValidCigar("3M5GA")
  expect_true(is)
  expect_false(isnt)
})

test_that("Tests shiftZeros()", {
  shifted = shiftZeros(runif(100) * 10, 2)
  expect_is(shifted, "numeric")
  expect_gt(min(shifted), 0.5)
})

test_that("Tests expandGRanges()", {
  IRange = IRanges(5000:5002, 8000:8002)
  GRange = GRanges(c("chr1", "chr1", "chr2"), IRange, strand = c("+", "-", "+"))
  expanded = expandGRanges(GRange)
  expect_is(expanded, "GRanges")
  for (i in 1:3) {
    expect_lt(expanded@ranges@start[i], GRange@ranges@start[i])
    expect_gt(expanded@ranges@width[i], GRange@ranges@width[i])
  }
})

test_that("Tests subSampleRle()", {
  rleobj = Rle(c(1, 1, 1, 2, 2, 3, 2, 2))
  subSample = subSampleRle(rleobj, 1:6)
  expect_is(subSample, "numeric")
})
