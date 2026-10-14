context("CellRangerMulti: Flex v2 chemistry from read lengths")

## Flex-v2-R1 needs the probe barcode at R1[40:50], Flex-v2-RNA-R2 at R2[79:89]
## (CellRanger 10.1 chemistry_defs.json). The old rule looked at R1 only
## (R1 > 28 -> Flex-v2-R1), which pinned an invalid chemistry for R1 29-49.

test_that("Read 1 config and full-length runs get Flex-v2-R1", {
  expect_equal(flexV2Chemistry(54, 50), "Flex-v2-R1")
  expect_equal(flexV2Chemistry(50, 50), "Flex-v2-R1")
  expect_equal(flexV2Chemistry(54, 90), "Flex-v2-R1")
  expect_equal(flexV2Chemistry(151, 151), "Flex-v2-R1")
  expect_equal(flexV2Chemistry(54, 50, hasAntibody = TRUE), "Flex-v2-R1")
})

test_that("Read 2 config gets Flex-v2-RNA-R2, or auto when a protein library is present", {
  expect_equal(flexV2Chemistry(28, 90), "Flex-v2-RNA-R2")
  expect_equal(flexV2Chemistry(49, 89), "Flex-v2-RNA-R2")
  expect_true(is.na(flexV2Chemistry(28, 90, hasAntibody = TRUE)))
})

test_that("reads too short for either barcode position stop with the lengths", {
  expect_error(flexV2Chemistry(40, 50), "R1 = 40 nt and R2 = 50 nt")
  expect_error(flexV2Chemistry(28, 88), "R2 >= 89")
  expect_error(flexV2Chemistry(49, 50), "R1 >= 50")
})

test_that("maxFastqReadLength takes the longest of the first reads", {
  fq <- tempfile(fileext = ".fastq.gz")
  con <- gzfile(fq, "w")
  writeLines(c("@a", strrep("A", 28), "+", strrep("F", 28),
               "@b", strrep("C", 54), "+", strrep("F", 54)), con)
  close(con)
  expect_equal(maxFastqReadLength(fq), 54)
  expect_equal(maxFastqReadLength(fq, n = 1), 28)
})
