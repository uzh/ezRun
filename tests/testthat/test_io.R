context("Tests functions in io.R")

test_that("Tests ezIntString() and mioString()", {
  num = 4.5e7
  intString = ezIntString(num)
  expect_is(intString, "character")
  expect_identical(as.numeric(intString), num)
  mioString = mioString(num)
  expect_is(mioString, "character")
  expect_identical(as.numeric(mioString) * 1e6, num)
})

test_that("Tests ezWrite.table() and ezRead.table()", {
  m1 = matrix(seq(0.01, 1, 0.01), 10)
  file = "exampleTable"
  colnames(m1) = letters[1:10]
  rownames(m1) = 1:10
  ezWrite.table(m1, file, sep = "\t")
  expect_true(file.exists(file))
  m2 = ezRead.table(file)
  expect_is(m2, "data.frame")
  expect_equal(m1, as.matrix(m2))
  ezSystem(paste("rm -fr", file))
})

test_that("Tests ezValidFilename()", {
  file1 = ezValidFilename("example:filename")
  file2 = ezValidFilename("or?to/remove(them", replace = "")
  expect_false(any(grepl("[$%#!?/:;() '=]", c(file1, file2))))
})

test_that("Tests removeSuffix() and getSuffix()", {
  file = "example.file"
  prefix = "example"
  names(prefix) = "example.file"
  expect_identical(getSuffix(file), "file")
  expect_equal(removeSuffix("example.file"), prefix)
})

test_that("Tests ezIsAbsolutePath()", {
  is = ezIsAbsolutePath("/ewukhlwa/gawrg/g/awg.txt")
  isnt = ezIsAbsolutePath("s/ega/aweh/awe.html")
  expect_true(is)
  expect_false(isnt)
})
