context("Tests the functions in ngsTools.r")

test_that("Tests the functions strandValue() and strandName()", {
  values = c(plus="+", minus="-")
  names = c("+"="plus", "-"="minus")
  expect_identical(strandValue(1), values[1])
  expect_identical(strandValue(2), values[2])
  expect_identical(strandName(1), names[1])
  expect_identical(strandName(2), names[2])
})


