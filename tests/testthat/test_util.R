context("Tests the functions in util.r")


test_that("Tests lastVal", {
  x = 1:5
  expect_equal(lastVal(x), rev(x)[1])
})

test_that("Tests vennFromSets", {
  expect_error(vennFromSets(list(a=1)))
  expect_error(vennFromSets(list(a=1,b=2,c=3,d=4)))
  expect_error(vennFromSets(list(1)))
  expect_error(vennFromSets(list(1,2,3)))
})

test_that("Tests tableFromSets", {
  expect_error(tableFromSets(list(a=1)))
  expect_error(tableFromSets(list(a=1,b=2,c=3,d=4)))
  expect_error(tableFromSets(list(1)))
  expect_error(tableFromSets(list(1,2,3)))
})

