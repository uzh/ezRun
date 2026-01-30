context("Test handling of DE results: twoGroups.R")


test_that("test instantiation, saving and loading", {
  deResult = EzResult$new()
  deResult = EzResult$new(
    param = list(p1 = "p1", p2 = "p2"),
    rawData = list(a = 0, b = 1),
    result = list(u = 0, v = 1)
  )
  rdFile = tempfile(fileext = ".RData")
  deResult$saveToFile(file = rdFile)

  deResult2 = EzResult$new(file = rdFile)
  expect_true(deResult2$param$p1 == deResult$param$p1)
  deResult$param$p1 = "p1Mod"
  expect_false(deResult2$param$p1 == deResult$param$p1)
  deResult2 = EzResult$new(file = rdFile)
  expect_false(deResult2$param$p1 == deResult$param$p1)
})
