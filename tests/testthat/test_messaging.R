context("Functions from messaging.r")

param = ezParam()

test_that("Tests ezTime()", {
  time = ezTime()
  expect_is(time, "character")
  expect_equal(nchar(time), 20)
})

test_that("Tests ezJobStart() and ezWriteElapsed()", {
  job = ezJobStart("a job")
  elapsed = ezWriteElapsed(job)
  expect_is(job$name, "character")
  expect_is(job$start, "proc_time")
  expect_is(elapsed$name, "character")
  expect_is(elapsed$start, "proc_time")
})

test_that("Tests logMessage()", {
  start = logMessage("a method", param, "Starting")
  expect_null(start)
})
