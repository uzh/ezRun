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
  expect_equal(job$start[2], elapsed$start[2])
  expect_less_than(job$start[3], elapsed$start[3])
})

test_that("Tests logMessage()", {
  start = logMessage("a method", param, "Starting")
  expect_null(start)
})
