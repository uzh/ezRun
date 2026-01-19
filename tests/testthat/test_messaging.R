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



test_that("Tests ezLog()", {
  expect_error(ezLog("info message"), NA)
  expect_error(ezLog("debug message", level="debug"), NA)
  expect_error(ezLog("warn message", level="warn"), NA)
  expect_error(ezLog("error message", level="error"), NA)
  # expect_error(ezLog("fatal message", level="fatal"), NA) # fatal might stop execution or behave differently
  expect_error(ezLog("invalid level", level="invalid_level"), "Invalid log level: invalid_level")
})
