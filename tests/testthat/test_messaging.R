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
  capture_log <- function(code) {
    paste(capture.output(code, type = "message"), collapse = "\n")
  }

  expect_match(
    capture_log(ezLog("info message")),
    "INFO \\(ezRun\\) \\[.*\\] info message"
  )

  ezLogChangeLevel("debug")

  expect_match(
    capture_log(ezLog("debug message", level = "debug")),
    "DEBUG \\(ezRun\\) \\[.*\\] debug message"
  )

  expect_match(
    capture_log(ezLog("warn message", level = "warn")),
    "WARN \\(ezRun\\) \\[.*\\] warn message"
  )
  expect_match(
    capture_log(ezLog("error message", level = "error")),
    "ERROR \\(ezRun\\) \\[.*\\] error message"
  )
  expect_error(
    ezLog("invalid level", level = "invalid_level"),
    "Invalid log level: invalid_level"
  )
  ezLogChangeLevel("warn")
  expect_equal(capture_log(ezLog("info hidden", level = "info")), "")
  expect_match(
    capture_log(ezLog("warn visible", level = "warn")),
    "WARN \\(ezRun\\) \\[.*\\] warn visible"
  )
  ezLogChangeLevel("info")
})
