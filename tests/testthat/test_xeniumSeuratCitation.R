## citation() must cite rctd-py only when a sample actually ran it: the job falls back to
## spacexr when the rctd-py env is missing, and that is recorded only in each sample's
## log.txt ("RCTD engine: ..."), not in parameters.tsv.

rctdPyCited <- function(param) {
  any(grepl("rctd-py", EzAppXeniumSeurat$new()$citation(param), fixed = TRUE))
}

writeSampleLog <- function(resultDir, sample, engineLine) {
  dir.create(file.path(resultDir, sample), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("2026-11-11 20:00:00 Loaded 1000 cells",
               paste("2026-11-11 20:10:00 RCTD engine:", engineLine)),
             file.path(resultDir, sample, "log.txt"))
}

base <- list(rctdReference = "/srv/GT/databases/RCTD/ref.qs2", rctdEngine = "rctd-py")

test_that("without a result folder the rctdEngine parameter decides", {
  expect_true(rctdPyCited(base))
  expect_false(rctdPyCited(modifyList(base, list(rctdEngine = "spacexr"))))
  expect_false(rctdPyCited(modifyList(base, list(rctdReference = "None"))))
})

test_that("a run that fell back to spacexr does not cite rctd-py", {
  resultDir <- withr::local_tempdir()
  writeSampleLog(resultDir, "S1", "spacexr")
  writeSampleLog(resultDir, "S2", "spacexr")
  expect_false(rctdPyCited(modifyList(base, list(resultDir = resultDir))))
})

test_that("rctd-py is cited when any sample ran it", {
  resultDir <- withr::local_tempdir()
  writeSampleLog(resultDir, "S1", "rctd-py (gi_rctd-py_0.4.3, GPU)")
  writeSampleLog(resultDir, "S2", "spacexr")
  expect_true(rctdPyCited(modifyList(base, list(resultDir = resultDir))))
})

test_that("a result folder without engine lines falls back to the parameter", {
  resultDir <- withr::local_tempdir()
  dir.create(file.path(resultDir, "S1"))
  writeLines("2026-06-06 13:00:00 Loaded 1000 cells", file.path(resultDir, "S1", "log.txt"))
  expect_true(rctdPyCited(modifyList(base, list(resultDir = resultDir))))
})
