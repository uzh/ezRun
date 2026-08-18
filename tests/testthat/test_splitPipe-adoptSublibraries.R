context("splitPipe: adopting pre-built sublibraries")

## adoptSplitPipeSublibraries() lets `--mode comb` run over sublibrary
## directories produced by separate jobs. Everything it accepts is then fed
## straight into a combine that takes hours, so each way of being incomplete has
## to stop the run rather than silently produce a partial result.

makeFinishedSublib <- function(root, name, version = "1_8_2", finished = TRUE) {
  d <- file.path(root, name)
  dir.create(file.path(d, "process"), recursive = TRUE, showWarnings = FALSE)
  lines <- c("# Starting split-pipe", "# Some work happened")
  if (finished) {
    lines <- c(lines, "# All done split-pipe v1.8.2")
  }
  writeLines(lines, file.path(d, paste0("split-pipe_v", version, ".log")))
  d
}

test_that("finished sublibraries are symlinked into sublibs/<Name>", {
  src <- tempfile("src")
  dir.create(src)
  for (nm in c("1", "2")) makeFinishedSublib(src, nm)

  wd <- tempfile("wd")
  dir.create(wd)
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)

  dirs <- adoptSplitPipeSublibraries(src, c("1", "2"))

  expect_equal(dirs, file.path("sublibs", c("1", "2")))
  ## The combine command is handed these paths, so they must resolve to the
  ## originals and not to empty placeholders.
  expect_true(all(dir.exists(dirs)))
  expect_equal(Sys.readlink("sublibs/1"), normalizePath(file.path(src, "1")))
  expect_true(file.exists("sublibs/2/split-pipe_v1_8_2.log"))
})

test_that("a missing sublibrary directory stops the run", {
  src <- tempfile("src")
  dir.create(src)
  makeFinishedSublib(src, "1")

  wd <- tempfile("wd")
  dir.create(wd)
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)

  expect_error(
    adoptSplitPipeSublibraries(src, c("1", "2")),
    "no directory for sublibrary '2'"
  )
})

test_that("a sublibrary that never finished stops the run", {
  src <- tempfile("src")
  dir.create(src)
  makeFinishedSublib(src, "1")
  makeFinishedSublib(src, "2", finished = FALSE)

  wd <- tempfile("wd")
  dir.create(wd)
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)

  expect_error(
    adoptSplitPipeSublibraries(src, c("1", "2")),
    "sublibrary '2' did not finish"
  )
})

test_that("a directory with no split-pipe log stops the run", {
  src <- tempfile("src")
  dir.create(src)
  makeFinishedSublib(src, "1")
  dir.create(file.path(src, "2"))

  wd <- tempfile("wd")
  dir.create(wd)
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)

  expect_error(
    adoptSplitPipeSublibraries(src, c("1", "2")),
    "expected exactly one split-pipe log"
  )
})

test_that("the log is found for a split-pipe version other than 1.8.2", {
  src <- tempfile("src")
  dir.create(src)
  makeFinishedSublib(src, "1", version = "9_0_0")

  wd <- tempfile("wd")
  dir.create(wd)
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)

  expect_equal(adoptSplitPipeSublibraries(src, "1"), file.path("sublibs", "1"))
})
