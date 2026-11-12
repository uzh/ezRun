## ezMethodXeniumQC skips a region without metrics_summary.csv; when every region lacks it,
## it used to die later with "undefined columns selected". It must stop with a message that
## names the missing files instead.

mockXeniumInput <- function(dataRoot, samples) {
  for (s in samples) dir.create(file.path(dataRoot, "p1", s), recursive = TRUE)
  meta <- data.frame(`XeniumPath [Link]` = file.path("p1", samples),
                     row.names = samples, check.names = FALSE)
  list(meta = meta, getColumn = function(col) meta[[paste0(col, " [Link]")]])
}

test_that("no region with metrics_summary.csv gives a clear error", {
  root <- withr::local_tempdir()
  withr::local_dir(withr::local_tempdir())
  input <- mockXeniumInput(root, c("R1", "R2"))
  expect_error(
    suppressWarnings(ezMethodXeniumQC(input = input, param = list(dataRoot = root))),
    "metrics_summary.csv.*R1.*R2"
  )
})
