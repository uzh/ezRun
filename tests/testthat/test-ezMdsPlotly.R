test_that("ezMdsPlotly handles factor names that are not syntactic", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("plotly")
  set.seed(1)
  logSignal <- matrix(
    rnorm(600),
    ncol = 6,
    dimnames = list(NULL, paste0("S", 1:6))
  )
  design <- data.frame(
    `Time Point` = rep(c("d0", "d7"), 3),
    `Cell-Type` = rep(c("A", "B"), each = 3),
    check.names = FALSE
  )
  traceNames <- function(p) {
    unique(vapply(plotly::plotly_build(p)$x$data, `[[`, "", "name"))
  }
  for (nd in c(2, 3)) {
    p <- ezMdsPlotly(logSignal, design, ndim = nd, main = "mds")
    expect_true(all(c("d0", "d7") %in% unlist(strsplit(traceNames(p), "<br />"))))
    p <- ezMdsPlotly(logSignal, design[, 1, drop = FALSE], ndim = nd, main = "mds")
    expect_setequal(traceNames(p), c("d0", "d7"))
  }
})

test_that("ezMdsPlotly works with default ndim and a design column named samples", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("plotly")
  set.seed(1)
  logSignal <- matrix(
    rnorm(600),
    ncol = 6,
    dimnames = list(NULL, paste0("S", 1:6))
  )
  design <- data.frame(samples = rep(c("a", "b"), 3))
  p <- ezMdsPlotly(logSignal, design, main = "mds")
  built <- plotly::plotly_build(p)
  expect_equal(built$x$data[[1]]$type, "scatter3d")
  expect_setequal(vapply(built$x$data, `[[`, "", "name"), c("a", "b"))
  p <- ezMdsPlotly(logSignal, design, ndim = 2, main = "mds")
  traces <- plotly::plotly_build(p)$x$data
  texts <- setNames(lapply(traces, `[[`, "text"), vapply(traces, `[[`, "", "name"))
  expect_setequal(texts$a, c("S1", "S3", "S5"))
  expect_setequal(texts$b, c("S2", "S4", "S6"))
  expect_error(ezMdsPlotly(logSignal, design, ndim = 1, main = "mds"), "2D")
})
