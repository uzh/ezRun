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
