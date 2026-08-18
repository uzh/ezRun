## Shared helpers for the ScSeurat QC tests. testthat sources helper-*.R before
## any test file, so both test_scSeurat-doubletQc.R and
## test_scSeurat-madThresholds.R can use these.

## Install a stand-in scDblFinder for the duration of a test.
##
## ezRun attaches scDblFinder with library() and calls it unqualified, so a
## binding in globalenv shadows it: a namespace's lookup order is
## namespace -> imports -> base -> globalenv -> attached packages. That is what
## makes this work without touching the scDblFinder namespace.
localFakeScDblFinder <- function(fn, env = parent.frame()) {
  had <- exists("scDblFinder", envir = globalenv(), inherits = FALSE)
  old <- if (had) get("scDblFinder", envir = globalenv()) else NULL
  assign("scDblFinder", fn, envir = globalenv())
  withr::defer(
    {
      if (had) {
        assign("scDblFinder", old, envir = globalenv())
      } else {
        rm("scDblFinder", envir = globalenv())
      }
    },
    envir = env
  )
}

## a returnType="table" stand-in: rownames are the cells, plus one artificial
## doublet row, exactly as scDblFinder returns
fakeDoubletTable <- function(cellNames) {
  data.frame(
    score = c(seq_along(cellNames) / (length(cellNames) + 1), 0.99),
    class = c(
      rep_len(c("singlet", "doublet"), length(cellNames)),
      "doublet"
    ),
    type = c(rep("real", length(cellNames)), "artificial"),
    row.names = c(cellNames, "artificial1"),
    stringsAsFactors = FALSE
  )
}

## Capture whatever ezLog() emits while evaluating `code`. ezLog goes through the
## logger package, so grab both streams rather than betting on which one.
captureEzLog <- function(code) {
  msg <- NULL
  out <- capture.output(
    msg <- capture.output(force(code), type = "message"),
    type = "output"
  )
  c(out, msg)
}
