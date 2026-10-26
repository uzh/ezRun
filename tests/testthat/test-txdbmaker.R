## GenomicFeatures 1.64 (R 4.6) no longer has makeTxDbFromGFF: it moved to txdbmaker. A bare call
## fails unless txdbmaker happens to be attached, and ChIPSeqPeakComparison's tryCatch turned that
## into a silently missing peak annotation.
test_that("every makeTxDbFromGFF call names txdbmaker", {
  files <- list.files(test_path("..", "..", "R"), pattern = "\\.R$", full.names = TRUE)
  skip_if(length(files) == 0, "package source tree not available")
  expect_gt(length(files), 50)   # positive control: the scan reads the package's R files
  bare <- unlist(lapply(files, function(f) {
    l <- readLines(f, warn = FALSE)
    hit <- grepl("makeTxDbFromGFF\\(", l) & !grepl("txdbmaker::makeTxDbFromGFF\\(", l) & !grepl("^\\s*#", l)
    if (any(hit)) paste0(basename(f), ":", which(hit))
  }))
  expect_identical(bare, NULL, label = paste(bare, collapse = ", "))
})
