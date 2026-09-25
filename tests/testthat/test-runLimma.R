test_that("runLimma returns statistics in the input gene order", {
  set.seed(1)
  x <- matrix(rpois(600, 100), nrow = 100, dimnames = list(paste0("g", 1:100), paste0("s", 1:6)))
  x["g50", 4:6] <- 2000                                   # one strongly up gene, not first
  grouping <- rep(c("ref", "trt"), each = 3)
  for (m in c("limma-trend", "voom")) {
    res <- runLimma(x, "trt", "ref", grouping, priorCount = 3, modelMethod = m)
    expect_identical(res$id, rownames(x), label = m)
    expect_identical(which.min(res$pval), 50L, label = m)
    expect_gt(res$log2FoldChange[50], 3)
  }
})
