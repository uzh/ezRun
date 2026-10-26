context("Enrichment result handling: go-analysis.R (ezFlagSignificantTerms, ezSignificantEnrichResult, ezWriteEnrichXlsx)")

makeTab <- function() {
  tab <- data.frame(ID = paste0("T", 1:6), pvalue = c(0.5, 0.0001, 0.03, 0.001, 0.2, 0.004))
  tab$p.adjust <- p.adjust(tab$pvalue, method = "BH")
  tab
}

test_that("ezFlagSignificantTerms flags p.adjust <= threshold and sorts", {
  res <- ezFlagSignificantTerms(makeTab(), 0.01)
  expect_false(is.unsorted(res$p.adjust))
  expect_identical(res$significant, res$p.adjust <= 0.01)
  ## the threshold is inclusive
  thr <- sort(res$p.adjust)[2]
  expect_equal(sum(ezFlagSignificantTerms(makeTab(), thr)$significant), 2L)
})

test_that("ezFlagSignificantTerms stops when rows were dropped after BH", {
  tab <- makeTab()
  expect_error(ezFlagSignificantTerms(tab[tab$pvalue <= 0.05, ], 0.05), "BH adjustment")
})

test_that("ezFlagSignificantTerms handles an empty table", {
  res <- ezFlagSignificantTerms(makeTab()[0, ], 0.05)
  expect_equal(nrow(res), 0L)
  expect_true("significant" %in% colnames(res))
})

test_that("ezSignificantEnrichResult and ezIsEmptyEnrichResult", {
  expect_null(ezSignificantEnrichResult(NULL, 0.05))
  expect_true(ezIsEmptyEnrichResult(NULL))
  expect_true(ezIsEmptyEnrichResult(NA))
  skip_if_not_installed("DOSE")
  er <- new("enrichResult", result = ezFlagSignificantTerms(makeTab(), 0.05))
  sig <- ezSignificantEnrichResult(er, 0.05)
  expect_true(all(sig@result$significant))
  expect_equal(nrow(sig@result), sum(er@result$significant))
})

test_that("ezWriteEnrichXlsx writes all terms, one sheet per table", {
  mf <- makeTab()[1:3, ]
  mf$p.adjust <- p.adjust(mf$pvalue, method = "BH")
  tabs <- list(BP = ezFlagSignificantTerms(makeTab(), 0.05),
               MF = ezFlagSignificantTerms(mf, 0.05))
  f <- tempfile(fileext = ".xlsx")
  on.exit(unlink(f))
  ezWriteEnrichXlsx(tabs, f)
  expect_setequal(openxlsx::getSheetNames(f), c("BP", "MF"))
  bp <- openxlsx::read.xlsx(f, sheet = "BP")
  expect_equal(nrow(bp), 6L)
  expect_equal(sum(bp$significant), sum(tabs$BP$significant))
})

test_that("ezRestrictGoList matches lapply(intersect)", {
  g2g <- list(g1 = c("GO:1", "GO:2", "GO:2", "GO:9"), g2 = "GO:9", g3 = character(0),
              g4 = c("GO:3", NA, "GO:1"), g5 = NULL)
  allGos <- c("GO:1", "GO:2", "GO:3")
  expect_identical(ezRestrictGoList(g2g, allGos),
                   lapply(g2g, function(x) intersect(x, allGos)))
  expect_identical(ezRestrictGoList(list(), allGos), list())
})

## Canary for clusterProfiler behaviour that ezEnricher/ezGSEA rely on: with
## pvalueCutoff = 1 every tested term is returned and p.adjust is BH over them,
## and GSEA is reproducible under a fixed seed.
test_that("clusterProfiler returns complete, reproducible tables with pvalueCutoff = 1", {
  skip_if_not_installed("clusterProfiler")
  set.seed(1)
  genes <- paste0("g", 1:2000)
  t2g <- data.frame(term = rep(paste0("S", 1:40), each = 25),
                    gene = as.vector(replicate(40, sample(genes, 25))))
  gl <- sort(setNames(rnorm(2000), genes), decreasing = TRUE)
  gl[t2g$gene[t2g$term == "S1"]] <- gl[t2g$gene[t2g$term == "S1"]] + 3
  gl <- sort(gl, decreasing = TRUE)

  runGsea <- function() {
    set.seed(42)
    suppressWarnings(clusterProfiler::GSEA(gl, pvalueCutoff = 1, TERM2GENE = t2g,
                                           verbose = FALSE))@result
  }
  g1 <- runGsea()
  expect_equal(nrow(g1), 40L)
  expect_silent(ezFlagSignificantTerms(g1, 0.05))
  expect_identical(g1, runGsea())

  ora <- clusterProfiler::enricher(names(gl)[1:200], universe = genes, TERM2GENE = t2g,
                                   pvalueCutoff = 1, qvalueCutoff = 1)@result
  expect_silent(ezFlagSignificantTerms(ora, 0.05))
})
