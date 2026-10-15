## clusterProfiler 4.21 (R 4.6) hands GSEA to enrichit 0.2.0, which keeps terms on the raw p-value
## only, and keeps every tested ORA term in @result. fdrThreshGSEA / fdrThreshORA are FDR cut-offs
## (the report says "fdr <="). Every tested term is stored (the Excel export lists them all); the
## `significant` column must mark exactly the terms with p.adjust at or below the cut-off, and a
## truly enriched term must be flagged.
goTestInput <- function(seed) {
  set.seed(seed)
  genes <- paste0("ENSG", sprintf("%05d", 1:400))
  ids <- AnnotationDbi::keys(GO.db::GO.db)[1:80]
  log2Ratio <- setNames(rnorm(400), genes)
  top <- names(sort(log2Ratio, decreasing = TRUE))[1:25]            # planted: the 25 most up genes
  up <- sample(setdiff(genes, top), 35)
  sets <- c(list(top), lapply(ids[-1], function(i) sample(genes, 25)))
  go2gene <- data.frame(term = rep(ids, each = 25), gene = unlist(sets))
  seqAnno <- data.frame(gene_id = genes, gene_name = paste0("G", 1:400), isPresent = TRUE,
                        log2Ratio = log2Ratio, pValue = runif(400), row.names = genes)
  list(seqAnno = seqAnno, presentGenes = genes, planted = ids[1],
       go2gene = list(BP = go2gene, MF = go2gene, CC = go2gene),
       selections = list(upGenes = c(top, up), downGenes = sample(genes, 60)))
}

test_that("GSEA flags as significant exactly the terms with p.adjust <= fdrThreshGSEA", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  input <- goTestInput(7)
  res <- ezGSEA(input, list(fdrThreshGSEA = 0.2, rankMetric = "log2Ratio"))
  for (onto in names(res)) {
    tab <- res[[onto]]@result
    expect_identical(tab$significant, tab$p.adjust <= 0.2, label = onto)
    expect_true(any(!tab$significant), label = paste(onto, "non-significant terms kept"))
    expect_true(input$planted %in% tab$ID[tab$significant], label = paste(onto, "planted term flagged"))
  }
})

test_that("ORA flags as significant exactly the terms with p.adjust <= fdrThreshORA", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  input <- goTestInput(11)
  res <- ezEnricher(input, list(fdrThreshORA = 0.05))
  for (onto in names(res)) for (sel in names(res[[onto]])) {
    tab <- res[[onto]][[sel]]@result
    expect_identical(tab$significant, tab$p.adjust <= 0.05, label = paste(onto, sel))
  }
  up <- res$BP$upGenes@result
  expect_true(input$planted %in% up$ID[up$significant], label = "planted term flagged for upGenes")
})

test_that("GSEA stops on an unknown rankMetric", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("GO.db")
  expect_error(ezGSEA(goTestInput(3), list(fdrThreshGSEA = 0.2, rankMetric = "signedpvalue")), "unsupported rankMetric")
})
