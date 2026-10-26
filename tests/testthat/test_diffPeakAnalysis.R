context("DiffPeakAnalysis helpers")

makeTestDds <- function(nPerGroup = c(A = 3, B = 3, C = 3), nPeaks = 600,
                        shiftB = 0, seed = 1) {
  set.seed(seed)
  groups <- rep(names(nPerGroup), nPerGroup)
  samples <- paste0("s", seq_along(groups))
  mu <- rgamma(nPeaks, shape = 2, scale = 50)
  counts <- sapply(seq_along(groups), function(i) {
    m <- mu
    if (groups[i] == "B") {
      ## 10% of peaks up 4-fold, plus an optional global shift
      m[seq_len(nPeaks / 10)] <- m[seq_len(nPeaks / 10)] * 4
      m <- m * 2^shiftB
    }
    rnbinom(nPeaks, mu = m, size = 20)
  })
  colnames(counts) <- samples
  peakIds <- paste0("peak_", seq_len(nPeaks))
  featureCounts <- data.frame(
    Geneid = peakIds, Chr = "chr1", Start = seq_len(nPeaks) * 1000,
    End = seq_len(nPeaks) * 1000 + 300, Strand = "+", Length = 301,
    counts, check.names = FALSE
  )
  grouping <- setNames(groups, samples)
  generateDESeqDS(featureCounts, c("Geneid", "Chr", "Start", "End", "Strand", "Length"),
    grouping)
}

test_that("default fit is a two-group fit with DESeq2 size factors", {
  dds <- makeTestDds()
  out <- runDiffPeakDESeq(dds, "B", "A")
  expect_equal(out$nFitSamples, 6)
  expect_false(out$noReplicates)
  expect_equal(out$normMethod, "DESeq2")
  tbl <- makeDiffPeakTable(out$res, out$noReplicates, 1, 0.05)
  ## the 60 simulated up peaks should dominate the candidates
  expect_gt(sum(tbl$candidate & tbl$direction == "up" & tbl$peakId %in% paste0("peak_", 1:60)), 40)
  expect_true(all(is.na(tbl$log2FoldChange_shrunk)))
})

test_that("fitAllSamples keeps all groups and returns the same contrast", {
  dds <- makeTestDds()
  out <- runDiffPeakDESeq(dds, "B", "A", fitAllSamples = TRUE)
  expect_equal(out$nFitSamples, 9)
  expect_setequal(levels(out$dds$group), c("A", "B", "C"))
  expect_equal(levels(out$dds$group)[1], "A")
  two <- runDiffPeakDESeq(dds, "B", "A")
  expect_gt(cor(out$res$log2FoldChange, two$res$log2FoldChange, use = "complete"), 0.95)
})

test_that("size-factor methods differ under a global shift", {
  dds <- makeTestDds(shiftB = 1)
  counts <- DESeq2::counts(dds)[, dds$group %in% c("A", "B")]
  sfDeseq <- diffPeakSizeFactors(counts, "DESeq2")
  sfLib <- diffPeakSizeFactors(counts, "readsInPeaks")
  sfTmm <- diffPeakSizeFactors(counts, "TMM")
  expect_equal(exp(mean(log(sfLib))), 1)
  expect_equal(exp(mean(log(sfTmm))), 1)
  expect_equal(unname(sfLib), unname(colSums(counts) / exp(mean(log(colSums(counts))))))
  expect_error(diffPeakSizeFactors(counts, "bogus"), "unsupported normMethod")
  outLib <- runDiffPeakDESeq(dds, "B", "A", normMethod = "readsInPeaks")
  expect_equal(unname(DESeq2::sizeFactors(outLib$dds)), unname(sfLib))
})

test_that("normalisation diagnostic reproduces the used method and shifts others", {
  dds <- makeTestDds(shiftB = 1)
  out <- runDiffPeakDESeq(dds, "B", "A")
  tbl <- makeDiffPeakTable(out$res, out$noReplicates, 1, 0.05)
  diag <- diffPeakNormDiagnostic(out$dds, "B", "A", tbl, 1, 0.05)
  b <- diag$balance
  expect_equal(b$lfcOffset[b$method == "DESeq2"], 0)
  expect_equal(b$up[b$method == "DESeq2"],
    sum(tbl$candidate & tbl$direction == "up" & !is.na(tbl$padj)))
  ## DESeq2 removes the global 2-fold shift, reads-in-peaks keeps part of it
  expect_lt(b$lfcOffset[b$method == "readsInPeaks"], 0)
  expect_setequal(colnames(diag$sizeFactors),
    c("sample", "group", "totalReadsInPeaks", "used", "DESeq2", "readsInPeaks", "TMM"))
})

test_that("shrinkage and lfc test add columns without changing candidates", {
  dds <- makeTestDds()
  out <- runDiffPeakDESeq(dds, "B", "A", lfcShrinkType = "ashr")
  expect_equal(out$lfcShrinkType, "ashr")
  expect_true(any(!is.na(out$res$log2FoldChange_shrunk)))
  base <- runDiffPeakDESeq(dds, "B", "A")
  expect_equal(out$res$log2FoldChange, base$res$log2FoldChange)
  lt <- runDiffPeakDESeq(dds, "B", "A", lfcTest = TRUE, lfcThreshold = 1)
  expect_true(lt$lfcTest)
  expect_lte(sum(lt$res$padj < 0.05, na.rm = TRUE), sum(base$res$padj < 0.05, na.rm = TRUE))
})

test_that("no-replicate comparisons still run", {
  dds <- makeTestDds(nPerGroup = c(A = 1, B = 1))
  out <- runDiffPeakDESeq(dds, "B", "A", lfcShrinkType = "apeglm")
  expect_true(out$noReplicates)
  expect_equal(out$lfcShrinkType, "none")
  tbl <- makeDiffPeakTable(out$res, TRUE, 1, 0.05)
  expect_true(all(tbl$candidate == (abs(tbl$log2FoldChange) >= 1)))
})

test_that("topDiffPeaks ranks by padj and caps", {
  tbl <- data.frame(peakId = paste0("p", 1:5), candidate = TRUE,
    direction = c("up", "up", "down", "up", "down"),
    padj = c(0.01, 0.001, 0.02, 0.03, 0.5), log2FoldChange = c(2, 3, -2, 1.5, -4))
  top <- topDiffPeaks(tbl, "up", 2)
  expect_equal(top$peakId, c("p2", "p1"))
  topNr <- topDiffPeaks(tbl, "down", 5, noReplicates = TRUE)
  expect_equal(topNr$peakId, c("p5", "p3"))
})

test_that("HOMER known results are parsed", {
  known <- data.frame(
    "Motif Name" = c("HIF2a(bHLH)/785_O-HIF2a-ChIP-Seq(GSE34871)/Homer", "CTCF(Zf)/CD4+-CTCF-ChIP-Seq/Homer"),
    "Consensus" = c("GCGTGC", "AYAGTGCCMYCTRGTGGCCA"),
    "P-value" = c("1e-22", "1e-1"),
    "Log P-value" = c(-50.66, -2.3),
    "q-value (Benjamini)" = c("0.0000", "1.0000"),
    "# of Target Sequences with Motif(of 2000)" = c(135, 40),
    "% of Target Sequences with Motif" = c("6.75%", "2.00%"),
    "# of Background Sequences with Motif(of 49912)" = c(1287, 998),
    "% of Background Sequences with Motif" = c("2.58%", "2.00%"),
    check.names = FALSE
  )
  p <- parseHomerKnown(known, nTargets = 2000)
  expect_equal(p$motif, c("HIF2a(bHLH)", "CTCF(Zf)"))
  expect_equal(p$tf, c("HIF2a", "CTCF"))
  expect_equal(p$family, c("bHLH", "Zf"))
  expect_equal(p$log10P[1], 50.66 / log(10))
  expect_equal(p$enrichment[1], 6.75 / 2.58)
  expect_equal(resolveHomerMotifSet("Mus_musculus/GENCODE/GRCm39/Annotation/x"), "vertebrates")
  expect_equal(resolveHomerMotifSet("Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/x"), "all")
})

test_that("gene ids prefer Ensembl gene ids and strip versions", {
  tbl <- data.frame(geneName = c("Egln3", "Bnip3", ""),
    `Entrez ID` = c("ENSMUSG00000035105.5", "ENSMUSG00000078566", NA), check.names = FALSE)
  ids <- diffPeakGeneIds(tbl)
  expect_equal(ids$keyType, "ENSEMBL")
  expect_equal(ids$ids, c("ENSMUSG00000035105", "ENSMUSG00000078566", NA))
  ids2 <- diffPeakGeneIds(tbl[, "geneName", drop = FALSE])
  expect_equal(ids2$keyType, "SYMBOL")
  expect_equal(ids2$ids, c("Egln3", "Bnip3", NA))
})
