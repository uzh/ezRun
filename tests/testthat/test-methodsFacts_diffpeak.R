## methods_facts() of EzAppChIPSeqPeakComparison and EzAppDiffPeakAnalysis: gates follow the
## run's parameters / input dataset, and every literal a fact states is still in the code.

## deparse breaks long calls across lines: collapse all whitespace so fixed matches survive
codeOf <- function(f) gsub("\\s+", " ", paste(deparse(f, width.cutoff = 500L), collapse = " "))
has <- function(facts, pattern, fixed = TRUE) any(grepl(pattern, facts, fixed = fixed))
## one input row whose file columns point at files under a temporary dataRoot (missing = paths not created)
withInput <- function(param, cols, missing = character(0)) {
  root <- tempfile("root"); dir.create(root)
  vals <- stats::setNames(as.list(paste0("f", seq_along(cols))), cols)
  for (i in seq_along(cols)) if (!cols[i] %in% missing) file.create(file.path(root, vals[[i]]))
  param$dataRoot <- root
  attr(param, "input") <- as.data.frame(vals, check.names = FALSE)
  param
}
mm <- "Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-03"

## ---- DiffPeakAnalysis -------------------------------------------------------
diffRun <- list(refBuild = mm, grouping = "Condition", sampleGroup = "KO_Hypo", refGroup = "KO_Norm",
                grouping2 = "", annotationMethod = "homer", normMethod = "DESeq2", lfcThreshold = 1,
                fdrThreshold = 0.05, lfcTest = FALSE, lfcShrink = "apeglm", fitAllSamples = FALSE,
                runMotifs = TRUE, motifDeNovo = FALSE, runGoOra = TRUE, enrichMaxPeaks = 2000)

test_that("DiffPeakAnalysis: unconditional facts only for list()", {
  f <- EzAppDiffPeakAnalysis$new()$methods_facts(list())
  expect_gte(length(f), 4); expect_lte(length(f), 12)
  expect_true(has(f, "blind design ~ 1"))
  expect_true(has(f, "unshrunken"))
  expect_true(has(f, "Enrichr buttons only send"))
  expect_true(has(f, "varianceStabilizingTransformation"))
  expect_false(has(f, "findMotifsGenome|enrichGO|annotatePeak|normMethod|lfcShrink", fixed = FALSE))
  expect_true(has(f, "When the input dataset had BigWig tracks"))      # input unknown
})

test_that("DiffPeakAnalysis: gated facts follow the run's parameters", {
  app <- EzAppDiffPeakAnalysis$new()
  f <- app$methods_facts(diffRun)
  expect_lte(length(f), 12)
  expect_true(has(f, "only the samples of sampleGroup and refGroup (fitAllSamples false)"))
  expect_true(has(f, "design ~ group (no second factor)"))
  expect_true(has(f, "median-of-ratios"))
  expect_true(has(f, "log2 fold change = 0 (lfcTest false)"))
  expect_true(has(f, "at least 1 in absolute value"))
  expect_true(has(f, "below 0.05"))
  expect_true(has(f, "apeglm-shrunken"))
  expect_true(has(f, "log2FoldChange_shrunk"))
  expect_true(has(f, "HOMER annotatePeaks.pl"))
  expect_true(has(f, "motif set vertebrates"))
  expect_true(has(f, "top 2000 up and down candidates"))
  expect_true(has(f, "-nomotif"))
  expect_true(has(f, "enrichGO (org.Mm.eg.db)"))
  expect_true(has(f, "When the input dataset had BigWig tracks"))       # input unknown

  off <- app$methods_facts(modifyList(diffRun, list(lfcShrink = "none", runMotifs = FALSE, runGoOra = FALSE,
                                                    fitAllSamples = TRUE, grouping2 = "Batch", lfcTest = TRUE,
                                                    normMethod = "TMM", annotationMethod = "chipseeker",
                                                    fdrThreshold = 0.01, lfcThreshold = 0.5, enrichMaxPeaks = 500)))
  expect_false(has(off, "apeglm|ashr|log2FoldChange_shrunk", fixed = FALSE))
  expect_true(has(off, "No log2 fold change shrinkage was computed (lfcShrink none)"))
  expect_false(has(off, "findMotifsGenome"))
  expect_false(has(off, "enrichGO"))
  expect_true(has(off, "No GO over-representation was computed (runGoOra false)"))
  expect_true(has(off, "every sample with a label in the grouping column"))
  expect_true(has(off, "~ group + Batch when"))
  expect_true(has(off, "|log2 fold change| > 0.5 (lfcTest true, altHypothesis greaterAbs)"))
  expect_true(has(off, "edgeR TMM"))
  expect_true(has(off, "-1000 to +1000 bp"))
  expect_true(has(off, "below 0.01"))
  expect_true(has(off, "at least 0.5 in absolute value"))

  expect_true(has(app$methods_facts(modifyList(diffRun, list(lfcShrink = "ashr"))), "ashr-shrunken"))
  expect_true(has(app$methods_facts(modifyList(diffRun, list(motifDeNovo = TRUE))), "de novo motif discovery also ran"))
  expect_true(has(app$methods_facts(modifyList(diffRun, list(annotationMethod = "chippeakanno"))), "nearestStart"))
  zf <- app$methods_facts(modifyList(diffRun, list(refBuild = "Danio_rerio/Ensembl/GRCz11/Annotation/Release_110")))
  expect_false(has(zf, "enrichGO"))
  expect_true(has(zf, "only mouse and human"))
  expect_true(has(zf, "motif set vertebrates"))
  expect_true(has(app$methods_facts(modifyList(diffRun, list(refBuild = "Arabidopsis_thaliana/Ensembl/TAIR10"))),
                  "motif set all"))
  ## BigWig profiles only when the input dataset has a BigWig column
  expect_true(has(app$methods_facts(withInput(diffRun, c("Name", "Count [Link]", "BigWig"))), "BigWig signal profiles"))
  expect_false(has(app$methods_facts(withInput(diffRun, c("Name", "Count [Link]"))), "BigWig"))
  expect_false(has(app$methods_facts(withInput(diffRun, c("Name", "Count [Link]", "BigWig"), missing = "BigWig")), "BigWig"))
})

test_that("DiffPeakAnalysis: every literal a fact states is still in the code", {
  f <- paste(EzAppDiffPeakAnalysis$new()$methods_facts(withInput(diffRun, c("Name", "Count [Link]", "BigWig"))), collapse = "\n")
  for (x in c("rounded to integers", "(-size 200)", "fewer than 20", "50,000 with seed 42", "p-value cutoff 0.05, q-value cutoff 0.2",
              "gene sets of 10 to 500 genes", "fewer than 10 genes", "up to 500 candidate gene names",
              "50 bp bins over plus or minus 2 kb", "top 1000 up and down", "blind = FALSE with replicates"))
    expect_match(f, x, fixed = TRUE)
  expect_match(codeOf(generateDESeqDS), "round(countData)", fixed = TRUE)
  r <- codeOf(runDiffPeakDESeq)
  for (s in c("keep <- as.character(dds$group) %in% c(sampleGroup, refGroup)", "~group + grouping2",
              "design(blind) <- ~1", "nbinomWaldTest", "useLfcTest <- isTRUE(lfcTest) && !noReplicates",
              "altHypothesis = \"greaterAbs\"", "!noReplicates && lfcShrinkType %in% c(\"apeglm\", \"ashr\")",
              "setDiffPeakSizeFactors(dds, normMethod)"))
    expect_match(r, s, fixed = TRUE)
  s <- codeOf(diffPeakSizeFactors)
  for (x in c("estimateSizeFactorsForMatrix", "method = \"TMM\"", "exp(mean(log(sizeFactors)))", "readsInPeaks = libSize"))
    expect_match(s, x, fixed = TRUE)
  sh <- codeOf(shrinkDiffPeakLfc)
  expect_match(sh, "type = \"apeglm\"", fixed = TRUE); expect_match(sh, "type = \"ashr\"", fixed = TRUE)
  t <- codeOf(makeDiffPeakTable)
  expect_match(t, "abs(log2FoldChange) >= lfcThreshold", fixed = TRUE); expect_match(t, "padj < fdrThreshold", fixed = TRUE)
  a <- codeOf(annotateConsensusPeaks)
  for (x in c("annotatePeaks.pl", "\"-gtf\"", "tssRegion = c(-1000, 1000)", "output = \"nearestStart\"", "multiple = FALSE"))
    expect_match(a, x, fixed = TRUE)
  expect_equal(formals(runHomerKnownMotifs)[c("minPeaks", "maxBackground", "size")], list(minPeaks = 20, maxBackground = 50000, size = 200))
  h <- codeOf(runHomerKnownMotifs)
  expect_match(h, "set.seed(42)", fixed = TRUE); expect_match(h, "\"-nomotif\"", fixed = TRUE)
  expect_match(h, "\"-size\", size", fixed = TRUE)
  expect_match(codeOf(topDiffPeaks), "order(cand$padj, -abs(cand$log2FoldChange)", fixed = TRUE)
  g <- codeOf(diffPeakGoOra)
  for (x in c("Mus_musculus = \"org.Mm.eg.db\"", "Homo_sapiens = \"org.Hs.eg.db\"", "pAdjustMethod = \"BH\"",
              "pvalueCutoff = 0.05", "qvalueCutoff = 0.2", "minGSSize = 10", "maxGSSize = 500", "length(genes) < 10"))
    expect_match(g, x, fixed = TRUE)
  expect_identical(formals(diffPeakGoOra)$ont, "BP")
  expect_match(codeOf(diffPeakVST), "blind = isTRUE(noReplicates)", fixed = TRUE)
  expect_equal(formals(bigwigPeakProfiles)[c("maxPeaks", "extend", "binSize")], list(maxPeaks = 1000, extend = 2000, binSize = 50))
  expect_match(codeOf(diffPeakNormDiagnostic), "p.adjust(p, method = \"BH\")", fixed = TRUE)
  m <- codeOf(ezMethodDiffPeakAnalysis)
  expect_match(m, "isTRUE(param$runMotifs)", fixed = TRUE); expect_match(m, "isTRUE(param$runGoOra)", fixed = TRUE)
  expect_match(m, "intersect(c(\"BigWig\", \"BigWigFile\"), input$colNames)", fixed = TRUE)
  qmd <- paste(readLines(system.file("templates", "DiffPeakAnalysis.qmd", package = "ezRun")), collapse = "\n")
  expect_match(qmd, "maxGenes = 500L", fixed = TRUE)
  expect_match(qmd, "https://maayanlab.cloud/Enrichr/enrich", fixed = TRUE)
  ## the shrunken estimate appears only in the setup flag and the results table, never in a plot or ranking
  shrunkLines <- grep("_shrunk", strsplit(qmd, "\n")[[1]], value = TRUE)
  expect_length(shrunkLines, 5)
  expect_true(all(grepl("hasShrunk|displayCols|is\\.na\\(peakTable\\$log2FoldChange_shrunk\\)|\"log2FoldChange_shrunk\", \"lfcSE\"", shrunkLines)))
})

## ---- ChIPSeqPeakComparison --------------------------------------------------
chipRun <- list(refBuild = mm, paired = TRUE, markType = "narrow", topN = 5, rankBy = "signalValue",
                minSamplesForConsensus = 1, minOverlapBp = 1, profileExtend = 2000, profileBinSize = 50,
                quantifyFrom = "bigwig", normalization = "CPM", peakFormat = "auto", useBlacklist = FALSE,
                blacklistFile = "", runDifferentialBinding = FALSE)
chipCols <- c("Name", "CalledPeaks [File]", "BED [File]", "BigWigFile [File]", "BAM [File]", "Condition [Factor]")

test_that("ChIPSeqPeakComparison: unconditional facts only for list()", {
  f <- EzAppChIPSeqPeakComparison$new()$methods_facts(list())
  expect_gte(length(f), 4); expect_lte(length(f), 12)
  expect_true(has(f, "called no peaks"))
  expect_true(has(f, "only picked the loci drawn as per-sample coverage tracks"))
  expect_true(has(f, "base-pair Jaccard"))
  expect_true(has(f, "No differential binding test"))
  expect_true(has(f, "When BAM files were given"))                      # input unknown
  expect_false(has(f, "blacklist|CPM|good at", fixed = FALSE))
})

test_that("ChIPSeqPeakComparison: gated facts follow the run's parameters and input", {
  app <- EzAppChIPSeqPeakComparison$new()
  f <- app$methods_facts(withInput(chipRun, chipCols))
  expect_lte(length(f), 12)
  expect_true(has(f, "markType (narrow)"))
  expect_true(has(f, "No blacklist filtering was applied (useBlacklist false)"))
  expect_true(has(f, "FRiP acceptable at 0.01, good at 0.05"))
  expect_true(has(f, "mean BigWig score"))
  expect_true(has(f, "multiplied by one million (CPM"))
  expect_true(has(f, "topN (5)"))
  expect_true(has(f, "highest-signalValue"))
  expect_true(has(f, "at least 1 sample(s)"))
  expect_true(has(f, "FRiP was the number of BAM records"))
  expect_false(has(f, "When BAM files were given"))

  b <- app$methods_facts(withInput(modifyList(chipRun, list(useBlacklist = TRUE, blacklistFile = "/x/bl.bed", markType = "broad",
                                                            quantifyFrom = "bam", normalization = "quantile", topN = 7,
                                                            minSamplesForConsensus = 3, minOverlapBp = 10)), chipCols))
  expect_true(has(b, "Peaks overlapping the regions of blacklistFile were removed"))
  expect_false(has(b, "No blacklist"))
  expect_true(has(b, "markType (broad)"))
  expect_false(has(b, "FRiP acceptable"))                                # broad marks have no FRiP threshold
  expect_true(has(b, "NSC acceptable at 1.02, good at 1.05"))
  expect_true(has(b, "summarizeOverlaps"))
  expect_true(has(b, "normalizeQuantiles"))
  expect_true(has(b, "topN (7)"))
  expect_true(has(b, "at least 3 sample(s) by at least 10 bp"))
  expect_true(has(app$methods_facts(modifyList(chipRun, list(useBlacklist = TRUE, blacklistFile = ""))),
                  "No blacklist filtering was applied (no blacklistFile)"))
  ## no BAM in the input: no BAM-based QC facts, BigWig used even with quantifyFrom bam
  nb <- app$methods_facts(withInput(modifyList(chipRun, list(quantifyFrom = "bam")), setdiff(chipCols, "BAM [File]")))
  expect_false(has(nb, "FRiP|NSC|summarizeOverlaps", fixed = FALSE))
  expect_true(has(nb, "mean BigWig score"))
  ## no BigWig in the input: no BigWig correlation, BAM counts
  nw <- app$methods_facts(withInput(chipRun, setdiff(chipCols, "BigWigFile [File]")))
  expect_false(has(nw, "BigWig score"))
  expect_true(has(nw, "summarizeOverlaps"))
  ## BAM column listed but its files absent: no BAM metrics
  expect_false(has(app$methods_facts(withInput(chipRun, chipCols, missing = "BAM [File]")), "FRiP"))
})

test_that("ChIPSeqPeakComparison: every literal a fact states is still in the code", {
  f <- paste(EzAppChIPSeqPeakComparison$new()$methods_facts(withInput(chipRun, chipCols)), collapse = "\n")
  for (x in c("Peaks, CalledPeaks, BED or MACS", "first 300,000 records", "10 bp bins with shifts up to 1000 bp",
              "first 100,000 records", "1000 bp genome bins", "multiplied by one million", "-3000 to +3000 bp"))
    expect_match(f, x, fixed = TRUE)
  m <- codeOf(ezMethodChIPSeqPeakComparison)
  for (x in c("c(\"Peaks\", \"CalledPeaks\", \"BED\", \"MACS\")", "flagQcMetrics(qcWide, markType",
              "isTRUE(param$useBlacklist) && nzchar(param$blacklistFile", "topSignificantPeaks(peaks, n = topN, rankBy = rankBy)",
              "genomeFootprint(peaks, seqinfo)", "lapply(peaks, cumulativeSignalCurve, rankBy = rankBy)",
              "quantFrom == \"bam\" && hasBam", "makeTxDbFromGFF", "error = function(e) signalMat"))
    expect_match(m, x, fixed = TRUE)
  expect_false(grepl("runDifferentialBinding|ctrlFiles\\[", m))           # read nowhere / never used
  hp <- codeOf(harmonizePeaks)
  for (x in c("off-reference contig", "zero-width", "out-of-bounds", "blacklisted")) expect_match(hp, x, fixed = TRUE)
  expect_false(grepl("keepStandardChromosomes", codeOf(body(harmonizePeaks)), fixed = TRUE))
  fr <- codeOf(body(computeFrip))
  expect_match(fr, "idxstatsBam", fixed = TRUE); expect_match(fr, "countBam", fixed = TRUE)
  expect_match(fr, "reduce(peaks, ignore.strand = TRUE)", fixed = TRUE)
  expect_false(grepl("paired", fr, fixed = TRUE))                       # read pairs are not collapsed
  expect_match(codeOf(libraryComplexity), "repeat", fixed = TRUE)
  expect_equal(formals(crossCorrelationMetrics)[c("binSize", "maxShift", "nSub")], list(binSize = 10L, maxShift = 1000L, nSub = 3e5))
  expect_equal(formals(alignmentStats)$nSub, 1e5)
  expect_match(codeOf(alignmentStats), "mapq[mapped] >= 30", fixed = TRUE)
  expect_match(codeOf(.scanBamSubsample), "yieldSize = nSub", fixed = TRUE)
  expect_equal(formals(bigwigCorrelation)$binSize, 1000L)
  expect_false(grepl("minOverlapBp", codeOf(body(pairwiseJaccard)), fixed = TRUE))
  bc <- codeOf(buildConsensus)
  for (x in c("reduce(pooled, ignore.strand = TRUE)", "minoverlap = minOverlapBp", "occCount >= minSamples"))
    expect_match(bc, x, fixed = TRUE)
  q <- codeOf(quantifyRegions)
  expect_match(q, "mode = \"Union\"", fixed = TRUE); expect_match(codeOf(.scoreRegionsBw), "viewMeans", fixed = TRUE)
  n <- codeOf(normalizeSignal)
  expect_match(n, "1e+06", fixed = TRUE); expect_match(n, "limma::normalizeQuantiles", fixed = TRUE)
  expect_match(codeOf(topSignificantPeaks), "overlapsAny(pooled[idx], chosen", fixed = TRUE)
  expect_equal(formals(annotateConsensus)$tssRegion, quote(c(-3000, 3000)))
  ## the flag thresholds quoted come from the app's threshold table
  thr <- qcThresholdTable("narrow")
  expect_equal(thr$min[thr$metric == "FRiP"], 0.01); expect_equal(thr$good[thr$metric == "FRiP"], 0.05)
})
