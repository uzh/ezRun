## methods_facts() of the count apps: CountQC, CountSpacer, MageckCount, SplitPipe.
## Gated facts follow the typed param; anchors check that each literal a fact states
## is still in the code of the function the fact names.

codeOf <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = "\n")
templateOf <- function(name) paste(readLines(system.file("templates", name, package = "ezRun")), collapse = "\n")
has <- function(facts, pattern, fixed = TRUE) any(grepl(pattern, facts, fixed = fixed))

test_that("CountQC: unconditional facts, GO only with runGO, run values carried", {
  app <- EzAppCountQC$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 4)
  expect_lte(length(f0), 12)
  expect_true(has(f0, "is not a threshold and removed no genes"))
  expect_true(has(f0, "nSampleClusters sets this number of gene clusters"))
  expect_true(has(f0, "ezRun default 10"))
  expect_false(has(f0, "GO over-representation"))                 # param unknown: no GO fact

  on <- list(runGO = TRUE, backgroundExpression = 4, topGeneSize = 50, selectByFtest = FALSE,
             normMethod = "logMean", sigThresh = 12, useSigThresh = TRUE, minSignal = 5,
             maxGenesForClustering = 1500, highVarThreshold = 0.7, nSampleClusters = 5,
             minGenesForClustering = 30, minCountFisher = 3, pValThreshFisher = 1e-4)
  f <- app$methods_facts(on)
  expect_true(has(f, "backgroundExpression (4) is not a threshold"))
  expect_true(has(f, "topGeneSize (50)"))
  expect_true(has(f, "highest standard deviation of the log2 signal across samples (selectByFtest false)"))
  expect_true(has(f, "count exceeded sigThresh (12)"))
  expect_true(has(f, "cut into 5 gene clusters"))
  expect_true(has(f, "up to 1500 genes"))
  expect_true(has(f, "at or below 0.7"))
  expect_true(has(f, "geometric mean of its counts over the genes present in all samples"))
  go <- grep("GO over-representation", f, value = TRUE)
  expect_length(go, 1)
  expect_match(go, "hypergeometric method and no gene-length bias correction", fixed = TRUE)
  expect_match(go, "p below 1e-04 (pValThreshFisher)", fixed = TRUE)
  expect_true(has(f, "and the GO analysis do not use them"))

  off <- modifyList(on, list(runGO = FALSE, selectByFtest = TRUE, useSigThresh = FALSE, normMethod = "quantile"))
  g <- app$methods_facts(off)
  expect_false(has(g, "GO", fixed = TRUE))
  expect_true(has(g, "smallest one-way ANOVA F-test p-value"))
  expect_true(has(g, "count exceeded 0 (useSigThresh false)"))
  expect_true(has(g, "the quantile method (normMethod)"))
  expect_false(has(g, "geometric mean of its counts"))

  ## Enrichr table only when the dataset is not smRNA (the chunk's eval)
  sm <- on; attr(sm, "input") <- data.frame(featureLevel = "smRNA")
  ge <- on; attr(ge, "input") <- data.frame(featureLevel = "gene")
  expect_false(has(app$methods_facts(sm), "Enrichr"))
  expect_true(has(app$methods_facts(ge), "Enrichr"))
  expect_true(has(f0, "Enrichr"))                                  # input unknown: stated
})

test_that("CountQC facts are anchored in the code they describe", {
  expect_match(codeOf(ezMethodCountQC), 'qmdFile = "CountQC.qmd"', fixed = TRUE)
  expect_match(codeOf(loadCountDataset), "presentFlag = counts > sigThresh", fixed = TRUE)
  expect_match(codeOf(ezLogmeanScalingFactor), "isAllPresent", fixed = TRUE)
  expect_match(codeOf(ezLogmeanScalingFactor), "target <- exp(mean(means))", fixed = TRUE)
  expect_match(codeOf(clusterPheatmap), 'method = "ward.D2"', fixed = TRUE)
  expect_match(codeOf(clusterPheatmap), "cutree(clusterInfo$tree_row, nClusters)", fixed = TRUE)
  expect_match(codeOf(clusterPheatmap), 'scale = "none"', fixed = TRUE)
  expect_match(codeOf(goClusterResults), 'method = "Hypergeometric"', fixed = TRUE)
  expect_match(codeOf(goClusterResults), "normalizedAvgSignal = NULL", fixed = TRUE)
  expect_match(codeOf(goClusterResults), 'ontologies = c("BP", "MF", "CC")', fixed = TRUE)
  expect_match(codeOf(ezGoseq), 'p.adjust(pvalues, method = "fdr")', fixed = TRUE)
  expect_match(codeOf(ezGoseq), "lengths(go2GenesList) >= param$minCountFisher", fixed = TRUE)
  expect_match(codeOf(.getGoTermsAsTd), "maxNumberOfTerms = 40", fixed = TRUE)
  expect_match(codeOf(ezMdsPlotly), "plotMDS(logSignal, plot = FALSE)", fixed = TRUE)
  qmd <- templateOf("CountQC.qmd")
  for (a in c("assays(rawData)$log2signal <- log2(assays(rawData)$signal +", "param$backgroundExpression)",
              "method = param$normMethod) + param$minSignal",
              "rowMeans(isPresentCond) >= 0.5", "!seqAnno$IsControl",
              'anova(fit)["Group", "Pr(>F)"]', "decreasing = TRUE), param$topGeneSize)",
              'hclust(d, method = "ward.D2")', "as.dist(1 - cor(",
              "if (ncol(rawData) > 3) {", "param$maxGenesForClustering",
              "<= param$highVarThreshold", "length(varFeatures) > param$minGenesForClustering",
              "hcl.colors(param$nSampleClusters", "universeProbeIds = rownames(seqAnno)",
              "maxGenes <- 500", "/srv/GT/databases/HRT/Human_Mouse_Common.csv", "n = maxGenes / 2",
              "unique(dataset$featureLevel) != 'smRNA'"))
    expect_match(qmd, a, fixed = TRUE, label = a)
  expect_no_match(qmd, "glmFit|DESeq|exactTest|t\\.test|lmFit")     # no DE test in CountQC
})

test_that("CountSpacer: unconditional facts, guessed patterns and spacer length follow the param", {
  app <- EzAppCountSpacer$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 4)
  expect_lte(length(f0), 12)
  expect_true(has(f0, "bowtie (version 1, not Bowtie2)"))
  expect_true(has(f0, "raw read counts"))
  expect_false(has(f0, "inferred from the per-position base composition"))   # param unknown

  p <- list(guessPatterns = TRUE, leftPattern = "", rightPattern = "", spacerLength = 0L,
            maxMismatch = 1, minReadLength = 18L, diffToLogMeanThreshold = 2, dictPath = "LibX")
  f <- app$methods_facts(p)
  expect_true(has(f, "inferred from the per-position base composition"))
  expect_true(has(f, "the library's single sgRNA length (spacerLength 0)"))
  expect_true(has(f, "up to 1 mismatches (maxMismatch)"))
  expect_true(has(f, "shorter than 18 bases (minReadLength)"))
  expect_true(has(f, "mean plus or minus 2 (diffToLogMeanThreshold)"))
  given <- modifyList(p, list(leftPattern = "G", rightPattern = "GTTTTAGAGCTA", spacerLength = 20L))
  g <- app$methods_facts(given)
  expect_false(has(g, "inferred from the per-position base composition"))
  expect_true(has(g, "a fixed spacer length of 20 (spacerLength)"))
  expect_false(has(app$methods_facts(modifyList(p, list(guessPatterns = FALSE))), "inferred from the per-position"))
})

test_that("CountSpacer facts are anchored in the code they describe", {
  m <- codeOf(ezMethodCountSpacer)
  for (a in c("final.csv$", "invert = TRUE", "ezMethodFastpTrim", 'getColumn("Read1")',
              "-f -p", "pattern = \"ebwt$\"", "width(reads) >= param$minReadLength",
              "table(result_bowtie$target)", "dict[is.na(dict$Count), \"Count\"] = 0",
              "log2(1 + sort(dict$Count[!dict$isControl]))", "param$diffToLogMeanThreshold",
              "dict2$Count > 2^lowerCutOff", "[[\"#sgRNAs > lowerCutOff\"]] < 2",
              "tapply(dict2$Count, dict2$TargetID, sum)", 'qmdFile = "CountSpacer.qmd"',
              "length(seqLengths) == 1"))
    expect_match(m, a, fixed = TRUE, label = a)
  gf <- formals(guessFlankingPatterns)
  expect_equal(gf$nSample, 1e5); expect_equal(gf$constFreq, 0.9); expect_equal(gf$varFreq, 0.5)
  expect_equal(gf$minSpacerRun, 10L); expect_equal(gf$maxFlank, 12L)
  t <- codeOf(twoPatternReadFilter)
  expect_match(t, "max.mismatch = maxMismatch", fixed = TRUE)
  expect_match(t, "spStart <- rightStart - spacerLength", fixed = TRUE)
  expect_match(t, "spEnd <- leftEnd + spacerLength", fixed = TRUE)
  qmd <- templateOf("CountSpacer.qmd")
  expect_match(qmd, "mappingRate   <- 100 * mapped / max(validReads, 1)", fixed = TRUE)
  expect_match(qmd, "giniVal    <- .gini(res$Count)", fixed = TRUE)
})

test_that("MageckCount: facts, cmdOptions carried, no normalization or test claimed", {
  app <- EzAppMageckCount$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 4)
  expect_lte(length(f0), 12)
  expect_true(has(f0, "ezRun sets no normalization method"))
  expect_true(has(f0, "rendered by ezRun"))
  expect_true(has(f0, "no mageck test or mle"))
  expect_false(has(f0, "size factor|normalized to|normalised to", fixed = FALSE))
  expect_true(has(app$methods_facts(list(cmdOptions = "--sgrna-len 19")), "cmdOptions '--sgrna-len 19'"))
  expect_true(has(app$methods_facts(list(cmdOptions = "")), "cmdOptions was empty"))
  expect_false(has(app$methods_facts(list(cmdOptions = "")), "cmdOptions '"))
})

test_that("MageckCount facts are anchored in the code they describe", {
  m <- codeOf(ezMethodMageckCount)
  for (a in c("mageck2 count", "\"-l\"", "--control-sgrna", "--fastq", "\"-n\"", 'getFullPaths("Read1")',
              "hasCtrl <-", 'qmdFile = "MageckCountQC.qmd"', "tryCatch(", "file.exists(summaryFile)"))
    expect_match(m, a, fixed = TRUE, label = a)
  expect_no_match(m, "norm-method|trim-5|sgrna-len|mageck2 test|mageck2 mle|fastp")
  r <- codeOf(getMageckReference)
  expect_match(r, "_MAGeCK_Ctrl", fixed = TRUE)
  expect_match(r, "prepareMageckLibrary", fixed = TRUE)
  p <- codeOf(prepareMageckLibrary)
  expect_match(p, 'c("true", "t", "1", "yes", "y")', fixed = TRUE)
  expect_match(p, 'sep = "_"', fixed = TRUE)
  expect_match(p, 'c("ID", "Sequence", "GeneSymbol")', fixed = TRUE)
  qmd <- templateOf("MageckCountQC.qmd")
  expect_match(qmd, "GiniIndex", fixed = TRUE)
})

test_that("SplitPipe: sublibrary, combine, reference and sample facts follow param and input", {
  app <- EzAppSplitPipe$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 4)
  expect_lte(length(f0), 12)
  expect_true(has(f0, "one sublibrary of the same experiment, not a biological sample"))
  expect_true(has(f0, "When the input had more than one sublibrary"))
  expect_true(has(f0, "write10xCounts"))
  expect_true(has(f0, "rendered by ezRun"))

  p <- list(kit = "WT", chemistry = "v4", sampleWells = "all-well A1-A12", sampleLoadingTable = "",
            sublibDir = "", transcriptTypes = c("protein_coding", "rRNA"), saveAnndata = FALSE,
            cmdOptions = "", cores = 8)
  attr(p, "input") <- data.frame(Name = c("S1", "S2"))
  f <- app$methods_facts(p)
  expect_true(has(f, "--chemistry v4 --kit WT"))
  expect_true(has(f, "The 2 sublibraries were merged with split-pipe --mode comb"))
  expect_true(has(f, "filtered to the transcript types protein_coding, rRNA"))
  expect_true(has(f, "--sample 'all-well A1-A12'"))
  expect_false(has(f, "--save_anndata"))
  expect_true(has(f, "sets no option for"))
  expect_false(has(f, "sublibDir"))

  one <- p; attr(one, "input") <- data.frame(Name = "S1")
  expect_true(has(app$methods_facts(one), "single sublibrary"))
  expect_false(has(app$methods_facts(one), "--mode comb"))
  adopt <- modifyList(p, list(sublibDir = "/srv/GT/analysis/pX/sublibs"))
  a <- app$methods_facts(adopt)
  expect_true(has(a, "were not run in this job"))
  expect_false(has(a, "--chemistry v4"))
  expect_true(has(app$methods_facts(modifyList(p, list(sampleWells = "A A1-A6+B A7-A12"))), "--sample 'A A1-A6', --sample 'B A7-A12'"))
  expect_true(has(app$methods_facts(modifyList(p, list(sampleWells = ""))), "--yes_allwell"))
  expect_true(has(app$methods_facts(modifyList(p, list(sampleLoadingTable = "/x/tab.xlsx"))), "sample loading table"))
  expect_true(has(app$methods_facts(modifyList(p, list(saveAnndata = TRUE))), "--save_anndata"))
})

test_that("SplitPipe facts are anchored in the code they describe", {
  m <- codeOf(ezMethodSplitPipe)
  for (a in c("split-pipe --mode comb", "nSublib > 1", "writeTenxMatrices(resultDir)",
              "makeSplitPipeReport(resultDir, param)", "adoptSplitPipeSublibraries"))
    expect_match(m, a, fixed = TRUE, label = a)
  r <- codeOf(runSplitPipeSublibraries)
  for (a in c("split-pipe --mode all", "--chemistry ", "--kit ", "--genome_dir ", "--save_anndata", "param$cmdOptions"))
    expect_match(r, a, fixed = TRUE, label = a)
  g <- codeOf(getParseReference)
  for (a in c("split-pipe --mode mkref", "gtfByTxTypes", "_Parse_SC_", "--genome_name ", "refFields[3]"))
    expect_match(g, a, fixed = TRUE, label = a)
  s <- codeOf(getParseSampleArgs)
  for (a in c("--yes_allwell", "--samp_sltab", "--samp_list", "\"--sample\"", "[;+]"))
    expect_match(s, a, fixed = TRUE, label = a)
  expect_match(codeOf(adoptSplitPipeSublibraries), "All done split-pipe", fixed = TRUE)
  cv <- codeOf(convertParseToTenx)
  for (a in c("DropletUtils::write10xCounts", 'version = "3"', "gene.id = genes$gene_id", "Matrix::t(m)"))
    expect_match(cv, a, fixed = TRUE, label = a)
  expect_match(codeOf(writeTenxMatrices), "raw_feature_bc_matrix", fixed = TRUE)
  expect_match(codeOf(makeSplitPipeReport), 'qmdFile = "SplitPipe.qmd"', fixed = TRUE)
  qmd <- templateOf("SplitPipe.qmd")
  expect_match(qmd, "_analysis_summary.html", fixed = TRUE)
  expect_match(qmd, "agg_sample_summary.csv", fixed = TRUE)
})
