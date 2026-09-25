###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodEdger <- function(input = NA, output = NA, param = NA) {
  require(withr)
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  defer(setwd(cwd))

  stopifnot(param$sampleGroup != param$refGroup)

  input <- cleanupTwoGroupsInput(input, param)
  param$groupingName <- param$grouping
  param$grouping <- input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1) {
    param$grouping2Name <- param$grouping2
    param$grouping2 <- input$getColumn(param$grouping2)
    groupNum <- as.numeric(param$grouping2)
    if (all(!is.na(groupNum))) {
      param$grouping2 <- setNames(groupNum, names(param$grouping2))
    }
  }

  rawData <- loadCountDataset(input, param)
  if (isError(rawData)) {
    writeErrorReport("00index.html", param = param, error = rawData$error)
    return("Error")
  }

  deResult <- twoGroupCountComparison(rawData)
  if (isError(deResult)) {
    writeErrorReport("00index.html", param = param, error = deResult$error)
    return("Error")
  }

  makeRmdReport(
    output = output,
    param = param,
    deResult = deResult,
    rmdFile = "twoGroups.Rmd",
    reportTitle = param$comparison
  )
  rmStatus <- file.remove(list.files(pattern = "enrichr-.*rds"))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodEdger(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run a differential expression analysis with the application edgeR on two groups.
EzAppEdger <-
  setRefClass(
    "EzAppEdger",
    contains = "EzApp",
    methods = list(
      ## edgeR (glm/exactTest) and DESeq2/limma are mutually exclusive alternate
      ## testMethod choices; ashr gated on useLfcShrink (deseq2 path only); RUVSeq
      ## gated on runRUV; clusterProfiler/GO.db/Enrichr gated on runGO (Enrichr also
      ## on doEnrichr's organism/featureLevel check) -- listed regardless of gating.
      citation = function() {
        c(
          "Robinson, M.D., McCarthy, D.J. & Smyth, G.K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26(1), 139-140 (2010). https://doi.org/10.1093/bioinformatics/btp616",
          "Ritchie, M.E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47 (2015). https://doi.org/10.1093/nar/gkv007",
          "Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8",
          "Stephens, M. False discovery rates: a new deal. Biostatistics 18(2), 275-294 (2017). https://doi.org/10.1093/biostatistics/kxw041",
          "Risso, D., Ngai, J., Speed, T.P. & Dudoit, S. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature Biotechnology 32(9), 896-902 (2014). https://doi.org/10.1038/nbt.2931",
          "Wu, T. et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation 2(3), 100141 (2021). https://doi.org/10.1016/j.xinn.2021.100141",
          "Bioconductor. GO.db: A set of annotation maps describing the entire Gene Ontology. R package. https://doi.org/10.18129/B9.bioc.GO.db",
          "Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128",
          "Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377"
        )
      },
      ## edgeR defaults quoted here were checked against edgeR 4.10.1 formals()
      ## and function bodies (R 4.6.0 system lib).
      methods_facts = function(param = list()) {
        known <- length(param) > 0
        runGO <- isTRUE(as.logical(param$runGO))
        ## twoGroupCountComparison: a NULL testMethod means glm; deTest is used only for glm
        testMethod <- if (!known) NA else if (is.null(param$testMethod)) "glm" else param$testMethod
        glm <- identical(testMethod, "glm")
        exact <- identical(testMethod, "exactTest")
        g2 <- ezIsSpecified(param$grouping2)
        ## twoGroupCountComparison: robust = ezIsSpecified(param$robust) && param$robust
        robust <- ezIsSpecified(param$robust) && isTRUE(as.logical(param$robust))
        baselines <- ezIsSpecified(param$sampleGroupBaseline) && ezIsSpecified(param$refGroupBaseline)
        c(
          ## ngsio.R:117-127 presentFlag = counts > sigThresh (EZ_PARAM_DEFAULTS sigThresh 10);
          ## twoGroupCountComparison (twoGroups.R:87-95, 143-147)
          "A gene was called present in a sample when its count exceeded sigThresh (ezRun default 10); every gene passing the gene-level transcriptTypes filter (genes whose type is in transcriptTypes, applied after any transcript-to-gene summing) was fitted and tested and has a p-value; presence in at least half of the samples of the sample group or of the reference group only decides which genes enter the reported FDR, the Benjamini-Hochberg adjustment of the edgeR p-values computed by ezRun over the present genes.",
          ## runGlm (twoGroups.R:338-350) / runEdger (twoGroups.R:283-284)
          paste0("Normalization factors were computed with edgeR calcNormFactors using the normMethod method on all genes",
                 if (glm) ", on the samples of the compared groups only" else if (exact) ", on all samples of the dataset",
                 "."),
          ## runGlm (twoGroups.R:344-360); ezMethodEdger numeric grouping2 (app-edgerTwoGroups.R:22-25)
          if (glm) paste0("Samples outside the sample group, the reference group and any sampleGroupBaseline or refGroupBaseline were removed, and the model used the no-intercept design ",
                          if (g2) "~0 + group + grouping2 (additive, no interaction), with a grouping2 made only of numbers entered as a continuous covariate." else "~0 + group."),
          ## runGlm contrastsIndices (twoGroups.R:376-380, 385-389); prior.count = backgroundExpression (twoGroups.R:117, 128)
          paste0("Log2 fold changes are for the sample group over the reference group",
                 if (glm) paste0(" (the glm contrast ",
                                 if (baselines) "(sample group minus sampleGroupBaseline) minus (reference group minus refGroupBaseline)" else "sample group minus reference group",
                                 ")"),
                 ", computed by edgeR with a prior count equal to backgroundExpression (prior.count)."),
          ## runGlm deTest == "QL" (twoGroups.R:373-381); glmQLFit.DGEList legacy = FALSE, dispersion = NULL
          if (glm && identical(param$deTest, "QL")) "Genes were tested with the quasi-likelihood F-test (glmQLFit, glmQLFTest) in edgeR's default non-legacy mode, in which glmQLFit estimates its own negative-binomial dispersion from the most abundant genes, so the estimateDisp dispersions were not used by this test.",
          ## runGlm deTest == "LR" (twoGroups.R:363-368, 382-391); estimateDisp defaults trend.method = "locfit", tagwise = TRUE
          if (glm && identical(param$deTest, "LR")) paste0(
            "Dispersions were estimated with estimateDisp on the design matrix (common, locfit-trended and tagwise) and genes were tested with glmFit and glmLRT using the tagwise dispersions; robust dispersion estimation (estimateGLMRobustDisp) was ",
            if (robust) "used (robust true)." else "not used (robust false)."),
          ## runEdger (twoGroups.R:283-297); exactTest dispersion = "auto"
          if (exact) "Dispersions were estimated with estimateDisp over all groups of the dataset and the reference and sample groups were compared with the edgeR exact test using the tagwise dispersions.",
          ## runEdger (twoGroups.R:288-292) / runGlm (twoGroups.R:363-371)
          "When the sample group or the reference group had fewer than two samples, the negative-binomial dispersion was fixed at 0.1 instead of being estimated (used by the exact test and the likelihood-ratio test; the QL fit estimates its own).",
          ## compileEnrichmentInput / ezEnricher / ezGSEA (go-analysis.R:153-157, 568-600) -- shared with DESeq2
          if (runGO) "When GO annotation is available, genes with p-value at or below pValThreshGO and absolute log2 ratio above log2RatioThreshGO were tested for GO over-representation (BP, MF, CC; up-regulated, down-regulated and both separately) with clusterProfiler enricher, using the present genes as universe, Benjamini-Hochberg adjustment, gene sets of 10 to 500 genes, cut-off fdrThreshORA and at least 3 genes per term.",
          ## ezGSEA (go-analysis.R:617-660) -- shared with DESeq2
          if (runGO) "GO gene set enrichment analysis (BP, MF, CC) was run with clusterProfiler GSEA on the present genes ranked by rankMetric (log2 ratio, -log10 p-value, or signed -log10 p-value), with Benjamini-Hochberg adjustment, gene sets of 10 to 500 genes and cut-off fdrThreshGSEA."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodEdger
        name <<- "EzAppEdger"
        appDefaults <<- rbind(
          testMethod = ezFrame(
            Type = "character",
            DefaultValue = "glm",
            Description = "which test method in edgeR to use: glm or exactTest"
          ),
          normMethod = ezFrame(
            Type = "character",
            DefaultValue = "TMM",
            Description = "edgeR's norm method: TMM, upperquartile, RLE, or none"
          ),
          useRefGroupAsBaseline = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "should the log-ratios be centered at the reference samples"
          ),
          onlyCompGroupsHeatmap = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Only show the samples from comparison groups in heatmap"
          ),
          priorCount = ezFrame(
            Type = "numeric",
            DefaultValue = 10,
            Description = "prior count to be added to shrink the log-fold-changes"
          ),
          deTest = ezFrame(
            Type = "character",
            DefaultValue = "QL",
            Description = "edgeR's differential expression test method: QL or LR"
          ),
          runGfold = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "should gfold run"
          ),
          doPrecomputeEnrichr = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "should enrichr be precomputed"
          ),
          rankMetric = ezFrame(
              Type = "character",
              DefaultValue = 'log2Ratio',
              Description = "how to rank genes for GSEA"
          )
        )
      }
    )
  )
