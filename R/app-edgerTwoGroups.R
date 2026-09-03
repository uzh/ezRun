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
          "Bioconductor. GO.db: A set of annotation maps describing the entire Gene Ontology. R package version 3.23.1. https://doi.org/10.18129/B9.bioc.GO.db",
          "Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128",
          "Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377"
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
