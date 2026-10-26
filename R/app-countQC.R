###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCountQC = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  setwdNew(basename(output$getColumn("Report")))
  dataset <- input$meta
  if (param$useFactorsAsSampleName) {
    dataset$Name = rownames(dataset)
    rownames(dataset) = addReplicate(apply(
      ezDesignFromDataset(dataset),
      1,
      paste,
      collapse = "_"
    ))
  }
  input$meta <- dataset
  rawData <- loadCountDataset(input, param)

  if (isError(rawData)) {
    writeErrorReport(htmlFile, param = param, error = rawData$error)
    return("Error")
  }

  ## The analysis below runs inside the CountQC.qmd render (a quarto subprocess),
  ## whose messages do not reach the job log, so its methods are logged here.
  if (ncol(rawData) >= 2) {
    ezLogMethods(
      "A feature was called present in a sample when its count exceeded ",
      if (isTRUE(param$useSigThresh)) param$sigThresh else "0",
      "."
    )
    if (ezIsSpecified(param$normMethod) && length(param$normMethod) == 1) {
      ezLogMethods(
        switch(
          param$normMethod,
          none = "The counts were not normalized",
          logMean = "The counts were scaled per sample by the geometric mean of the features present in all samples",
          median = "The counts were scaled per sample by the median of the features present in all samples",
          quantile = "The counts were quantile normalized (limma normalizeQuantiles)",
          vsn = "The counts were normalized with vsn",
          paste0("The counts were normalized with the ", param$normMethod, " method")
        ),
        if (ezIsSpecified(param$minSignal)) paste0(", and a pseudo-count of ", param$minSignal, " was added"),
        "."
      )
    }
    ezLogMethods(
      "Sample correlation, clustering and MDS used the features present in at least half of the samples ",
      "of at least half of the conditions",
      if (!is.null(rowData(rawData)$IsControl)) ", excluding control sequences",
      ", on log2(normalized signal + ", param$backgroundExpression, ")."
    )
    ezLogMethods(
      "The top-gene analyses used the ", param$topGeneSize, " features with ",
      if (isTRUE(param$selectByFtest)) {
        "the smallest ANOVA F-test p-value across the conditions"
      } else {
        "the highest standard deviation of the log2 signal across samples"
      },
      "."
    )
    ezLogMethods(
      "Sample correlations were Pearson correlations of the log2 signal",
      if (ncol(rawData) > 3) "; samples were clustered hierarchically (Ward.D2) on 1 minus the Pearson correlation",
      "."
    )
    if (ncol(rawData) > 3) {
      ezLogMethods(
        "Up to ", param$maxGenesForClustering, " features with the highest standard deviation of the log2 signal (above ",
        param$highVarThreshold, ") were row-centred and clustered hierarchically (Euclidean distance, Ward.D2) into ",
        param$nSampleClusters, " clusters."
      )
      ezLogMethods(
        "The number of feature clusters was assessed for 2 to 10 clusters with the Dunn, Calinski-Harabasz, ",
        "silhouette, entropy, Davies-Bouldin and gap statistics."
      )
      if (isTRUE(param$runGO)) {
        ezLogMethods(
          "Each feature cluster was tested for GO over-representation (BP, MF, CC) with goseq (hypergeometric, ",
          "without length-bias correction) against all annotated genes of the count table",
          if (isTRUE(param$includeGoParentAnnotation)) ", with each gene also annotated to the parent terms of its GO terms",
          "; terms with an unadjusted p-value below ", param$pValThreshFisher, " and at least ",
          param$minCountFisher, " cluster genes were reported."
        )
      }
      ezLogMethods(
        "MDS was computed with limma plotMDS (default settings) on the log2 signal."
      )
    }
  }

  metadata(rawData)$output <- output
  makeQuartoReport(
    rawData = rawData,
    qmdFile = "CountQC.qmd",
    reportTitle = "CountQC",
    colour = isTRUE(param$colour),
    number = isTRUE(param$number)
  )

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodCountQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
EzAppCountQC <-
  setRefClass(
    "EzAppCountQC",
    contains = "EzApp",
    methods = list(
      ## ezRun unconditional; goseq (goClusterResults -> ezGoseq, hypergeometric test, in the
      ## report) + GO.db unless runGO is false (default TRUE; doGo also needs GO annotation,
      ## data-dependent); RUVSeq only when runRUV is true (loadCountDataset; no Ruby app
      ## declares it). pheatmap/WGCNA-dendrogram/cluster-validation stats deliberately
      ## excluded: no citable paper, or only a plotting helper is used rather than the method.
      citation = function(param = list()) {
        runGO <- !isFALSE(as.logical(param$runGO))
        c(
          "Rehrauer, H. et al. ezRun: An R meta-package for the analysis of Next Generation Sequencing Data. https://github.com/uzh/ezRun",
          if (runGO) "Young, M.D., Wakefield, M.J., Smyth, G.K. & Oshlack, A. Gene ontology analysis for RNA-seq: accounting for selection bias. Genome Biology 11, R14 (2010). https://doi.org/10.1186/gb-2010-11-2-r14",
          if (runGO) "Carlson, M. GO.db: A set of annotation maps describing the entire Gene Ontology. R package. https://doi.org/10.18129/B9.bioc.GO.db",
          if (isTRUE(as.logical(param$runRUV))) "Risso, D., Ngai, J., Speed, T.P. & Dudoit, S. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature Biotechnology 32(9), 896-902 (2014). https://doi.org/10.1038/nbt.2931"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodCountQC
        name <<- "EzAppCountQC"
        appDefaults <<- rbind(
          runGO = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "whether to run the GO analysis"
          ),
          nSampleClusters = ezFrame(
            Type = "numeric",
            DefaultValue = 6,
            Description = "Number of SampleClusters, default value 6"
          ),
          selectByFtest = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "select topGenes by Test instead of SD"
          ),
          topGeneSize = ezFrame(
            Type = "numeric",
            DefaultValue = 100,
            Description = "number of genes to consider in gene clustering, mds etc"
          ),
          colour = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "coloured tab hierarchy in the report"
          ),
          number = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "numbered tab hierarchy in the report"
          )
        )
      }
    )
  )
