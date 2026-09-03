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
      ## GO enrichment (runGO) actually runs via goseq + GO.db -- go-analysis.R only
      ## require()'s GOstats, never calls it; the earlier citation named the wrong
      ## tool. RUVSeq gated on runRUV. pheatmap/WGCNA-dendrogram/cluster-validation
      ## stats (gap/silhouette/DB-index) deliberately excluded: no citable paper, or
      ## only a plotting helper is used rather than the method itself.
      citation = function() {
        c(
          "Rehrauer, H. et al. ezRun: An R meta-package for the analysis of Next Generation Sequencing Data. https://github.com/uzh/ezRun",
          "Young, M.D., Wakefield, M.J., Smyth, G.K. & Oshlack, A. Gene ontology analysis for RNA-seq: accounting for selection bias. Genome Biology 11, R14 (2010). https://doi.org/10.1186/gb-2010-11-2-r14",
          "Bioconductor. GO.db: A set of annotation maps describing the entire Gene Ontology. R package version 3.23.1. https://doi.org/10.18129/B9.bioc.GO.db",
          "Risso, D., Ngai, J., Speed, T.P. & Dudoit, S. Normalization of RNA-seq data using factor analysis of control genes or samples. Nature Biotechnology 32(9), 896-902 (2014). https://doi.org/10.1038/nbt.2931"
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
