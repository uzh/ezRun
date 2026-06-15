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
  makeRmdReport(
    rawData = rawData,
    rmdFile = "CountQC.Rmd",
    reportTitle = "CountQC",
    selfContained = TRUE
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
          )
        )
      }
    )
  )
