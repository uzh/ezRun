###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Deseq2
##' @seealso \code{\link{EzAppDeseq2}}
ezMethodDeseq2 = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  param$testMethod = "deseq2"
  param$normMethod = ""
  param$runGfold = FALSE ## no need to compute moderated ratios; deseq2 does this already
  if (!is.null(param$markOutliers) && param$markOutliers){
    stop("DESeq2 does not support marking outliers because marked outliers would still be used in dispersion estimates")
  }
  cwd = getwd()
  on.exit(setwd(cwd))
  setwdNew(basename(output$getColumn("Report")))
  stopifnot(param$sampleGroup != param$refGroup)
  
  input = cleanupTwoGroupsInput(input, param)
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$batch) && length(param$batch) == 1){
    param$batch = input$meta[[param$batch]]
  }
  
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport(htmlFile, param=param, error=rawData$error)
    return("Error")
  }
  
  result = twoGroupCountComparison(rawData, param)
  if (isError(result)){
    writeErrorReport(htmlFile, param=param, error=rawData$error)
    return("Error")
  }
  result$featureLevel = rawData$featureLevel
  result$countName = rawData$countName
  
  writeNgsTwoGroupReport(input$meta, result, htmlFile, param=param, rawData=rawData)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodDeseq2()
##' @seealso \code{\link{ezMethodDeseq2}}
EzAppDeseq2 <-
  setRefClass("EzAppDeseq2",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDeseq2
                  name <<- "EzAppDeseq2"
                }
              )
  )
