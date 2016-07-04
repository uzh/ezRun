###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodDeseq2 = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  if (!is.null(param$markOutliers) && param$markOutliers){
    stop("DESeq2 does not support marking outliers because marked outliers would still be used in dispersion estimates")
  }
  setwdNew(basename(output$getColumn("Report")))
  stopifnot(param$sampleGroup != param$refGroup)
  
  input = cleanupTwoGroupsInput(input, param)
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1){
    param$grouping2 = input$getColumn(param$grouping2)
  }
  
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport(htmlFile, param=param, error=rawData$error)
    return("Error")
  }
  
  deResult = twoGroupCountComparison(rawData, param)
  if (isError(deResult)){
    writeErrorReport(htmlFile, param=param, error=deResult$error)
    return("Error")
  }
  
  writeNgsTwoGroupReport(rawData$dataset, deResult, output, htmlFile)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodDeseq2(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
EzAppDeseq2 <-
  setRefClass("EzAppDeseq2",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDeseq2
                  name <<- "EzAppDeseq2"
                  appDefaults <<- rbind(testMethod=ezFrame(Type="character",  DefaultValue="deseq2",  Description="which test method in DESeq to use: deseq2"),
                                        runGfold=ezFrame(Type="logical", DefaultValue=FALSE, Description="no need to compute moderated ratios; deseq2 does this already"))
                }
              )
  )
