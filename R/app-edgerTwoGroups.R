###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodEdger = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
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
    writeErrorReport(htmlFile, param=param, error=result$error)
    return("Error")
  }
  result$featureLevel = rawData$featureLevel
  result$countName = rawData$countName
  
  writeNgsTwoGroupReport(input$meta, result, output, htmlFile, param=param, rawData=rawData)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodEdger
##' @templateVar htmlArg , htmlFile="00index.html"
##' @description Use this reference class to run a differential expression analysis with the application edgeR on two groups.
EzAppEdger <-
  setRefClass("EzAppEdger",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodEdger
                  name <<- "EzAppEdger"
                  appDefaults <<- rbind(testMethod=ezFrame(Type="character",  DefaultValue="glm",  Description="which test method in edgeR to use: glm or exactTest"),
                                        normMethod=ezFrame(Type="character", DefaultValue="TMM", Description="edgeR's norm method: TMM, upperquartile, RLE, or none"))
                }
              )
  )
