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
##' @templateVar method ezMethodEdger(input=NA, output=NA, param=NA, htmlFile="00index.html")
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
                                        normMethod=ezFrame(Type="character", DefaultValue="TMM", Description="edgeR's norm method: TMM, upperquartile, RLE, or none"),
                                        useRefGroupAsBaseline=ezFrame(Type="logical", DefaultValue=FALSE, Description="should the log-ratios be centered at the reference samples"))
                  
                }
              )
  )
