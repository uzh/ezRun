###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodEdger = function(input=NA, output=NA, param=NA){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))
  
  stopifnot(param$sampleGroup != param$refGroup)
  
  input = cleanupTwoGroupsInput(input, param)
  param$groupingName = param$grouping
  param$grouping = input$getColumn(param$grouping)
  if (ezIsSpecified(param$grouping2) && length(param$grouping2) == 1){
    param$grouping2 = input$getColumn(param$grouping2)
  }
  
  rawData = loadCountDataset(input, param)
  if (isError(rawData)){
    writeErrorReport("00index.html", param=param, error=rawData$error)
    return("Error")
  }
  
  deResult = twoGroupCountComparison(rawData)
  if (isError(deResult)){
    writeErrorReport("00index.html", param=param, error=deResult$error)
    return("Error")
  }

  makeRmdReport(output=output, param=param, deResult=deResult, rmdFile="twoGroups.Rmd")

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
                                        useRefGroupAsBaseline=ezFrame(Type="logical", DefaultValue=FALSE, Description="should the log-ratios be centered at the reference samples"),
                                        onlyCompGroupsHeatmap=ezFrame(Type="logical", DefaultValue=FALSE, Description="Only show the samples from comparison groups in heatmap"),
                                        priorCount=ezFrame(Type="numeric", DefaultValue=10, Description="prior count to be added to shrink the log-fold-changes"),
                                        deTest=ezFrame(Type="character", DefaultValue="QL", Description="edgeR's differential expression test method: QL or LR"),
                                        runGfold=ezFrame(Type="logical", DefaultValue=FALSE, Description="should gfold run"),
					                              doPrecomputeEnrichr=ezFrame(Type="logical", DefaultValue=FALSE, Description="should enrichr be precomputed")
                                        )
                }
              )
  )
