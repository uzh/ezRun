###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title The object that represents a result of an expression analysis
##' @field param \code{EzParam} object that was used during the result generation
##' @field rawData list with the rawData that was used during the result generation
##' @field result list holding the actual results, in particular, the log2 ration and the p-value
##' @section Functions:
##' \itemize{
##'   \item{\code{saveToFile(file): }}{Saves the object to the given filepath. Saved files can be loaded with \code{EzResult$new(file="mystoredResult.RData"}}
##' }
##' @template roxygen-template
##' @examples
##' deResult = EzResult$new()
##' deResult = EzResult$new(param=list(p1="p1", p2="p2"), rawData=list(a=0, b=1), result=list(u=0, v=1))
##' rdFile = tempfile(fileext = ".RData")
##' deResult$saveToFile(file = rdFile)
##' deResult2 = EzResult$new(file=rdFile)
EzResult <-
  setRefClass("EzResult",
              fields=c("param", "rawData", "result", "se", "sceset"),
              methods=list(
                initialize = function(paramNew=list(), rawDataNew=list(), 
                                      resultNew=list(),
                                      file=NULL, seNew=SummarizedExperiment::SummarizedExperiment(), 
                                      scesetNew=list()){
                  param <<- paramNew
                  rawData <<- rawDataNew
                  result <<- resultNew
                  sceset <<- scesetNew
                  se <<- seNew
                  if (!is.null(file)){
                    stopifnot(length(paramNew) == 0 && length(rawDataNew) == 0 && length(resultNew) == 0)
                    stopifnot(file.exists(file))
                    if (grepl("RData$", file)){
                    load(file = file) ## loads
                    param <<- param
                    rawData <<- rawData
                    result <<- result
                    se <<- se
                    sceset <<- sceset
                    }
                    if (grepl("rds$", file)){
                      se <<- readRDS(file)
                    }
                  }
                },
                saveToFile = function(file){
                  save(param, rawData, result, se, sceset, file=file)
                }
              )
  )

makeSummarizedExperiment = function(param, rawData, result){
  require(SummarizedExperiment)
  if (is.null(rawData)){
    return(NULL)
  }
  assayList = list(counts=rawData$counts)
  if (!is.null(rawData$signal)){
    assayList$countsNorm = rawData$signal
  } else {
    if (!is.null(result$xNorm)){
      assayList$countsNorm = result$xNorm
    }
  }

  SummarizedExperiment(assays=assayList,
                       rowData=rawData$seqAnno,
                       colData=ezDesignFromDataset(rawData$dataset, param),
                       metadata=param)
  ## DE results can be represented
  ## --- as mcols in the rowData
  ## --- as an element in the metaData list (subsetting does not work)
  ## --- as a separate SummarizedExperiment object
}
