###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title The object that represents a result of an expression analysis
##' @description 
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
              fields=c("param", "rawData", "result"),
              methods=list(
                initialize = function(paramNew=list(), rawDataNew=list(), resultNew=list(),
                                      file=NULL){
                  param <<- paramNew
                  rawData <<- rawDataNew
                  result <<- resultNew
                  if (!is.null(file)){
                    stopifnot(length(paramNew) == 0 && length(rawDataNew) == 0 && length(resultNew) == 0)
                    stopifnot(file.exists(file))
                    load(file = file) ## loads
                    param <<- param
                    rawData <<- rawData
                    result <<- result
                  }
                },
                saveToFile = function(file){
                  save(param, rawData, result, file=file)
                }
              )
  )
