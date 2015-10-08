###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @author Opitz, Lennart
##' @template method-template
##' @templateVar methodName MEME
##' @seealso \code{\link{EzAppMEME}}
ezMethodMEME = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  db = param$motifDB
  sampleName = input$getNames()
  cmd = paste(MEME,"-oc",sampleName,"-index-name",paste(sampleName,"_meme-chip.html",sep=""),"-time 300 -order 1",db, opt, 
              input$getFullPaths(param, "PeakSequences"))
  ezSystem(cmd)
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodMEME()
##' @seealso \code{\link{ezMethodMEME}}
EzAppMEME <-
  setRefClass("EzAppMEME",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  runMethod <<- ezMethodMEME
                  name <<- "EzAppMEME"
                }
              )
  )
