###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMEME = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  db = param$motifDB
  sampleName = input$getNames()
  cmd = paste(MEME,"-oc",sampleName,"-index-name",paste0(sampleName,"_meme-chip.html"),"-time 300 -order 1",db, opt, 
              input$getFullPaths(param, "PeakSequences"))
  ezSystem(cmd)
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodMEME
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMEME <-
  setRefClass("EzAppMEME",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMEME
                  name <<- "EzAppMEME"
                }
              )
  )
