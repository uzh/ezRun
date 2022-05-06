###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNanoPlot = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  opt = param$cmdOptions
  sampleName = input$getNames()
  cmd = paste("NanoPlot", "-t", param$cores, "-p", paste0(sampleName, "."), "--title", sampleName, "-o", NanoPlot_Result, "--fastq", input$getFullPaths("Read1"))
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodNanoPlot()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppNanoPlot <-
  setRefClass("EzAppNanoPlot",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodNanoPlot
                  name <<- "EzAppNanoPlot"
                }
              )
)
              
              
