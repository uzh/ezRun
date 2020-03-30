###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodCountSpacer = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  sampleName = input$getNames()
  trimmedInput = ezMethodTrim(input = input, param = param)
  cmd = paste("python", "/usr/local/ngseq/src/countSpacer/count_spacers_LO.py", "-f", trimmedInput$getColumn("Read1"), "-o", 
              paste0(sampleName, "_counts.csv"), "-i", param[['dictPath']], "--keyStart", as.numeric(param[['keyStart']]), "--keyEnd", 
              as.numeric(param[['keyEnd']]), opt)
  ezSystem(cmd)
  
  ezSystem(paste('mv statistics.txt', paste0(sampleName,'_statistics.txt')))
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodFlash(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppCountSpacer <-
  setRefClass("EzAppCountSpacer",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCountSpacer
                  name <<- "EzAppCountSpacer"
                }
              )
  )
