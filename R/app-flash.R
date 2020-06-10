###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodFlash = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  sampleName = input$getNames()
  stopifnot((param$paired))
  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  cmd = paste("flash",trimmedInput$getColumn("Read1"), trimmedInput$getColumn("Read2"),
              "-o",sampleName,'-t',ezThreads(),opt,"1> ",paste0(sampleName,"_flash.log"))
  ezSystem(cmd)
  cmd = paste0('pigz ',sampleName,'.extendedFrags.fastq')
  ezSystem(cmd)
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodFlash(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppFlash <-
  setRefClass("EzAppFlash",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodFlash
                  name <<- "EzAppFlash"
                }
              )
  )
