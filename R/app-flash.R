###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodFlash = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  sampleName = input$getNames()
  param$fastpCompression = 9
  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  if(!param$skipFlash){
    stopifnot((param$paired))
    cmd = paste("flash",trimmedInput$getColumn("Read1"), trimmedInput$getColumn("Read2"),
              "-o",sampleName,'-t',ezThreads(),opt,"1>> ",paste0(sampleName,"_preprocessing.log"))
    ezSystem(cmd)
    cmd = paste0('pigz --best ',sampleName,'.extendedFrags.fastq')
    ezSystem(cmd)
    cmd = paste('mv', paste0(sampleName,'.extendedFrags.fastq.gz'), paste0(sampleName,'.R1.fastq.gz'))
    ezSystem(cmd)
  } else {
      ezSystem(paste('mv', paste0(sampleName,'-trimmed_R1.fastq.gz'), paste0(sampleName,'.R1.fastq.gz')))
      if(param$paired){
        ezSystem(paste('mv',paste0(sampleName, '-trimmed_R2.fastq.gz'), paste0(sampleName,'.R2.fastq.gz')))
      }
  }
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
                  appDefaults <<- rbind(
                      skipFlash = ezFrame(
                          Type = "logical",
                          DefaultValue = FALSE,
                          Description = "run or skip flash"
                      )
                  )
                }
              )
  )
