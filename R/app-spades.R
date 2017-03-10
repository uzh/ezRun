###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodSpades = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  opt = param$cmdOptions
  sampleName = input$getNames()
  ##stopifnot((param$paired))
  trimmedInput = ezMethodTrim(input = input, param = param)
  if (param$paired){
    read1 = trimmedInput$getColumn("Read1")
    read2 = trimmedInput$getColumn("Read2")
    readOpt = paste("-1", read1, "-2", read2)
    cmd = paste(SPADES, readOpt, param$spadesBasicOpt,
                param$spadesPipeOpt, "-m", param$ram, "-o", "spades", '-t', ezThreads(), opt, "1> ", paste0(sampleName,"_spades.log"))
    ezSystem(cmd)
  } else {
    read1 = trimmedInput$getColumn("Read1")
    readOpt = paste("-s", read1)
    cmd = paste(SPADES, readOpt, param$spadesBasicOpt,
                param$spadesPipeOpt, "-m", param$ram, "-o", "spades", '-t', ezThreads(), opt, "1> ", paste0(sampleName,"_spades.log"))
    ezSystem(cmd)
  }
  #pathFasta = file.path(param$outputDir, sampleName, "scaffolds.fasta")
  #ezSystem(paste("mv", pathFasta, basename(output$getColumn("Fasta"))))
  ezSystem(paste("mv", "spades/scaffolds.fasta", basename(output$getColumn("Fasta"))))
  cmd = paste0('mv trimmomatic.err ', sampleName, '_trimmomatic.log')
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppSpades <-
  setRefClass("EzAppSpades",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSpades
                  name <<- "EzAppSpades"
                  appDefaults <<- rbind(spadesBasicOpt = ezFrame(Type="character",  DefaultValue="",  Description="spades basic options: --sc --meta --rna --plasmid. Default is empty for genome assembly without MDA"),
                                        spadesPipeOpt = ezFrame(Type="character",  DefaultValue="--careful",  Description="spades pipeline options: --only-assembler --careful"))
                }
              )
)
              
              
