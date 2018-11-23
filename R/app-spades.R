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
  if (param$gneomeType == "pacbioSmrtCell"){
    basicOpts <- paste(param$spadesBasicOpt, "--meta")
  } else {
    basicOpts <- param$spadesBasicOpt
  }
  if (param$paired){
    read1 = trimmedInput$getColumn("Read1")
    read2 = trimmedInput$getColumn("Read2")
    readOpt = paste("-1", read1, "-2", read2)
    cmd = paste("spades.py", readOpt, basicOpts,
                param$spadesPipeOpt, "-m", param$ram, "-o", "spades", '-t', ezThreads(), opt, "1> ", paste0(sampleName,"_spades.log"))
    ezSystem(cmd)
  } else {
    read1 = trimmedInput$getColumn("Read1")
    readOpt = paste("-s", read1)
    cmd = paste("spades.py", readOpt, basicOpts,
                param$spadesPipeOpt, "-m", param$ram, "-o", "spades", '-t', ezThreads(), opt, "1> ", paste0(sampleName,"_spades.log"))
    ezSystem(cmd)
  }
  wddir <- "."
  sfile <- file.path(wddir, "spades/scaffolds.fasta")
  cfile <- file.path(wddir, "spades/contigs.fasta")
  if (file.exists(sfile)){
  	ezSystem(paste("cp", "spades/scaffolds.fasta", basename(output$getColumn("Draft"))))
  }else if (file.exists(cfile)){
  	ezSystem(paste("cp", "spades/contigs.fasta", basename(output$getColumn("Draft"))))
  }
  ezSystem(paste("mv", "spades", sampleName))
  #cmd = paste0('mv ', sampleName, '_preprocessing.log ', sampleName, '_trimmomatic.log')
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
              
              
