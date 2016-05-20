###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodTrinity = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  
  if (ezIsSpecified(param$samples)){
    input$subset(param$samples)
  }
  
  trimmedInput = ezMethodTrim(input = input, param = param)
  param$dataRoot = ""
        
  if (param$paired){
    read1 = paste(trimmedInput$getColumn("Read1"), collapse=",")
    read2 = paste(trimmedInput$getColumn("Read2"), collapse=",")
    readOpt = paste("--left", read1, "--right", read2)
#    readOpt = paste("--left", reads1, "--right", reads2)
    libOpt = switch(param$strandMode, sense="--SS_lib_type FR", antisense="--SS_lib_type RF", both="")
  } else {
    read1 = paste(trimmedInput$getColumn("Read1"), collapse=",")
    readOpt = paste("--single", read1)
#    readOpt = paste("--single", reads1)
    libOpt = switch(param$strandMode, sense="--SS_lib_type F", antisense="--SS_lib_type R", both="")
  }
  
  cmd = paste(TRINITY, "-seqType fq", readOpt,
              "--max_memory", paste0(param$ram, "G"), "--bflyCalculateCPU", ## "--bflyHeapSpaceMax", paste0(round(as.numeric(param$ram)/4), "G"),
              "--CPU", ezThreads(),
              libOpt,
              param$trinityOpt,
              "--output", "trinity", ">", "trinity.stdout")
  ezSystem(cmd)
  ezSystem(paste("mv", "trinity/Trinity.fasta", basename(output$getColumn("Fasta"))))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodTrinity(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
##' @seealso \code{\link{ezMethodTrim}}
EzAppTrinity <-
  setRefClass("EzAppTrinity",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodTrinity
                  name <<- "EzAppTrinity"
                  appDefaults <<- rbind(trinityOpt = ezFrame(Type="character",  DefaultValue="--min_kmer_cov 2",  Description="trinity commandline options"))
                }
              )
  )
