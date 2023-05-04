###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodPreqc = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  sampleName = input$getNames()
  cores=param$cores
  ##stopifnot((param$paired))
  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  if (param$paired){
    read1 = trimmedInput$getColumn("Read1")
    read2 = trimmedInput$getColumn("Read2")
    readOpt = paste(read1, read2)
    cmd = paste("sga preprocess --pe-mode", param$peMode, "--dust", paste0("--dust-threshold=",param$dustThreshold), readOpt, ">", paste0(sampleName,".fastq"))
    ezSystem(cmd)
  } else {
    read1 = trimmedInput$getColumn("Read1")
    readOpt = paste(read1)
    cmd = paste("sga preprocess --pe-mode 0", "--dust", paste0("--dust-threshold=",param$dustThreshold), readOpt, ">", paste0(sampleName,".fastq"))
    ezSystem(cmd)
  }
  cmd=paste("sga index -a ropebwt --no-reverse -t", cores, paste0(sampleName, ".fastq")) 
  ezSystem(cmd)
  cmd=paste("sga preqc -t", cores, paste0(sampleName, ".fastq"), ">" , paste0(sampleName, ".preqc"))
  ezSystem(cmd)
  cmd=paste("sga-preqc-report.py -o", sampleName, paste0(sampleName,  ".preqc"), paste0("/usr/local/ngseq/src/sga/src/examples/preqc/", param$example, ".preqc"))
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppPreqc <-
  setRefClass("EzAppPreqc",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodPreqc
                  name <<- "EzAppPreqc"
                  appDefaults <<- rbind(peMode = ezFrame(Type="integer",  DefaultValue="1",  Description="mode of paired-end reads"),
                                        dustThreshold = ezFrame(Type="numeric",  DefaultValue="4.0",  Description="filter out reads that have a dust score higher than 4.0"),
		  			example = ezFrame(Type="character", DefaultValue="human", Description="pre-computed readset"))
                }
              )
)


