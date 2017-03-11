###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodCanu = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  opt = param$cmdOptions
  sampleName = input$getNames()
  SMRT_Path = input$getFullPaths("Reads")
  SMRT_File = basename(input$getColumn("Reads")))
  ezSystem(paste("cp -r", SMRT_Path, "."))
  ezSystem(paste("mkdir", "smrt_input"))
  ezSystem(paste("tar -zxf", SMRT_File, "--strip-components=4 -C smrt_input"))
  readFile = file.path(dirname(.), "smrt_input", "Analysis_Results", "*.subreads.fastq") 
  ezSystem(paste("cat", readFile, ">", paste0(sampleName,".fastq")))
  fixOpt = paste("useGrid=false", "gnuplotTested=true")
  cmd = paste(CANU, "-p", sampleName, "-d", sampleName, paste0("genomeSize=", param$canuGenomeSize, "m"), fixOpt, opt, 
                param$canuReadOpt, paste0(sampleName,".fastq"), "1> ", paste0(sampleName,"_canu.log"))
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppCanu <-
  setRefClass("EzAppCanu",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCanu
                  name <<- "EzAppCanu"
                  appDefaults <<- rbind(canuReadOpt = ezFrame(Type="character",  DefaultValue="-pacbio-raw",  Description="input read types: -pacbio-raw, -pacbio-corrected, -nanopore-raw, -nanopore-corrected. Default is pacbio raw data"),
                                        canuGenomeSize = ezFrame(Type="integer",  DefaultValue="5",  Description="estimated genome size in Mb"))
                }
              )
)
              
              
