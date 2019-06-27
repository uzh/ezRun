###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodKraken = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  opt = param$cmdOptions
  sampleName = input$getNames()
  trimmedInput = ezMethodTrim(input = input, param = param)
  dbOpt <- param$krakenDBOpt
  conOpt <- param$krakenConfidenceOpt
  phredOpt <- param$krakenPhredOpt
  if (param$paired){
    read1 = trimmedInput$getColumn("Read1")
    read2 = trimmedInput$getColumn("Read2")
    readOpt = paste(read1, read2)
    cmd = paste("/usr/local/ngseq/stow/kraken2-2.0.7-beta/bin/kraken2 -db", paste0("/srv/GT/databases/kraken2/",dbOpt), "--paired", "--confidence", conOpt, "--minimum-base-quality", phredOpt,
           "--output", paste0(sampleName,".txt"), "--report", paste0(sampleName,".report.txt"), "--threads", ezThreads(), param$cmdOptions, readOpt, "1> ", paste0(sampleName,".log"))
    ezSystem(cmd)
    cmd = paste("/usr/local/ngseq/stow/KronaTools-2.7/bin/ktImportTaxonomy -q 2 -t 3", paste0(sampleName,".txt"), "-o", paste0(sampleName,".html"))
    ezSystem(cmd)
  } else {
    read1 = trimmedInput$getColumn("Read1")
    readOpt = paste(read1)
    cmd = paste("/usr/local/ngseq/stow/kraken2-2.0.7-beta/bin/kraken2 -db", paste0("/srv/GT/databases/kraken2/",dbOpt), "--confidence", conOpt, "--minimum-base-quality", phredOpt,
           "--output", paste0(sampleName,".txt"), "--report", paste0(sampleName,".report.txt"), "--threads", ezThreads(), param$cmdOptions, readOpt, "1> ", paste0(sampleName,".log"))
    ezSystem(cmd)
    cmd = paste("/usr/local/ngseq/stow/KronaTools-2.7/bin/ktImportTaxonomy -q 2 -t 3", paste0(sampleName,".txt"), "-o", paste0(sampleName,".html"))
    ezSystem(cmd)
  }
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppKraken <-
  setRefClass("EzAppKraken",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodKraken
                  name <<- "EzAppKraken"
                  appDefaults <<- rbind(krakenDBOpt = ezFrame(Type="character",  DefaultValue="bacteria",  Description="kraken database options: viruses bacteria. Default is bacteria"),
                                        krakenConfidenceOpt = ezFrame(Type="numeric",  DefaultValue="0.0",  Description="Confidence score threshold (default: 0.0); must be in [0, 1]."),
                                        krakenPhredOpt = ezFrame(Type="integer", DefaultValue="0",  Description="minimum Phred quality, default 0;"))
                }
              )
)
              
              
