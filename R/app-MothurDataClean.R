###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurDataClean = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  require(rmarkdown)
  mothurExe = "/usr/local/ngseq/src/mothur-1.39.5/mothur"
  mothurBatchFile = "/home/grusso/Rcodes/giancarlo/genericScipts/mothurMiSeqSOPtest.batch"
  mothurInput = "datasetNoHeaderForMothur.tsv"
  
  setwdNew(basename(output$getColumn("Report")))
  ## remove header to creat input files for mothur and copy it 
  cmdNoHead = paste("grep -v Name", "dataset.tsv", ">", mothurInput)
  ezSystem(cmdNoHead)
  cmdMothur = paste(mothurExe,mothurBatchFile)
  
}
 
##' @template app-template
##' @templateVar method ezMethodMothurDataClean()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurDataClean <-
  setRefClass("EzAppMothurDataClean",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurDataClean
                  name <<- "EzAppMothurDataClean"
                  appDefaults <<- rbind(cutOff = ezFrame(Type="real",  DefaultValue="0,03",Description="Cut-off for OTU clustering.")
                  )
                }
              )
  )
