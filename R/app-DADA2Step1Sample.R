###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodDADA2Step1Sample = function(input=NA, output=NA, param=NA, 
                                     htmlFile="00index.html"){
  
  require(dada2)
  require(purrr)
  dataset = input$meta
  sampleName = input$getNames() 
  minLen <- param$minLen
  maxLen <- param$maxLen
  isPaired <- param$paired
  ### read fastq files and prepare inputs for DADA2
  ### are reads paired? should they be joined? 
  file1PathInDataset <- input$getFullPaths("Read1")
  if(isPaired){
    file2PathInDataset <- input$getFullPaths("Read2")
  DADA2mainSeqTabObj <- DADA2CreateSeqTab(sampleName,maxLen,file1PathInDataset,
                                          file2PathInDataset)
  }else{
    DADA2mainSeqTabObj <- DADA2CreateSeqTab(sampleName,maxLen,file1PathInDataset,minLen)
  }
  
  ## rename output files
  ## Files needed for the report 
  #1) 
  DADA2mainSeqTabObjFileName <- basename(output$getColumn("RObjectWithSeqTab"))
  saveRDS(DADA2mainSeqTabObj,DADA2mainSeqTabObjFileName)
}

##' @template app-template
##' @templateVar method ezMethodDADA2Step1Sample()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppDADA2Step1Sample <-
  setRefClass("EzAppDADA2Step1Sample",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDADA2Step1Sample
                  name <<- "EzAppDADA2Step1Sample"
                  appDefaults <<- rbind(minLen = ezFrame(Type="integer",  DefaultValue="290",Description="Min length"),     
                                        maxLen= ezFrame(Type="integer",  DefaultValue="330",Description="Max length")
                  )
                }
              )
  )












