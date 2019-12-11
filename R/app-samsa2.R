###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodSamsa2 = function(input=NA, output=NA, param=NA, 
                              htmlFile="00index.html"){
  ### metatranscriptomics assemby with Samsa2, annotation with RefSeq
  
  library(plyr)
  library(dplyr)
  library(Biostrings)
  
  sampleName = input$getNames()
  file1PathInDatset <- input$getFullPaths("Read1")
  fastqName1 <- paste0(sampleName,".R1.fastq")
  cpCmd1 <- paste0("gunzip -c ", file1PathInDatset, "  > ", fastqName1)
  ezSystem(cpCmd1)  
  if(param$paired){
    file2PathInDatset <- input$getFullPaths("Read2")
    fastqName2 <- paste0(sampleName,".R2.fastq")
    cpCmd2 <- paste0("gunzip -c ", file2PathInDatset, "  > ", fastqName2)
    ezSystem(cpCmd2)
  } 
  ##update template
  updateTemplateScriptCmd <- paste("cp",
                                   file.path(METAGENOMICS_ROOT,SAMSA2_TEMPLATE_SCRIPT),
                                    SAMSA2_TEMPLATE_SCRIPT)
  ezSystem(updateTemplateScriptCmd)
  ## run script
  metaspadesToBeExecCmd <- paste("bash",SAMSA2_TEMPLATE_SCRIPT)
  ezSystem(metaspadesToBeExecCmd)
  
  ## place output files

  oldAnnFile <- list.files("step_4_output",pattern = "RefSeq_annotated",full.names = T)
  newAnnFile <- basename(output$getColumn("annotationFileRefSeq"))
  ezSystem(paste("mv",oldAnnFile,newAnnFile))
}

##' @template app-template
##' @templateVar method ezMethodSamsa2()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppSamsa2 <-
  setRefClass("EzAppSamsa2",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSamsa2
                  name <<- "EzAppSamsa2"
                  appDefaults <<- rbind(useSubsystemDB = ezFrame(Type="logical",
                                                                     DefaultValue="RefSeq",
                                                                     Description="database")
                                        
                                        
                  )
                }
              )
  )

