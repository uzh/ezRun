###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMegahit = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  ### run UPARSE and reformat output for phylotseq 
  
  library(plyr)
  library(dplyr)
  
  sampleName = input$getNames()
  file1PathInDatset <- input$getFullPaths("Read1")
  fastqName1 <- paste0(sampleName,".R1.fastq")
  fastqName1Trimmed  <- paste0(sampleName,".R1.trimmed.fastq") 
  cpCmd1 <- paste0("gunzip -c ", file1PathInDatset, "  > ", fastqName1)
  ezSystem(cpCmd1)  
  if(param$paired){
    file2PathInDatset <- input$getFullPaths("Read2")
    fastqName2 <- paste0(sampleName,".R2.fastq")
    fastqName2Trimmed  <- paste0(sampleName,".R2.trimmed.fastq")
    cpCmd2 <- paste0("gunzip -c ", file2PathInDatset, "  > ", fastqName2)
    ezSystem(cpCmd2)
    inputString <- paste0("-1", fastqName1Trimmed, "-2", fastqName2Trimmed)
    pairedString="YES"
  } else {
    inputString <- paste("-r",fastqName1Trimmed)
    pairedString="NO"
  }
  megahitTemplScript <- MEGAHIT_TEMPLATE_SCRIPT
  megahitToBeExec <- paste0("megahit.",sampleName,".sh")
  ##update template
  updateTemplateScriptCmd <- paste0("sed -e s/\"SAMPLE_NAME\"/", sampleName, "/g",
                                    " -e s/\"INPUT_FILE_STRING\"/", inputString, "/g ",
                                    " -e s/\"KMER_LIST\"/", param$kmerList, "/g ",
                                    " -e s/\"ARE_READ_PAIRED\"/", pairedString, "/g ",
                                    megahitTemplScript, " >",
                                    megahitToBeExec)
  ezSystem(updateTemplateScriptCmd)
  ## run script
  megahitCmd <- paste("bash",megahitToBeExec)
  ezSystem(megahitCmd)
  
  ## rename output files
  #1) 
  oldContigFile <- "megahitResults/final.contigs.fa"
  newContigFile <- basename(output$getColumn("contigFile"))
  ezSystem(paste("mv",oldContigFile,newContigFile))
  #2) 
  oldProdigalFile <- "prodigalAnnotation.gff"
  newProdigalFile <- basename(output$getColumn("prodigalPredictionFile"))
  ezSystem(paste("mv",oldProdigalFile,newProdigalFile))
  #3) 
  oldDiamAnnFile <- "annotatedProteins.tsv"
  newDiamAnnFile <- basename(output$getColumn("diamondAnnotationFile"))
  ezSystem(paste("mv",oldDiamAnnFile,newDiamAnnFile))
}

##' @template app-template
##' @templateVar method ezMethodMegahit()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMegahit <-
  setRefClass("EzAppMegahit",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMegahit
                  name <<- "EzAppMegahit"
                  appDefaults <<- rbind(fastqErrorMax = ezFrame(Type="numeric",
                                                                DefaultValue="1",
                                                                Description="Max EE in fastx_truncate")
                                        
                  )
                }
              )
  )

