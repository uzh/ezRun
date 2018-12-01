###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMetaspades = function(input=NA, output=NA, param=NA, 
                           htmlFile="00index.html"){
  ### de novo metagenome assemby with Metaspades, gene prediction with prodigal and annotation with diamond
  
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
    inputString <- paste("\"-1", fastqName1Trimmed, "-2", fastqName2Trimmed,"\"")
    pairedString="YES"
  } else {
    inputString <- paste("-s",fastqName1Trimmed)
    pairedString="NO"
  }
  metaspadesTemplScript <- METASPADES_TEMPLATE_SCRIPT
  metaspadesToBeExec <- paste0("metaspades.",sampleName,".sh")
  ##update template
  updateTemplateScriptCmd <- paste0("sed -e s/\"SAMPLE_NAME\"/", sampleName, "/g",
                                    " -e s/\"INPUT_FILE_STRING\"/", inputString, "/g ",
                                    " -e s/\"KMER_LIST\"/", param$metaspadesKmerList, "/g ",
                                    " -e s/\"ARE_READ_PAIRED\"/", pairedString, "/g ",
                                    metaspadesTemplScript, " >",
                                    metaspadesToBeExec)
  ezSystem(updateTemplateScriptCmd)
  ## run script
  metaspadesToBeExecCmd <- paste("bash",metaspadesToBeExec)
  ezSystem(metaspadesToBeExecCmd)
  
  ## rename output files
  #1) 
  oldContigFile <- "metaspadesResults/contigs.fa"
  newContigFile <- basename(output$getColumn("contigFile"))
  ezSystem(paste("mv",oldContigFile,newContigFile))
  #2) 
  oldProdigalFile <- "prodigalAnnotation.gff"
  newProdigalFile <- basename(output$getColumn("prodigalPredictionFile"))
  ezSystem(paste("mv",oldProdigalFile,newProdigalFile))
  #3) 
  oldIPSAnnFile <- "interProScanOut.gff3"
  newIPSAnnFile <- basename(output$getColumn("interproscanFile"))
  ezSystem(paste("mv",oldIPSAnnFile,newIPSAnnFile))
  #4) 
  oldKrakenLablesFile <- "interProScanOut.gff3"
  newKrakenLablesFile <- basename(output$getColumn("krakenLabelsFile"))
  ezSystem(paste("mv",oldKrakenLablesFile,newKrakenLablesFile))
}

##' @template app-template
##' @templateVar method ezMethodMetaspades()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMetaspades <-
  setRefClass("EzAppMetaspades",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMetaspades
                  name <<- "EzAppMetaspades"
                  appDefaults <<- rbind(metaspadesKmerList = ezFrame(Type="character",
                                                                  DefaultValue="69,79,89",
                                                                  Description="Comma-separated list of k-mer for the assembly."),
                                         ### TO DO: allow diamond option
                                        diamondEvalue = ezFrame(Type="numeric",
                                                                DefaultValue="0.05",
                                                                Description="Blast e-value cut-off."),
                                        diamondMaxSeqs = ezFrame(Type="integer",
                                                                 DefaultValue="30",
                                                                 Description="Blast maximum number of sequences to report.")
                                        
                                        
                  )
                }
              )
  )

