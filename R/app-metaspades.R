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
    inputStringAss <- paste("\"-1", fastqName1, "-2", fastqName2,"\"")
    inputStringBwt2 <- paste("\"-1", fastqName1, "-2", fastqName2,"\"")
  } else {
    inputStringAss <- paste("-s",fastqName1)
    inputStringBwt2 <- paste("-U",fastqName1)
  }
  ##update template
  updateTemplateScriptCmd <- paste0("sed -e s/\"SAMPLE_NAME\"/", sampleName, "/g",
                                    " -e s/\"INPUT_FILE_STRING_BWT\"/", inputStringBwt2, "/g ",
                                    " -e s/\"INPUT_FILE_STRING_ASS\"/", inputStringAss, "/g ",
                                    " -e s/\"KMER_LIST\"/", param$metaspadesKmerList, "/g ",
                                    " -e s/\"NUM_CPU\"/", param$cores, "/g ",
                                    " -e s/\"ANNOTATION\"/", param$annotation, "/g ",
                                    " -e s/\"MAXIMUM_SEQS\"/", param$diamondMaxSeqs, "/g ",
                                    " -e s/\"E_VALUE_CUTOFF\"/", param$diamondEvalue, "/g ",
                                    file.path(METAGENOMICS_ROOT,METASPADES_TEMPLATE_SCRIPT), " >",
                                    METASPADES_TEMPLATE_SCRIPT)
  ezSystem(updateTemplateScriptCmd)
  ## run script
  metaspadesToBeExecCmd <- paste("bash",METASPADES_TEMPLATE_SCRIPT)
  ezSystem(metaspadesToBeExecCmd)

  ## place output files
  #1) contigs
  oldContigFile <- "metaspadesResults/contigs.fasta"
  newContigFile <- basename(output$getColumn("contigFile"))
  ezSystem(paste("cp",oldContigFile,newContigFile))
  #2) 
  oldProdigalFile <- "prodigalAnnotation.gff"
  newProdigalFile <- basename(output$getColumn("prodigalPredictionFile"))
  ezSystem(paste("cp",oldProdigalFile,newProdigalFile))
  #3) 
  oldIPSAnnFile <- "interProScanOut.gff3"
  newIPSAnnFile <- basename(output$getColumn("interproscanFile"))
  ezSystem(paste("cp",oldIPSAnnFile,newIPSAnnFile))
  #4)   
  krakenFile <- "kraken.labels"
  binFiles <- list.files("maxbinOut", pattern = "fasta", full.names = T)
  binSummaryFile <- summaryMetagenomeBins(krakenFile,binFiles)
  binSummaryFileName <- basename(output$getColumn("binSummaryFile"))
  write.table(binSummaryFile,binSummaryFileName,
              row.names = F, col.names = T, quote = F,sep = "\t")
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

