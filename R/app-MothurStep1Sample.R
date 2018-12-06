###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurStep1Sample = function(input=NA, output=NA, param=NA, 
                                        htmlFile="00index.html"){
  
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  library(scales)
  dataset = input$meta
  sampleName = input$getNames() 
  ### read fastq files and prepare inputs for Mothur
  ### are reads paired? should they be joined? 
  file1PathInDatset <- input$getFullPaths("Read1")
  unzipped1 <- paste0(sampleName,".R1.fastq")
  cpCmd <- paste0("gunzip -c ", file1PathInDatset, "  > ", unzipped1)
  ezSystem(cpCmd)
  if(param$paired){
    contigString = "make" 
    file2PathInDatset <- input$getFullPaths("Read2")
    unzipped2 <- paste0(sampleName,".R2.fastq")
    cpCmd2 <- paste0("gunzip -c ", file2PathInDatset, "  > ", unzipped2)
    ezSystem(cpCmd2)
    filesNameMothur <- cbind(input$getNames(),unzipped1,unzipped2)
    filesNameMothurFile <- paste(sampleName,"files",sep = ".")
    write.table(filesNameMothur,filesNameMothurFile,row.names = F,col.names = F, quote = F)
    initialFastaSuffix = "trim.contigs.fasta"
    initialGroupSuffix = "contigs.groups"
  }else{
    contigString = "###make"
    singleReadFileName <- paste0(sampleName,".R1")
    fastqFileToRead <- readDNAStringSet(paste0(singleReadFileName,".fastq"),format = "fastq")
    readNames <- ldply(strsplit(names(fastqFileToRead)," "), function(x)x[1])$V1
    groupFile <- data.frame(readNames, sampleName, stringsAsFactors = F)
    groupFileName <- paste0(sampleName,".R1.groups")
    write.table(groupFile,groupFileName, col.names = F, row.names = F, quote = F)
    fastaOutName <- paste0(singleReadFileName,".fasta")
    fastaFileToWrite <- writeXStringSet(fastqFileToRead, fastaOutName)
    initialFastaSuffix = "R1.fasta"
    initialGroupSuffix = "R1.groups"
  }
  
  ### is there at least a mock sample for the error estimate? The error estimates for the Non-mock samples will be ignored downstream
  if(param$mockSample){
    if (input$getColumn("Mock") == "Yes") {
      copyRefCmd <- paste("cp", param$referenceFasta,"./", sep = " ")
      ezSystem(copyRefCmd)
      mockString = "seq.error" 
      refString = param$referenceFasta
      oldErrFile <- paste(sampleName,
                          "unique.good.good.good.filter.unique.precluster.pick.pick.error.count",
                          sep = ".")
    }else{
      mockString = "###seq.error"  
      oldErrFile <- "mockErrorFileNotForDownstreamAnalysis.txt"
      write("There is no relevant info in this file",oldErrFile)
    }
  }else{
    mockString = "###seq.error"
  }
  ###

  ### cp silva reference locally
  cpSilvaRefCmd <- "cp /srv/GT/databases/silva/silva.bacteria.forMothur.fasta silva.bacteria.fasta"
  ezSystem(cpSilvaRefCmd)
  
  ### update batch file  with parameters and run mothur: step 1, identify region
  updateBatchCmd1 <- paste0("sed -e s/\"MIN_LEN\"/", param$minLen, "/g",
                           " -e s/\"MAX_LEN\"/", param$maxLen, "/g",
                           " -e s/\"INITIAL_SUFFIX_FASTA\"/", initialFastaSuffix, "/g",
                           " -e s/\"INITIAL_SUFFIX_GROUPS\"/", initialGroupSuffix, "/g",
                           " -e s/\"Mothur\"/", sampleName,"/g",
                           " -e s/\"###make\"/", contigString,"/g ",
                           file.path(METAGENOMICS_ROOT,FINAL_MOTHUR_WORKFLOW_TEMPLATE_STEP1), 
                           " >",
                           FINAL_MOTHUR_WORKFLOW_TEMPLATE_STEP1)
  ezSystem(updateBatchCmd1)
  cmdMothur1 = paste(MOTHUR_EXE,FINAL_MOTHUR_WORKFLOW_TEMPLATE_STEP1)
  ezSystem(cmdMothur1)
 
  
  ## rename output files
  ## Files needed for the report 
  #1) 
  oldRawSummary <- paste(sampleName,"summary",
                         sep = ".")
  newRawSummary <- basename(output$getColumn("RawDataSummary"))
  ezSystem(paste("mv",oldRawSummary,newRawSummary))
  
  #2) 
  oldDeduppedFileName <- paste(sampleName,"unique.summary",
                               sep = ".")
  newDeduppedFileName <- basename(output$getColumn("DeduppedSummary"))
  ezSystem(paste("mv",oldDeduppedFileName,newDeduppedFileName))
  
  #3) 
  oldLengthAndHomopFilteredFileName <- paste(sampleName,"unique.good.summary",
                                             sep = ".")
  newLengthAndHomopFilteredFileName <- basename(output$getColumn("LenAndHomopSummary"))
  ezSystem(paste("mv",oldLengthAndHomopFilteredFileName,newLengthAndHomopFilteredFileName))
  

  ## file needed for the second step 
  #5) 
  oldMappedFilteredFileName <- paste(sampleName,"unique.good.align",
                                     sep = ".")
  newMappedFilteredFileName <- basename(output$getColumn("alignedFile"))
  ezSystem(paste("mv",oldMappedFilteredFileName,newMappedFilteredFileName))
  
  oldMappedFilteredFileName <- paste(sampleName,"groups",
                                     sep = ".")
  newMappedFilteredFileName <- basename(output$getColumn("groupFile"))
  ezSystem(paste("mv",oldMappedFilteredFileName,newMappedFilteredFileName))

}

##' @template app-template
##' @templateVar method ezMethodMothurStep1Sample()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurStep1Sample <-
  setRefClass("EzAppMothurStep1Sample",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurStep1Sample
                  name <<- "EzAppMothurStep1Sample"
                  appDefaults <<- rbind(minLen = ezFrame(Type="integer",  DefaultValue="290",Description="Min length"),     
                                        maxLen= ezFrame(Type="integer",  DefaultValue="330",Description="Max length"),
                                        referenceFasta = ezFrame(Type="character",  DefaultValue="",Description="Mock reference seqeuences.")
                  )
                }
              )
  )
