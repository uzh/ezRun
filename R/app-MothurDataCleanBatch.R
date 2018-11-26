###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurDataCleanBatch = function(input=NA, output=NA, param=NA, 
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
  cpCmd <- paste0("gunzip -c ", file1PathInDatset, "  > ", sampleName,".R1",".fastq")
  ezSystem(cpCmd)
  if(param$paired){
    contigString = "make" 
    file2PathInDatset <- input$getFullPaths("Read2")
    cpCmd2 <- paste0("gunzip -c ", file2PathInDatset, "  > ", sampleName,".R2",".fastq")
    ezSystem(cpCmd2)
    initialFastaSuffix = "trim.contigs.fasta"
    initialGroupSuffix = "contigs.groups"
  }else{
    contigString = "###make"
    singleReadFileName <- paste0(sampleName,".R1")
    fastqFileToRead <- readDNAStringSet(paste0(singleReadFileName,".fastq"),format = "fastq")
    readNames <- ldply(strsplit(names(fastqFileToRead)," "), function(x)x[1])$V1
    groupFile <- data.frame(readNames, singleReadFileName, stringsAsFactors = F)
    groupFileName <- paste0(sampleName,".R1.groups")
    write.table(groupFile,groupFileName, col.names = F, row.names = F, quote = F)
    fastaOutName <- paste0(singleReadFileName,".fasta")
    fastaFileToWrite <- writeXStringSet(fastqFileToRead, fastaOutName)
    initialFastaSuffix = "R1.fasta"
    initialGroupSuffix = "R1.groups"
  }
  
  ### is there at least a mock sample for the error estimate? The error estimates for the Non-mock samples will be ignored downstream
  if(param$mockSample){
    mockString = "seq.error" 
  }else{
    mockString = "###seq.error" 
  }
  ###
  projNum <- dirname(param$resultDir)
  
  ### cp silva reference locally
  cpSilvaRefCmd <- "cp /srv/GT/databases/silva/silva.bacteria.forMothur.fasta silva.bacteria.fasta"
  ezSystem(cpSilvaRefCmd)
  
  ### update batch file  with parameters and run mothur: step 1, identify region
  updateBatchCmd <- paste0("sed -e s/\"MIN_LEN\"/", param$minLen, "/g",
                                   " -e s/\"MAX_LEN\"/", param$maxLen, "/g",
                                   " -e s/\"INITIAL_SUFFIX_FASTA\"/", initialFastaSuffix, "/g",
                                   " -e s/\"INITIAL_SUFFIX_GROUPS\"/", initialGroupSuffix, "/g",
                                   " -e s/\"Mothur\"/", sampleName,"/g",
                                   " -e s/\"###make\"/", contigString,"/g ",
                           MOTHUR_DATA_CLEAN_BATCH_TEMPLATE_STEP1, " >",
                           MOTHUR_DATA_CLEAN_BATCH_STEP1)
  ezSystem(updateBatchCmd)
  cmdMothur = paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_STEP1)
  ezSystem(cmdMothur)
  ### extract region
  summaryFileToExtractRegion <- paste(sampleName,"unique.good.good.summary", sep = ".")
  summaryOfAlign <- read.delim(summaryFileToExtractRegion, stringsAsFactors = FALSE, header = T)
  regionStartCoord <- quantile(summaryOfAlign$start, probs = seq(0, 1, 0.025))["95%"]
  regionEndCoord <- quantile(summaryOfAlign$end, probs = seq(0, 1, 0.025))["5%"]

  ### update batch file  with parameters and run mothur: step 2 and, precluster, remove non-bacterial reads and generate final cluster
  updateBatchCmd_step2_3 <- paste0("sed -e s/\"CUTOFF_TAXON\"/", param$cutOffTaxonomy, "/g",
                                   " -e s/\"CUTOFF_CLUST\"/", param$cutOffCluster, "/g",
                                   " -e s/\"DIFFS\"/", param$diffs, "/g",
                                   " -e s/\"START_COORD\"/", regionStartCoord, "/g",
                                   " -e s/\"END_COORD\"/", regionEndCoord, "/g",
                                   " -e s/\"Mothur\"/" ,sampleName, "/g",
                                   " -e s/\"###seq.error\"/", mockString,"/g ",
                                   MOTHUR_DATA_CLEAN_BATCH_TEMPLATE_STEP2_3,
                                   " >", MOTHUR_DATA_CLEAN_BATCH_STEP2_3)
  ezSystem(updateBatchCmd_step2_3)
  cmdMothur_step2_3= paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_STEP2_3)
  ezSystem(cmdMothur_step2_3)
  
  ## rename output files
  #1) 
  oldReadsToCountFileName <- paste(sampleName,
                                   "unique.good.good.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table",
                                   sep = ".")
  newReadsToCountFileName <- basename(output$getColumn("ReadsCountTable"))
  ezSystem(paste("mv",oldReadsToCountFileName,newReadsToCountFileName))
  #2) 
  oldPreclusterFileName <- paste(sampleName,
                                 "unique.good.good.good.filter.unique.precluster.pick.pick.fasta",
                                 sep = ".")
  newPreclusterFileName <- basename(output$getColumn("PreClusteredFastaFile"))
  ezSystem(paste("mv",oldPreclusterFileName,newPreclusterFileName))
  #3) 
  oldReadsToTaxonomyFileName <- paste(sampleName,
                                      "unique.good.good.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy",
                                      sep = ".")
  newReadsToTaxonomyFileName <- basename(output$getColumn("ReadsTaxonomyFile"))
  ezSystem(paste("mv",oldReadsToTaxonomyFileName,newReadsToTaxonomyFileName))
  #4) 
  oldOTUsToCountFileName <- paste(sampleName,
                                  "unique.good.good.good.filter.unique.precluster.pick.pick.opti_mcc.shared",
                                  sep = ".")
  newOTUsToCountFileName <- basename(output$getColumn("OTUsCountTable"))
  ezSystem(paste("mv",oldOTUsToCountFileName,newOTUsToCountFileName))

  #5) 
  oldErrFile <- paste(sampleName,
                      "unique.good.good.good.filter.unique.precluster.pick.pick.error.count",
                      sep = ".")
  newErrFile <- basename(output$getColumn("ErrorFile"))
  ezSystem(paste("mv",oldErrFile,newErrFile))
  #6) 
  oldStepConvFile <- paste(sampleName,
                           "unique.good.good.good.filter.unique.precluster.pick.pick.opti_mcc.steps",
                           sep = ".")
  newStepConvFile <- basename(output$getColumn("stepConvergence"))
  ezSystem(paste("mv",oldStepConvFile,newStepConvFile))
  #7) 
  oldOTUsToTaxFileName <- paste(sampleName,
                                "unique.good.good.good.filter.unique.precluster.pick.pick.opti_mcc",
                                params$cutOffCluster, "cons.taxonomy",
                                sep = ".")
  newOTUsToTaxFileName <- basename(output$getColumn("OTUsToTaxonomyFile"))
  ezSystem(paste("mv",oldOTUsToTaxFileName,newOTUsToTaxFileName))
}

##' @template app-template
##' @templateVar method ezMethodMothurDataCleanBatch()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurDataCleanBatch <-
  setRefClass("EzAppMothurDataCleanBatch",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurDataCleanBatch
                  name <<- "EzAppMothurDataCleanBatch"
                  appDefaults <<- rbind(cutOffTaxonomy = ezFrame(Type="integer",  DefaultValue="80",Description="Cut-off for taxonomy assignment"),
                                        diffs= ezFrame(Type="integer",  DefaultValue="2",Description="Differences allowed in the pre.cluster step.Should be 1 every 100 bases"),
                                        minLen = ezFrame(Type="integer",  DefaultValue="290",Description="Min length"),     
                                        maxLen= ezFrame(Type="integer",  DefaultValue="330",Description="Max length"),
                                        cutOffCluster = ezFrame(Type="numeric",  DefaultValue="0.03",Description="Cut-off for OTU clustering."),
                                        referenceFasta = ezFrame(Type="character",  DefaultValue="",Description="Mock reference seqeuences.")
                  )
                }
              )
  )
