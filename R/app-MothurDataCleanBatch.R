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
    contigString = "###make.contigs" 
    file2PathInDatset <- input$getFullPaths("Read2")
    cpCmd2 <- paste0("gunzip -c ", file2PathInDatset, "  > ", sampleName,".R2",".fastq")
    ezSystem(cpCmd2)
  }else{
    contigString = "###make" 
  }
    projNum <- dirname(param$resultDir)
  
  ### cp silva reference locally
  cpSilvaRefCmd <- "cp /srv/GT/databases/silva/silva.bacteria.forMothur.fasta silva.bacteria.fasta"
  ezSystem(cpSilvaRefCmd)
  
  

  ### update batch file Illumina with parameters and run mothur: step 1, identify region
  updateBatchCmd <- paste0("sed -e s/\"MIN_LEN\"/", param$minLen, "/g",
                                   " -e s/\"MAX_LEN\"/", param$maxLen, "/g",
                                   " -e s/\"Mothur\"/", sampleName,"/g",
                                   " -e s/\"###make\"/", contigString,"/g",
                           MOTHUR_DATA_CLEAN_BATCH_TEMPLATE_STEP1, " >", MOTHUR_DATA_CLEAN_BATCH_STEP1)
  ezSystem(updateBatchCmd)
  cmdMothur = paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_STEP1)
  ezSystem(cmdMothur)
  ### extract region
  summaryFileToExtractRegion <- paste(sampleName,"unique.good.good.summary", sep = ".")
  summaryOfAlign <- read.delim(summaryFileToExtractRegion, stringsAsFactors = FALSE, header = T)
  regionStartCoord <- quantile(summaryOfAlign$start, probs = seq(0, 1, 0.025))["5%"]
  regionEndCoord <- quantile(summaryOfAlign$end, probs = seq(0, 1, 0.025))["95%"]
  
  ### update batch file Illumina with parameters and run mothur: step 2, precluster and remove non-bacterial reads
  updateBatchCmd_step2 <- paste0("sed -e s/\"CUTOFF\"/", param$cutOff, "/g",
                                   " -e s/\"DIFFS\"/", param$diffs, "/g",
                                   " -e s/\"START_COORD\"/", regionStartCoord, "/g",
                                   " -e s/\"END_COORD\"/", regionEndCoord, "/g",
                                   " -e s/\"Mothur\"/" ,sampleName, "/g ",
                                   MOTHUR_DATA_CLEAN_BATCH_TEMPLATE_STEP2,
                                   " >", MOTHUR_DATA_CLEAN_BATCH_STEP2)
  ezSystem(updateBatchCmd_step2)
  cmdMothur_step2= paste(MOTHUR_EXE,MOTHUR_DATA_CLEAN_BATCH_STEP2)
  ezSystem(cmdMothur_step2)
  ## rename files 
  oldCountFileName <- paste(sampleName,
                            "unique.good.good.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table",
                            sep = ".")
  newCountFileName <- output$CountTable
  renameCountFileCmd <- ezSystem(paste("mv",oldCountFileName,newCountFileName))
    
  oldClusterFileName <- paste(sampleName,
                              "unique.good.good.good.filter.unique.precluster.pick.pick.fasta",
                              sep = ".")
  newClusterFileName <- output$PreClusteredFastaFile
  renameCLusterFileCmd <-  ezSystem(paste("mv",oldClusterFileName,newClusterFileName))
    
  oldTaxonomyFileName <- paste(sampleName,
                               "unique.good.good.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy",
                               sep = ".")
  newTaxonomyFileName <- output$TaxonomyFile
  renameTaxonomyFileCmd <-  ezSystem(paste("mv",oldTaxonomyFileName,newTaxonomyFileName))
  
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
                  appDefaults <<- rbind(cutOff = ezFrame(Type="integer",  DefaultValue="80",Description="Cut-off for taxonomy assignment"),
                                        diffs= ezFrame(Type="integer",  DefaultValue="2",Description="Differences allowed in the pre.cluster step.Should be 1 every 100 bases"),
                                        minLen = ezFrame(Type="integer",  DefaultValue="290",Description="Min length"),     
                                        maxLen= ezFrame(Type="integer",  DefaultValue="330",Description="Max length")
                  )
                }
              )
  )
