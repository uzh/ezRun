###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodSGSeq = function(input=NA, output=NA, param=NA){
  library(SGSeq)
  library(plyr)
  library(GenomicFeatures)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(ezRun)
  library(GenomicFeatures)
## create and save TxDb object of not there
  dataset = input$meta
  currentRoot <- unlist(strsplit(GENOMES_ROOT,":"))[1]
  txdbGtfOrigin <- file.path(currentRoot,param$refBuild,"Genes",param$refFeatureFile)
  txdbDatabase <- file.path(currentRoot,param$refBuild,"Genes","txdb.sqlite")
  if (!file.exists(txdbDatabase)) {
    txDbObject <- makeSGSeqtxdbObject(param$refBuild,param$refFeatureFile,txdbDatabase)
  } else {
    txDbObject <- loadDb(txdbDatabase)
  }
  ##

## Load bam file   
  bamFile = input$getFullPaths("BAM")
  localBamFile = getBamLocally(bamFile)
  if(localBamFile != bamFile){
    on.exit(file.remove(c(localBamFile, paste0(localBamFile, ".bai"))),
            add=TRUE)
  }
  
## create R object  
  sampleName = input$getNames()
  inputBamDataFrame <- data.frame(sample_name=sampleName,file_bam=bamFile, stringsAsFactors = FALSE)
  sampleInfoComplete <- getBamInfo(inputBamDataFrame)
  bamFile <- sampleInfoComplete$file_bam
  usedChr <- scanBamHeader(bamFile)[[1]]$text[names(scanBamHeader(sampleInfoComplete$file_bam[1])[[1]]$text) =="@SQ"]
  usedChr <- gsub("SN:", "", ldply(usedChr)$V1)
  txdb <- keepSeqlevels(txDbObject,usedChr)
  txf_ucsc <- convertToTxFeatures(txdb)
  restoreSeqlevels(txdb)
  sgfc_ucscF <- analyzeFeatures(sampleInfoComplete, features = txf_ucsc)
  sgfc_ucscV <- analyzeVariants(sgfc_ucscF)
  outRdataFile <- basename(output$getColumn("SgSeqRDataFile"))
  save(sgfc_ucscF,sgfc_ucscV, file = outRdataFile)
  
## create freq and count file 
  ### get var freqs
  varFreq <- data.frame(freq = variantFreq(sgfc_ucscV), stringsAsFactors = FALSE)
  outVarFreq <- basename(output$getColumn("SgSeqVarFreqFile"))
  write.table(varFreq, outVarFreq, col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  ### get counts
  varCounts <- getSGVariantCounts(rowRanges(sgfc_ucscV), sample_info = sampleInfoComplete)
  countValues <- counts(varCounts)
  vid <- variantID(varCounts)
  eid <- eventID(varCounts)
  inputForDexSeq <- data.frame(paste(eid,vid,sep = ":"), countValues, stringsAsFactors = FALSE)
  outVarFreq <- basename(output$getColumn("SgSeqCountFile"))
  write.table(inputForDexSeq, outVarFreq, row.names = FALSE,
              col.names = FALSE, quote = FALSE, sep = "\t")
}



##' @template app-template
##' @templateVar method ezMethodSGSeq(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppSGSeq <-
  setRefClass("EzAppSGSeq",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSGSeq
                  name <<- "EzAppSGSeq"
                  appDefaults <<- rbind(gtfFeatureType = ezFrame(Type="character",  DefaultValue="exon",  
                                                                 Description="which gtf feature types to use; with Ensembl GTF files; use 'transcript' to count also intronic reads"),
                                        upstreamFlanking=ezFrame(Type="integer", DefaultValue="250", 
                                                                 Description="the number of bases upstream of TSS"))
                }
              )
  )







