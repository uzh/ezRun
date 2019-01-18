###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNcpro = function(input=NA, output=NA, param=NA){
  setwdNew(basename(output$getColumn("Report")))
  param$readCountsBarplot = basename(output$getColumn("TrimCounts"))
  rownames(input$meta) = gsub('_','-',rownames(input$meta))  
  ncpro(input=input, dataset=input$meta, param=param)
  postProcessResults(dataset=input$meta, psInputFn=input$file, psReportDir=output$getColumn("Report"))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodNcpro(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppNcpro <-
  setRefClass("EzAppNcpro",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodNcpro
                  name <<- "EzAppNcpro"
                }
              )
  )

## TODO: make sure there's no conflict with input and refactor trimMirna
##' Analysis of small RNA sequences using ncpro
##' 
##' Trimming is done using ezMethodTrim which is similar to trimMirna(). 
##' @param input     list of input configuration
##' @param dataset   dataframe with description of input dataset
##' @param param     further configuration parameters
ncpro = function(input, dataset, param=NULL){
  samples = rownames(dataset)
  fqFiles = input$getFullPaths("Read1")
  names(fqFiles) = samples
  adapter = unique(dataset$Adapter1)
  stopifnot(length(adapter) == 1)
  ### # trimming using ezMethodTrim
  # set parameter defaults for trimming, to work similarly to trimMirna
  param[['trimAdapter']]            <-  TRUE
  param[['minTailQuality']]         <-  20
  param[['minAvgQuality']]          <-  4
  param[['minReadLength']]          <-  18
  refObjTrimResult <- ezMethodTrim(input, output=NA, param)
  buildName = param$ezRef["refBuildName"]
  ### # get names of trimmed read files from output of ezMethodTrim
  trimmedFastqFiles <- as.vector(refObjTrimResult$meta[,"Read1 [File]"])
  names(trimmedFastqFiles) <- samples
  ncproConfigFile = list.files(paste0(NCPRO_ANNOTATION_DIR,"/config-templates/"), pattern=paste0('-', buildName,'-'), full.names=T)[1]
  if(is.na(ncproConfigFile))
    stop(paste0("No ncpro config-template for Genome-Build ", param[['refBuild']]," available."))
  jobDir = getwd()
  workDir = file.path(jobDir, "ncpro")
  ezSystem(paste("ncPRO-deploy", "-o", workDir))
  x = readLines(ncproConfigFile)
  x = sub("^N_CPU.*", paste("N_CPU =", param[['cores']]), x)
  writeLines(x, file.path(workDir, "param-ncrna.txt"))
  rawDir = file.path(workDir, "rawdata")
  ## link the files to the raw directory
  for (sm in samples){
    fqFile = paste0(rawDir, "/", sm, ".fastq")
    ezSystem(paste("ln -s", trimmedFastqFiles[sm], fqFile))
  }
  # we now specify the absolute path of the index in the config file
  # refIndex = getBowtieReference(param)
  # Sys.setenv(BOWTIE_INDEXES=dirname(refIndex))
  setwd(workDir)
  
  workflowSteps = c("processRead","mapGenome","mapGenomeStat","processBam","mapAnnOverview",
                    "overviewRfam","generateNcgff","ncrnaProcess", "ncrnaTracks","sigRegion","html_builder")
  ezSystem(paste("ncPRO-seq", "-c", "param-ncrna.txt","-s processRead ",">ncpro.log"))
  for(i in 2:length(workflowSteps)){
    ezSystem(paste("ncPRO-seq", "-c", "param-ncrna.txt","-s",workflowSteps[i],">>ncpro.log"))
  }
  stopifnot(file.exists("report.html"))
  stopifnot(!grepl("^make.*Error", readLines("ncpro.log")))
  setwd(jobDir)
  readCounts = data.frame(row.names=samples)
  if (!is.null(dataset$"Read Count") && is.numeric(dataset$"Read Count") && all(dataset$"Read Count" > 0)){
    readCounts$untrimmed = dataset[samples, "Read Count"]    
  } else {
    readCounts$untrimmed = countReadsInFastq(fqFiles[samples])
  }
  readCounts$remaining = countReadsInFastq(trimmedFastqFiles)
  ezWrite.table(readCounts, "trimCounts.txt")
  readCounts$removed = readCounts$untrimmed - readCounts$remaining
  
#   plotCmd = expression({
#     par(mar=c(12, 4.1, 4.1, 2.1))  
#     barplot(t(as.matrix(readCounts[ , c("remaining", "removed")])), las=2, border=NA,
#             main="Read Counts after trimming", legend.text=TRUE, col=c("gray30", "gray"))
#   })
#   unusedLink = ezImageFileLink(plotCmd, file=param$readCountsBarplot, width=400 + nrow(readCounts) * 10, height=700)
  
  png(param$readCountsBarplot, width=400 + nrow(readCounts) * 10, height=700)
  par(mar=c(12, 4.1, 4.1, 2.1))  
  barplot(t(as.matrix(readCounts[ , c("remaining", "removed")])), las=2, border=NA,
          main="Read Counts after trimming", legend.text=TRUE, col=c("gray30", "gray"))
  dev.off()
  
  ezSystem(paste("pigz", paste(trimmedFastqFiles, collapse=" ")))
  ezSystem("rm -f data/*bed")
  ezSystem("rm -f data/*tmp")
  ezSystem("rm -f data/*fas")
  ezSystem("rm -f data/*gff")
  ezSystem("rm -f manuals")
  ezSystem("rm -f annotation")
  ezSystem("rm -rf rawdata")
  setwd(jobDir)
}


#' Postprocessing counts produced as results from ncpro 
#' 
#' \code{postProcessResults} takes counts for categories all, mature and precursor
#' and creates a separate result file for each sample. This splitting of count 
#' results is done in function splitCounts. In the tsv-formatted input metadata the  
#' path to the read files is replaced by the path to the count result files. This 
#' replacement is done in function modifyInput.
#' @param input   input parameters
#' @param psReportDir   directory where ncpro results are found 
postProcessResults <- function(dataset, psInputFn, psReportDir) {
  # setting directories to where count files are
  jobDir <- getwd()
  # extract samples
  vSamples <- rownames(dataset)
  # get some parameter settings specific for splitting counts
  lLocalParam <- lGetGlobalCountParam()
  # parameters to different count categories
  lCountParam <- lLocalParam$countParam
  # extract smRNA categories
  vCountCategories <-names(lCountParam)
  # split the counts using apply over the count categories
  lapply(vCountCategories, FUN=splitCounts, pvSamples=vSamples, plCountParam=lCountParam )
  # modify input file such that Reads column is replaced by counts column
  lapply(vCountCategories, FUN=modifyInput, pdataset=dataset, psInputFn=psInputFn, psReportDir=psReportDir)
  # reset directory back to jobDir
  setwd(jobDir)
}

#' split count result into different files and put them into result dir
#' 
#' \code{splitCounts} for a given count category, count results 
#' for all samples given in pvSamples are written to different 
#' files, according to parameters given in plCountParam
#' 
#' @param psCountType    for which count category {all, mature or precursor}
#'                       count results should be split
#' @param pvSamples      vector of sample names
#' @param plCountParam   list of parameters specifying existing result files
splitCounts <- function(psCountType, pvSamples, plCountParam) {
  # extract file, directory and colIdx parameter from the list plCountParam
  sInputFile <- plCountParam[[psCountType]]$inputFile
  sResultDir <- plCountParam[[psCountType]]$resultDir
  nSampleNameColIdx <- plCountParam[[psCountType]]$nSampleNameColIdx
  # check whether input file exists
  stopifnot(file.exists(sInputFile))
  # create result directory
  if (!dir.exists(sResultDir)) 
    dir.create(path = sResultDir)
  # read input file
  dfCountInput <- read.delim(file = sInputFile, stringsAsFactors = FALSE,check.names=F)
  # names of original count df
  vColNames <- names(dfCountInput)
  sSampleNameColName <- vColNames[nSampleNameColIdx]
  # for every sample, extract counts from original results
  for (sam in pvSamples) {
    # check that sample name occurs
    if (!(sam %in% vColNames)) {
      stop("ERROR in splitCounts, sample: ", sam, " not found in result file", sInputFile, "\n")  
    }
    # generate name of result file
    sResFn <- file.path(sResultDir, paste0(sam,".txt"))
    # combine counts of current sample into separate dataframe
    dfSamSplit <- cbind(dfCountInput[,nSampleNameColIdx], dfCountInput[,sam])
    # set the column names
    colnames(dfSamSplit) <- c('Identifier','matchCounts') #c(sSampleNameColName, sam)
    # write output to result file
    write.table(dfSamSplit, file = sResFn, quote = FALSE, sep = "\t", row.names = FALSE)
  }
}

#' input is modified to show the count results
#' 
#' \code{modifyInput} for a given count category, tsv-formatted 
#' metadata is modified such that paths to read files are replaced 
#' with paths to output count files
modifyInput <- function(psCountCategory, pdataset, psInputFn, psReportDir) {
  # localInput metadata dataframe
  dfDataSet <- pdataset
  # extract sample names
  vSampleNames <- rownames(dfDataSet)
  # generate resultpaths
  vResultPaths <- as.vector(sapply(vSampleNames, 
                                   function(sam){
                                     file.path(psReportDir, psCountCategory, 
                                               paste(sam,"txt", sep="."))}))
  
  #featureLevel (Wert='smRNA'), refFeatureFile (ohne Inhalt) und refBuild
  featureLevel <- rep("smRNA", nrow(dfDataSet))
  refFeatureFile <- rep("", nrow(dfDataSet))
  refBuild <- rep("", nrow(dfDataSet))
  # replace column "Read1 [File]" in dfDataSet with column "Count [File]"
  dfDataSet[,"Read1 [File]"] <- vResultPaths
  colnames(dfDataSet)[which(colnames(dfDataSet) == "Read1 [File]")] <- "Count [File]"
  vColNames <- colnames(dfDataSet)
  # cbind the new dataframe together
  dfDataSet <- cbind(vSampleNames,dfDataSet,featureLevel,refFeatureFile,refBuild)
  colnames(dfDataSet) <- c("Name", vColNames, "featureLevel", "refFeatureFile", "refBuild")
  # write new data frame to file
  write.table(dfDataSet, file = file.path('.',paste0(psCountCategory,'_dataset.tsv')), quote = FALSE, sep = "\t", row.names = FALSE)
}

#' specification of default values used for count splitting
lGetGlobalCountParam <- function(){
  return(list(countParam = list(allRNA         = list(inputFile = "ncpro/doc/all_samples_all_subfamcov.data", 
                                                      resultDir = "allRNA",
                                                      nSampleNameColIdx = 1),
                               mature_miRNA    = list(inputFile = "ncpro/doc/mature_miRNA_miRNA_e_+2_+2_all_samples_subfamcov.data", 
                                                      resultDir = "mature_miRNA",
                                                      nSampleNameColIdx = 1),
                               precursor_miRNA = list(inputFile = "ncpro/doc/precursor_miRNA_miRNA_all_samples_subfamcov.data",
                                                      resultDir = "precursor_miRNA",
                                                      nSampleNameColIdx = 1))))
}

