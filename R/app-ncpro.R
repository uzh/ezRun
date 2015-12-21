###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodNcpro = function(input=NA, output=NA, param=NA){
  setwdNew(basename(output$getColumn("Report")))
  param$readCountsBarplot = basename(output$getColumn("TrimCounts"))
  ncpro(input=input, dataset=input$meta, param=param)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodNcpro
##' @templateVar htmlArg )
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
##' Trimming is done using function trimMirna(). 
##' @param input     list of input configuration
##' @param dataset   dataframe with description of input dataset
##' @param param     further configuration parameters
ncpro = function(input, dataset, param=NULL){
  samples = rownames(dataset)
  fqFiles = input$getFullPaths(param, "Read1")
  names(fqFiles) = samples
  adapter = unique(dataset$Adapter1)
  stopifnot(length(adapter) == 1)
  jobList = lapply(fqFiles, function(fq){list(input=fq, output=file.path(getwd(), sub(".gz$", "", basename(fq))))})
  .myFunc = function(job, param){
    trimMirna(input=job$input, output=job$output, adapter=adapter, param=param)
  }
  buildName = param$ezRef["refBuildName"]
  trimmedFastqFiles = unlist(ezMclapply(jobList,.myFunc,param=param,mc.cores=as.numeric(param[['cores']]),mc.preschedule =FALSE, mc.set.seed=FALSE))
  ncproConfigFile = list.files(paste0(NCPRO_ANNOTATION_DIR,"/config-templates/"), pattern=paste0('-', buildName,'-'), full.names=T)[1]
  if(is.na(ncproConfigFile))
    stop(paste0("No ncpro config-template for Genome-Build ", param[['refBuild']]," available."))
  jobDir = getwd()
  workDir = file.path(jobDir, "ncpro")
  ezSystem(paste(file.path(NCPRO_DIR, "bin/ncPRO-deploy"), "-o", workDir))
  x = readLines(ncproConfigFile)
  x = sub("^N_CPU.*", paste("N_CPU =", param[['cores']]), x)
  writeLines(x, file.path(workDir, "param-ncrna.txt"))
  rawDir = file.path(workDir, "rawdata")
  ## link the files to the raw directory
  for (sm in samples){
    fqFile = paste0(rawDir, "/", sm, ".fastq")
    ezSystem(paste("ln -s", trimmedFastqFiles[sm], fqFile))
  }
  refIndex = getBowtieReference(param)
  Sys.setenv(BOWTIE_INDEXES=dirname(refIndex))
  setwd(workDir)
  
  workflowSteps = c("processRead","mapGenome","mapGenomeStat","processBam","mapAnnOverview",
                    "overviewRfam","generateNcgff","ncrnaProcess", "ncrnaTracks","sigRegion","html_builder")
  ezSystem(paste(file.path(NCPRO_DIR, "bin/ncPRO-seq"), "-c", "param-ncrna.txt","-s processRead ",">ncpro.log"))
  for(i in 2:length(workflowSteps)){
    ezSystem(paste(file.path(NCPRO_DIR, "bin/ncPRO-seq"), "-c", "param-ncrna.txt","-s",workflowSteps[i],">>ncpro.log"))
  }
  stopifnot(file.exists("report.html"))
  stopifnot(!grepl("^make.*Error", readLines("ncpro.log")))
  setwd(jobDir)
  readCounts = data.frame(row.names=samples)
  if (!is.null(dataset$"Read Count") && is.numeric(dataset$"Read Count") && all(dataset$"Read Count" > 0)){
    readCounts$untrimmed = dataset[samples, "Read Count"]    
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
