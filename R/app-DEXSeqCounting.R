###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


#' Run-method for ezApp EzAppDEXSeqCounting
#' 
ezMethodDEXSeqCounting <- function(input=NA, output=NA, param=NA){
  
  ### # check whether GFF formatted annotation is available
  sGtfFile <- param$ezRef@refFeatureFile
  sGffFile <- gsub("gtf$", "gff", sGtfFile)
  if (!file.exists(sGffFile))
    convertGtfToGff(psGtfFile = sGtfFile)
  ### # if there is no link in current working directory to the gff File,
  ### #  created it
  sGffBaseFn <- basename(sGffFile)
  if (!file.exists(sGffBaseFn)){
    sLnCmd <- paste("ln -s", sGffFile, sGffBaseFn)
  }
  ### # do the counting, get the bam files from input
  bamFile = input$getFullPaths(param, "BAM")

  ### # call counting routine
  vCountFiles <- sapply(bamFile, runCountSingleBam, sGffFile)
  
  ### # write count files back to input$file
  writeCountFilesToMeta(pvCountFiles = vCountFiles, input = input)
  
  ### # clean up, remove link to gff file
  if (file.exists(sGffBaseFn))
    unlink(sGffBaseFn)
  
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodDEXSeqCounting
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppDEXSeqCounting <-
  setRefClass("EzAppDEXSeqCounting",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDEXSeqCounting
                  name <<- "EzAppDEXSeqCounting"
                }
              ))


#' Convert annotation file from GTF to GFF
#' 
#' @param psGtfFile   name of the GTF annotation file
convertGtfToGff <- function(psGtfFile) {
  ### # assume that we are somewhere in a scratch working directory
  ### #  because conversion scripts want input in current directory,
  ### #  we create a link to the gtf file, but not, if the working 
  ### #  directory is the same as the directory of the annotation file
  cat(" * Converting GTF to GFF ...\n")
  sCurWd <- getwd()
  sAnnotDir <- dirname(psGtfFile)
  setwd(sAnnotDir)
  sGtfFn <- basename(psGtfFile)
  ### # convert using the python scripts
  sGffFn <- gsub("gtf$", "gff", sGtfFn)
  ### # check whether the path exists
  if (!exists("DEXSEQ_PREPARE")) {
    DEXSEQ_PREPARE <- lGetPyScriptPaths()$DEXSEQ_PREPARE
  }
  sPyConvCmd <- paste(DEXSEQ_PREPARE, sGtfFn, sGffFn)
  ezSystem(sPyConvCmd)
  cat("  ... created", sGffFn, "\n")
  setwd(sCurWd)
  invisible(TRUE)
}

#' Run counts for a single BAM file
#' 
runCountSingleBam <- function(psBamFile, psGffFile){
  sSamCmd <- paste(SAMTOOLS, "view -h", psBamFile)
  ### # run counting on sam file
  sCountBaseFn <- gsub("bam$", "count", basename(psBamFile))
  if (!exists("DEXSEQ_COUNT")){
    DEXSEQ_COUNT <- lGetPyScriptPaths()$DEXSEQ_COUNT
  }
  sPyCountCmd <- paste(sSamCmd, "|", DEXSEQ_COUNT, psGffFile, "-", sCountBaseFn)
  ezSystem(sPyCountCmd)
  sCountDir <- getwd()
  return(file.path(sCountDir, sCountBaseFn))
}


#' Write names of countfiles back into the input file
#' 
writeCountFilesToMeta <- function(pvCountFiles, input) {
  ### # add column with counts to the meta information
  input$meta$Count <- pvCountFiles
  ### # write the extended meta information back to the file
  write.table(input$meta, file = input$file, quote = FALSE, sep = "\t")
}


#' Get list with required python script paths
#' 
lGetPyScriptPaths <- function(){
  if (!exists("PYTHON_CMD")){
    PYTHON_CMD='PYTHONPATH="/usr/local/ngseq/lib/python/:/usr/local/ngseq/lib/python2.7:/usr/local/ngseq/lib/python2.7/dist-packages" /usr/local/ngseq/bin/python'
  }
  return(list(DEXSEQ_PREPARE = paste(PYTHON_CMD, file.path(system.file(package = "DEXSeq", "python_scripts"), "dexseq_prepare_annotation.py")),
              DEXSEQ_COUNT = paste(PYTHON_CMD, file.path(system.file(package = "DEXSeq", "python_scripts"), "dexseq_count.py"))))
  
}

