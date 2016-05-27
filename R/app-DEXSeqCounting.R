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
  ### # gff will be placed in actual working directory, hence no 
  ### #  soft links will be required
  sGffFile <- gsub("gtf$", "gff", basename(sGtfFile))
  if(ezIsSpecified(param$gff_file))
    sGffFile <- param$gff_file
  if (!file.exists(sGffFile))
    convertGtfToGff(psGtfFile = sGtfFile, psGffFile = sGffFile)

  ### # do the counting, get the bam files from input
  bamFile = input$getFullPaths("BAM")

  ### # determine extension for count files
  sCountfileExt <- 'count'
  if (ezIsSpecified(param$countfile_ext))
    sCountfileExt <- param$countfile_ext
  ### # call counting routine
  vCountFiles <- sapply(bamFile, runCountSingleBam, sGffFile, sCountfileExt)
  
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodDEXSeqCounting(input=NA, output=NA, param=NA)
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


#' Convert annotation file from GTF format to GFF
#' 
#' @description 
#' \code{convertGtfToGff} converts an annotation file from 
#' the GTF format into the GFF format which is required 
#' by the package \code{DEXSeq}. Input file name and the 
#' name of the file to be generated are both given as 
#' function parameters. The conversion is done by a python 
#' script that is given by the content of \code{DEXSEQ_PREPARE} 
#' which is either taken as a global variable or from the 
#' result of function \code{lGetPyScriptPaths}
#' 
#' @param psGtfFile   name of the GTF annotation file
#' @param psGffFile   name of the GFF file to be generated
convertGtfToGff <- function(psGtfFile, psGffFile) {
  cat(" * Converting GTF to GFF ...\n")
  ### # check whether the path exists
  if (!exists("DEXSEQ_PREPARE")) {
    DEXSEQ_PREPARE <- lGetPyScriptPaths()$DEXSEQ_PREPARE
  }
  sPyConvCmd <- paste(DEXSEQ_PREPARE, psGtfFile, psGffFile)
  ezSystem(sPyConvCmd)
  cat("  ==> created: ", psGffFile, "\n")
  invisible(TRUE)
}

#' Run counts for a single BAM file
#' 
runCountSingleBam <- function(psBamFile, psGffFile, psCountfileExt){
  sSamCmd <- paste(SAMTOOLS, "view -h", psBamFile)
  ### # run counting on sam file
  sCountBaseFn <- gsub("bam$", psCountfileExt, basename(psBamFile))
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

