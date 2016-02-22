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
  sapply(as.vector(bamFile), runCountSingleBam, sGffFile)
  
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
                  # appDefaults <<- rbind()
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
  sPyScrPath <- system.file(package = "DEXSeq", "python_scripts")
  sPyScrFn <- "dexseq_prepare_annotation.py"
  sGffFn <- gsub("gtf$", "gff", sGtfFn)
  sPyConvCmd <- paste(HTSEQ_PREFIX, file.path(sPyScrPath,sPyScrFn), sGtfFn, sGffFn)
  ezSystem(sPyConvCmd)
  cat("  ... created", sGffFn, "\n")
  setwd(sCurWd)
  invisible(TRUE)
}

#' Run counts for a single BAM file
#' 
runCountSingleBam <- function(psBamFile, psGffFile){
  sBamBaseFn <- basename(psBamFile)
  sSamBaseFn <- gsub("bam$", "sam", sBamBaseFn)
  ### # convert bamfile to sam file
  sSamCmd <- paste(SAMTOOLS, "view -h -o", sSamBaseFn, psBamFile)
  ezSystem(sSamCmd)
  ### # run counting on sam file
  sCountBaseFn <- gsub("sam$", "count", sSamBaseFn)
  sPyScrCountFn <- "dexseq_count.py"
  sPyScrPath <- system.file(package = "DEXSeq", "python_scripts")
  sPyCountCmd <- paste(HTSEQ_PREFIX, file.path(sPyScrPath,sPyScrCountFn), psGffFile, sSamBaseFn, sCountBaseFn)
  ezSystem(sPyCountCmd)
  ### # clean up and remove sam file
  if (file.exists(sSamBaseFn))
    unlink(sSamBaseFn)
}