###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodCellRanger(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppCellRangerATAC <-
  setRefClass("EzAppCellRangerATAC",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCellRangerATAC
                  name <<- "EzAppCellRangerATAC"
                }
              )
  )

ezMethodCellRangerATAC <- function(input=NA, output=NA, param=NA){
  sampleName = input$getNames()
  sampleDirs = strsplit(input$getColumn("RawDataDir"), ",")[[sampleName]]
  sampleDirs <- file.path(input$dataRoot, sampleDirs)
  if(all(grepl("\\.tar$", sampleDirs))){
    # This is new .tar folder
    lapply(sampleDirs, untar)
    sampleDirs <- sub("\\.tar$", "", basename(sampleDirs))
  }
  sampleDir <- paste(sampleDirs, collapse=",")
  cellRangerFolder = paste0(sampleName, "-cellRanger")
  
  refDir <- getCellRangerATACReference(param)
  message("Using the reference: ", refDir)
  
  cmd <- paste(CELLRANGERATAC, "count", paste0("--id=", cellRangerFolder),
               paste0("--reference=", refDir),
               paste0("--fastqs=", sampleDir),
               paste0("--sample=", sampleName),
               paste0("--localmem=", param$ram),
               paste0("--localcores=", param$cores))
  
  if(ezIsSpecified(param$cmdOptions)){
    cmd = paste(cmd, param$cmdOptions)
  }
  ezSystem(cmd)
  
  unlink(basename(sampleDirs), recursive=TRUE)
  file.rename(file.path(cellRangerFolder, "outs"),  sampleName)
  unlink(cellRangerFolder, recursive=TRUE)
  
  return("Success")
}

getCellRangerATACReference <- function(param){
  refDir <- file.path(param$ezRef["refBuildDir"], "Sequence")
  refDir <- list.files(refDir, pattern="cellranger-atac", 
                       full.names=TRUE, recursive=FALSE)
  if(length(refDir) == 0L){
    stop("Cell Ranger ATAC compatible genome reference is not available. Please download from https://support.10xgenomics.com/single-cell-atac/software/downloads/latest!")
  }
  refDir <- refDir[1]
  return(refDir)
}