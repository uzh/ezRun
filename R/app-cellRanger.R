###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##ToD: --localmem, --localcores, use opt-Parameters
ezMethodCellRanger = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  sampleName = input$getNames()
  sampleDir = input$getFullPaths("RawDataDir")
  
  refDir <- dirname(param$ezRef["refFeatureFile"])
  refDirs <- list.files(path=refDir, pattern="^10X_Ref", full.names = TRUE)
  if(length(refDirs) == 0){
    stop("No 10X_Ref folder found in", refDir)
  }
  if(length(refDirs) > 1){
    warning("Multiple 10X_Ref folders in ", refDir)
  }
  refDir <- refDirs[1]
  
  cmd = paste(CELLRANGER,"count", paste0("--id=", sampleName),
              paste0("--transcriptome=", refDir),
              paste0("--fastqs=", sampleDir),
              paste0("--sample=", sampleName),
              paste0("--localmem=",param$ram),
              paste0("--localcores=",param$cores))
  if(opt!=''){
    cmd = paste(cmd, opt)
  }
  ezSystem(cmd)
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodCellRanger(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppCellRanger <-
  setRefClass("EzAppCellRanger",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCellRanger
                  name <<- "EzAppCellRanger"
                }
              )
  )
