###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppCellBender <-
  setRefClass("EzAppCellBender", contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCellBender
                  name <<- "EzAppCellBender"
                  appDefaults <<- rbind(
                    cmdOptions=ezFrame(Type="character", DefaultValue="", Description="for -expected-cells and --total-droplets-included"),
                    gpuMode = ezFrame(
                      Type = "logical",
                      DefaultValue = FALSE,
                      Description = "run with cuda option"
                    )
                  )
                }
              )
  )

ezMethodCellBender <- function(input = NA, output = NA, param = NA) {
  sampleName = input$getNames()
  setwdNew(sampleName)
  
  # Initialize cmDir before tryCatch
  cmDir <- NULL
  
  # Try to get path for single modality case first
  tryCatch({
    if ("UnfilteredCountMatrix" %in% input$colNames) {
      cmDir <- input$getFullPaths("UnfilteredCountMatrix")
    } else {
      # Multi-modal case - construct path from ResultDir
      resultPath <- input$getColumn("ResultDir")
      cmDir <- file.path(param$dataRoot, resultPath, 
                         "multi/count/raw_feature_bc_matrix")
    }
  }, error = function(e) {
    stop(sprintf("Failed to construct path: %s", e$message))
  })
  
  if(is.null(cmDir)) {
    stop("Failed to get valid path for count matrix")
  }
  
  inputFile <- paste0(cmDir, '.h5')
  if(!file.exists(inputFile)) {
    stop('RawCountMatrix missing! Unsupported InputFormat.')
  }
  
  cmd <- paste("cellbender remove-background",
               "--input", inputFile,
               "--output cellbender.h5")
  
  if(param$cmdOptions!=''){
    cmd <- paste(cmd, cmdOptions)
  }
  
  if(param$gpuMode) {
    cmd <- paste(cmd, "--cuda")
  } else {
    cmd <- paste(cmd, '--cpu-threads', param$cores)
  }
  system(cmd)
  
  ##Post processing for Seurat
  cmd <- paste("ptrepack --complevel 5 cellbender_filtered.h5:/matrix cellbender_filtered_seurat.h5:/matrix")
  system(cmd)
  cmd <- paste("ptrepack --complevel 5 cellbender.h5:/matrix cellbender_raw_seurat.h5:/matrix")
  system(cmd)
  ##Clean Up:
  system('rm ckpt.tar.gz cellbender_posterior.h5 cellbender.h5 cellbender_filtered.h5')
  return("Success")
}
