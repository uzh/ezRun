###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSpaceRanger <-
  setRefClass("EzAppSpaceRanger",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSpaceRanger
                  name <<- "EzAppSpaceRanger"
                  appDefaults <<- rbind(scMode=ezFrame(Type="character",
                                                       DefaultValue="SC",
                                                       Description="Single cell or single nuclei?"),
                                        controlSeqs=ezFrame(Type="charVector",
                                                            DefaultValue="",
                                                            Description="control sequences to add"))
                }
              )
  )

ezMethodSpaceRanger <- function(input=NA, output=NA, param=NA){
  sampleName <- input$getNames()
  sampleDirs <- strsplit(input$getColumn("RawDataDir"), ",")[[sampleName]]
  sampleDirs <- file.path(input$dataRoot, sampleDirs)
  if(all(grepl("\\.tar$", sampleDirs))){
    # This is new .tar folder
    lapply(sampleDirs, untar)
    sampleDirs <- sub("\\.tar$", "", basename(sampleDirs))
  }
  sampleDir <- paste(sampleDirs, collapse=",")
  spaceRangerFolder <- paste0(sampleName, "-spaceRanger")

  refDir <- getCellRangerGEXReference(param)
  cmd <- paste(SPACERANGER, "count", paste0("--id=", spaceRangerFolder),
               paste0("--transcriptome=", refDir),
               paste0("--fastqs=", sampleDir),
               paste0("--sample=", sampleName),
               paste0("--localmem=", param$ram),
               paste0("--localcores=", param$cores),
               paste0("--image=", input$getFullPaths("image")),
               paste0("--slide=", input$getColumn("slide")),
               paste0("--area=", input$getColumn("area")))
  if(ezIsSpecified(param$cmdOptions)){
    cmd = paste(cmd, param$cmdOptions)
  }
  ezSystem(cmd)
  
  unlink(basename(sampleDirs), recursive=TRUE)
  file.rename(file.path(spaceRangerFolder, "outs"),  sampleName)
  unlink(spaceRangerFolder, recursive=TRUE)
  
  if(ezIsSpecified(param$controlSeqs)){
    unlink(refDir, recursive = TRUE)
  }
  
  require(DropletUtils)
  require(Matrix)
  require(readr)
  countMatrixFn <- list.files(path=file.path(sampleName, 'filtered_feature_bc_matrix'),
                              pattern="\\.mtx(\\.gz)*$", recursive=TRUE,
                              full.names=TRUE)
  sce <- read10xCounts(dirname(countMatrixFn), col.names=TRUE)
    
  cellPhase <- getCellCycle(sce, param$refBuild)
  write_tsv(cellPhase,
            path=file.path(dirname(countMatrixFn), "CellCyclePhase.txt"))
  
  return("Success")
}
