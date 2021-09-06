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
  sampleDirs <- getFastqDirs(input, "RawDataDir",sampleName)
  
  # decompress tar files if they are in tar format
  if (all(grepl("\\.tar$", sampleDirs)))
    sampleDirs <- deCompress(sampleDirs)
  
  sampleDirs <- normalizePath(sampleDirs)
  sampleDir <- paste(sampleDirs, collapse = ",")
  spaceRangerFolder <- str_sub(sampleName, 1, 45) %>% str_c("-spaceRanger")

  refDir <- getCellRangerGEXReference(param)
  cmd <- paste("spaceranger count", paste0("--id=", spaceRangerFolder),
               paste0("--transcriptome=", refDir),
               paste0("--fastqs=", sampleDir),
               paste0("--sample=", sampleName),
               paste0("--localmem=", param$ram),
               paste0("--localcores=", param$cores),
               paste0("--image=", input$getFullPaths("image")),
               paste0("--slide=", input$getColumn("slide")),
               paste0("--area=", input$getColumn("area")))
  tryCatch(
    {
    json_paths <- input$getFullPaths("loupe-alignment")
    cmd <- paste0(cmd, " --loupe-alignment=", json_paths)
    },
    error=function(e) {
      return()
    })
    
  if(ezIsSpecified(param$cmdOptions)){
    cmd <- paste(cmd, param$cmdOptions)
  }
  ezSystem(cmd)
  
  unlink(basename(sampleDirs), recursive=TRUE)
  file.rename(file.path(spaceRangerFolder, "outs"),  sampleName)
  unlink(spaceRangerFolder, recursive=TRUE)
  
  if(ezIsSpecified(param$controlSeqs)){
    unlink(refDir, recursive = TRUE)
  }
  
  return("Success")
}
getFastqDirs <- function(input, column, sampleName) {
  fastqDirs <- strsplit(input$getColumn(column), ",")[[sampleName]]
  fastqDirs <- file.path(input$dataRoot, fastqDirs)
  return(fastqDirs)
}

deCompress = function(fastqDirs){
  fastqDirs = sapply(fastqDirs, function(scTar){
    targetDir = basename(dirname(scTar))
    untar(scTar, exdir = targetDir, tar=system("which tar", intern=TRUE))
    return(targetDir)
  })
  return(fastqDirs)
}