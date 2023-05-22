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
                  appDefaults <<- rbind(controlSeqs=ezFrame(Type="charVector",
                                                            DefaultValue="",
                                                            Description="control sequences to add"))
                }
              )
  )

ezMethodSpaceRanger <- function(input=NA, output=NA, param=NA){

  sampleName <- input$getNames()
  
  sampleDirs <- getFastqDirs(input, "RawDataDir", sampleName)
  sampleNameFQ <- sub('.tar', '', basename(sampleDirs))
  finalSampleName <-sampleName
  
  # decompress tar files if they are in tar format
  if (all(grepl("\\.tar$", sampleDirs)))
    sampleDirs <- deCompress(sampleDirs)
  
  sampleDirs <- normalizePath(sampleDirs)
  sampleDir <- paste(sampleDirs, collapse = ",")
  spaceRangerFolder <- str_sub(sampleName, 1, 45) %>% str_c("-spaceRanger")
  spaceRangerFolder <- gsub('\\.', '_', spaceRangerFolder)
  if(sampleName != sampleNameFQ){
      sampleName <- sampleNameFQ
  }
  
  inputCols <- colnames(input$meta)
  
  refDir <- getCellRangerGEXReference(param)
  cmd <- paste("spaceranger count", paste0("--id=", spaceRangerFolder),
               paste0("--transcriptome=", refDir),
               paste0("--fastqs=", sampleDir),
               paste0("--sample=", sampleName),
               paste0("--localmem=", param$ram),
               paste0("--localcores=", param$cores),
               paste0("--slide=", input$getColumn("Slide")),
               paste0("--area=", input$getColumn("Area")))
  
  if('Image' %in% inputCols & grepl('tif$|tiff$|jpeg$|jpg$',input$getFullPaths("Image"))){
      cmd <- paste(cmd, paste0("--image=", input$getFullPaths("Image")))
  }
  
  if('CytaImage' %in% inputCols){
      cmd <- paste(cmd, paste0("--cytaimage=", input$getFullPaths("CytaImage")))
  }
  
  if(param$probesetFile!=''){
    myFile <- file.path('/srv/GT/databases/10x_Probesets/Visium',param$probesetFile)
    outputFile <- sub('.csv','_filtered.csv', basename(myFile))
    maxHeaderLine <- max(grep('#', readLines(myFile)))
    headerSection <- readLines(myFile, n = maxHeaderLine)
    headerSection[grep('reference_genome', headerSection)] = paste0('#reference_genome=',basename(refDir))
    probeInfo <- ezRead.table(myFile, sep = ',', row.names = NULL, skip = maxHeaderLine)
    annotation <- ezRead.table(file.path(refDir, 'star', 'geneInfo.tab'), row.names = NULL, skip = 1, header = FALSE)
    intersectionGenes <- intersect(annotation$V1, probeInfo$gene_id)
    probeInfo <- probeInfo[probeInfo$gene_id %in% intersectionGenes, ]
    writeLines(headerSection, outputFile)
    ezWrite.table(probeInfo, outputFile, sep = ',', row.names = FALSE, append = TRUE)
    cmd <- paste(cmd, paste0("--probe-set=", file.path(getwd(), outputFile)))
  }
  
  tryCatch(
    {
    json_paths <- input$getFullPaths("loupe-alignment")
        if(grepl('json$', basename(json_paths))){
            cmd <- paste0(cmd, " --loupe-alignment=", json_paths)
        }
    },
    error=function(e) {
      return()
    })
    
  if(ezIsSpecified(param$cmdOptions)){
    cmd <- paste(cmd, param$cmdOptions)
  }
  ezSystem(cmd)
  
  
  unlink(basename(sampleDirs), recursive=TRUE)
  file.rename(file.path(spaceRangerFolder, "outs"),  finalSampleName)
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