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
                                                            Description="control sequences to add"),
                                        keepAlignment = ezFrame(
                                            Type = "logical",
                                            DefaultValue = FALSE,
                                            Description = "keep bam/cram file produced by SpaceRanger"
                                        ),
                                        darkImage = ezFrame(
                                            Type = "logical",
                                            DefaultValue = FALSE,
                                            Description = "use dark image mode"
                                        )
                                        )
                }
              )
  )

ezMethodSpaceRanger <- function(input=NA, output=NA, param=NA){
  library(Seurat)
  sampleName <- input$getNames()
  
  sampleDirs <- getFastqDirs(input, "RawDataDir", sampleName)
  sampleNameFQ <- sub('.tar', '', basename(sampleDirs))
  finalSampleName <-sampleName
  
  # extract tar files if they are in tar format
  if (all(grepl("\\.tar$", sampleDirs))) {
      sampleDirs <- tarExtract(sampleDirs, prependUnique=TRUE)
  } else {
      stop("Require inputs to be provided in .tar files.")
  }
  
  sampleDirs <- normalizePath(sampleDirs)
  sampleDir <- paste(sampleDirs, collapse = ",")
  spaceRangerFolder <- str_sub(sampleName, 1, 45) %>% str_c("-spaceRanger")
  spaceRangerFolder <- gsub('\\.', '_', spaceRangerFolder)
  
  if(length(sampleNameFQ) == 1){
    if(sampleName != sampleNameFQ){
      sampleName <- sampleNameFQ
    }
  } else if(any(sampleNameFQ != sampleName)){
      #2.1 Fix FileNames
      cwd <- getwd()
      sampleNameFQ <- file.path(strsplit(sampleDir, ',')[[1]], sampleNameFQ)
      for (fileLevelDir in sampleNameFQ) {
              setwd(fileLevelDir)
              cmd <- paste('rename', 
                           paste0('s/^', basename(fileLevelDir),'/',sampleName, '/g'), 
                           paste0(basename(fileLevelDir),'*.gz'))
              ezSystem(cmd)
          }
          setwd(cwd)
      }
  
  inputCols <- colnames(input$meta)
  
  refDir <- getCellRangerGEXReference(param)
  cmd <- paste("spaceranger count", paste0("--id=", spaceRangerFolder),
               paste0("--transcriptome=", refDir),
               paste0("--fastqs=", sampleDir),
               paste0("--sample=", sampleName),
               paste0("--localmem=", param$ram),
               paste0("--localcores=", param$cores),
               if(grepl('^3', basename(param$SpaceRangerVersion))){paste0("--create-bam true")})
  
  if('Image' %in% inputCols && grepl('btf$|tif$|tiff$|jpeg$|jpg$',input$meta['Image']$Image)){
      if(!param$darkImage){
      cmd <- paste(cmd, paste0("--image=", input$getFullPaths("Image")))
      } else {
          if(file.exists(input$meta['Image']$Image)){
            cmd <- paste(cmd, paste0("--darkimage=", input$getFullPaths("Image")))
          } else {
              images <- paste(file.path(param$dataRoot,unlist(strsplit(input$meta['Image']$Image, ',')[[1]])),collapse=',')
              cmd <- paste(cmd, paste0("--darkimage=", images))
          }
    }
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
    if (ezIsSpecified(param$customProbesFile) && param$customProbesFile != '') {
        customProbes <- ezRead.table(file.path(param$dataRoot, param$customProbesFile), sep=',', row.names=NULL) %>%
            mutate(gene_id=ifelse(startsWith(gene_id, "Gene_"), gene_id, paste0("Gene_", gene_id)),
                   probe_id=ifelse(startsWith(probe_id, "Gene_"), probe_id, paste0("Gene_", probe_id)))
        probeInfo <- bind_rows(list(probeInfo, customProbes))
    }
    writeLines(headerSection, outputFile)
    ezWrite.table(probeInfo, outputFile, sep = ',', row.names = FALSE, append = TRUE)
    cmd <- paste(cmd, paste0("--probe-set=", file.path(getwd(), outputFile)))
  }
  
  if(param$panelFile!=''){
      cmd <- sub("--fastqs=.*--localmem=","--localmem=", cmd)
      myFile <- file.path('/srv/GT/databases/10x/Visium/panels',param$panelFile)
      cmd <- paste(cmd, "--feature-ref", myFile)
      
      panelSampleDirs <- getFastqDirs(input, "PanelRawDataDir", sampleName)
      panelSampleNameFQ <- sub('.tar', '', basename(panelSampleDirs))
      
      panelSampleDirs <- tarExtract(panelSampleDirs, prependUnique=TRUE)
      panelSampleDirs <- normalizePath(panelSampleDirs)
      
      librariesDS <- data.frame(fastqs = c(sampleDir, panelSampleDirs), sample = c(sampleName, panelSampleNameFQ), library_type = c('Gene Expression','Antibody Capture'))
      write_csv(librariesDS, 'libraries.csv')      
      cmd <- paste(cmd, "--libraries=libraries.csv")
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
  
  if(!grepl('--unknown-slide', cmd)){
    cmd <- paste(cmd,
           paste0("--slide=", input$getColumn("Slide")),
           paste0("--area=", input$getColumn("Area")))
  }
  
  ezSystem(cmd)
  
  
  unlink(basename(sampleDirs), recursive=TRUE)
  file.rename(file.path(spaceRangerFolder, "outs"),  finalSampleName)
  unlink(spaceRangerFolder, recursive=TRUE)
  
  if(ezIsSpecified(param$controlSeqs)){
    unlink(refDir, recursive = TRUE)
  }
  
  # Optional removal of the bam files
  if(!param$keepAlignment){
      print(ezSystem('find . -name "*.bam" -type f'))
      ezSystem('find . -name "*.bam" -type f -delete')
      ezSystem('find . -name "*.bam.bai" -type f -delete')
  } else {
      if(param$secondRef == ''){
        setwd(finalSampleName)
         refDir <- param$ezRef["refFastaFile"]
        bamFile <- 'possorted_genome_bam.bam'
        out <- tryCatch(ezSystem(paste('samtools view', '-T', refDir, '-@', param$cores, '-o', sub('.bam$', '.cram', bamFile), '-C', bamFile)), error = function(e) NULL)
        system('rm possorted_genome_bam.bam')
        setwd('..')
      }
  }
  
  cmDir <- file.path(finalSampleName, 'filtered_feature_bc_matrix')
  if(dir.exists(cmDir)){
    cts <- Read10X(cmDir, gene.column = 1)
  } else { ##visium HD data
      setwd(finalSampleName)
      system('ln -s binned_outputs/square_016um/filtered_feature_bc_matrix .')
      setwd('..')
      cmDir <- file.path(finalSampleName, 'binned_outputs/square_016um/filtered_feature_bc_matrix')
      cts <- Read10X(cmDir, gene.column = 1)
  }
  if(is.list(cts)){
      cts <- cts[['Gene Expression']]
      bulkData <- apply(cts,1,sum)
  } else {
    bulkData <- rowSums(data.frame(cts))
  }
  bulkData <- data.frame(Identifier = names(bulkData), matchCounts = bulkData)
  countFile <- paste0(finalSampleName,'-counts.txt')
  ezWrite.table(bulkData, file.path(finalSampleName, countFile), row.names = FALSE)
  return("Success")
}

getFastqDirs <- function(input, column, sampleName) {
  fastqDirs <- strsplit(input$getColumn(column), ",")[[sampleName]]
  fastqDirs <- file.path(input$dataRoot, fastqDirs)
  return(fastqDirs)
}
