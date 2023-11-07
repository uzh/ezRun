###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratCombine <-
  setRefClass("EzAppScSeuratCombine",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScSeuratCombine
                  name <<- "EzAppScSeuratCombine"
                  appDefaults <<- rbind(npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=30, 
                                                    Description="The maximal dimensions to use for reduction"),
                                        pcGenes = ezFrame(
                                          Type = "charVector",
                                          DefaultValue = "",
                                          Description = "The genes used in supvervised clustering"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        batchCorrection=ezFrame(Type="logical", 
                                                                DefaultValue="TRUE",
                                                                Description="Perform batch correction."),
                                        integrationMethod=ezFrame(Type="character", 
                                                                  DefaultValue="Classic", 
                                                                  Description="Choose integration method in Seurat (Classic or RPCA)"),
                                        enrichrDatabase=ezFrame(Type = "charVector", 
                                                                DefaultValue = "", 
                                                                Description="enrichR databases to search"),
                                        computePathwayTFActivity=ezFrame(Type="logical", 
                                                                 DefaultValue="TRUE",
                                                                 Description="Whether we should compute pathway and TF activities."),
                                        SCT.regress.CellCycle=ezFrame(
                                          Type = "logical", 
                                          DefaultValue = FALSE,
                                          Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
                                        ),
                                        DE.method=ezFrame(Type="charVector", 
                                                          DefaultValue="wilcox", 
                                                          Description="Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"),
                                        DE.regress=ezFrame(Type="charVector", 
                                                           DefaultValue="Batch", 
                                                           Description="Variables to regress out if the test LR is chosen"))
                }
              )
  )

ezMethodScSeuratCombine = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(rlist)
  library(HDF5Array)
  library(scanalysis)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(AUCell)
  library(enrichR)
  library(decoupleR)
  library(Azimuth)
  library(BiocParallel)
  
  BPPARAM <- MulticoreParam(workers = param$cores)
  register(BPPARAM)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  #the individual sce objects can be in hdf5 format (for new reports) or in rds format (for old reports)
  filePath <- file.path("/srv/gstore/projects", input$getColumn("SC Cluster Report"), 'scData.rds')
  filePath_course <- file.path("/srv/GT/analysis/course_sushi/public/projects", input$getColumn("SC Cluster Report"), 'scData.rds')
  
  if(!file.exists(filePath[1])) 
    filePath <- filePath_course
  names(filePath) <- input$getNames()
  
  # Load the data and prepare metadata for integration
  scDataList <- lapply(names(filePath), function(sm) {
    scData <- readRDS(filePath[sm])
    scData$Sample <- sm
    # If we have new information in the Condition column, add it to the dataset
    if (all(scData$Condition == "NA" | scData$Condition == "") || 
           (ezIsSpecified(param$overwriteCondition) && as.logical(param$overwriteCondition))) {
      if (any(startsWith(colnames(input$meta), "Condition"))) {
        scData$Condition <- unname(input$getColumn("Condition")[sm])
      } else {
        scData$Condition <- scData$Sample
      }
    }
    # Harmony will complain if the Condition is the same across all samples
    if (param$integrationMethod == "Harmony" && 
        length(unique(input$meta$`Condition`)) == 1) {
      scData$Condition <- scData$Sample
    } else if (ezIsSpecified(param$STACASAnnotationFile)) {
      clusterAnnoFn <- file.path(param$dataRoot, param$STACASAnnotationFile)
      clusterAnno <- readxl::read_xlsx(clusterAnnoFn) %>% 
        as_tibble() %>%
        dplyr::select(1:3) %>% # remove all other columns
        dplyr::rename(c("Sample"=1, "Cluster"=2, "ClusterLabel"=3)) %>%
        dplyr::filter(Sample == sm)
      labelMap <- as.character(clusterAnno$ClusterLabel)
      names(labelMap) <- as.character(clusterAnno$Cluster)
      scData$stacasLabelColumn <- unname(labelMap[as.character(scData$seurat_clusters)])
    }
    # Also add the other factors in the input dataset to the objects
    if (!is.null(param$harmonyFactors)) {
      harmonyFactors <- str_split(param$harmonyFactors, ",", simplify=TRUE)[1,]
      metaFactorNames <- paste0("har_", harmonyFactors) %>% str_replace(., " ", ".")
      names(metaFactorNames) <- harmonyFactors
      for (hf in harmonyFactors) {
        scData[[metaFactorNames[hf]]] <- unname(input$getColumn(hf)[sm])
      }
    }
    # Rename the cells and add original sample-level clusters back in
    scData <- RenameCells(scData, new.names = paste0(scData$Sample, "-", colnames(scData)))
    scData$sample_seurat_clusters <- paste0(scData$Sample, "-", sprintf("%02d", scData$seurat_clusters))
    return(scData)
  })
  
  # perform all of the analysis
  results <- seuratIntegrateDataAndAnnotate(scDataList, input, output, param)
  
  # save the markers
  writexl::write_xlsx(results$markers, path="posMarkers.xlsx")
  
  # Save some results in external files
  reportTitle <- 'SCReport - MultipleSamples based on Seurat'
  makeRmdReport(param=param, output=output, scData=results$scData, 
                enrichRout=results$enrichRout, TFActivity=results$TFActivity, 
                pathwayActivity=results$pathwayActivity, aziResults=results$aziResults,
                cells.AUC=results$cells.AUC, singler.results=results$singler.results,
                rmdFile = "ScSeuratCombine.Rmd", reportTitle = reportTitle) 
  return("Success")
  
}

seuratIntegrateDataAndAnnotate <- function(scDataList, input, output, param) {
  pvalue_allMarkers <- 0.05
  
  if(ezIsSpecified(param$chosenClusters)){
    for(eachSample in names(param$chosenClusters)){
      chosenCells <- names(Idents(scDataList[[eachSample]]))[Idents(scDataList[[eachSample]]) %in% param$chosenClusters[[eachSample]]]
      scDataList[[eachSample]] <- scDataList[[eachSample]][, chosenCells]
    }
  }
  
  scData_noCorrected <- cellClustNoCorrection(scDataList, param)
  if (param$batchCorrection) {
    scData_corrected = cellClustWithCorrection(scDataList, param)
    #in order to compute the markers we switch again to the original assay
    DefaultAssay(scData_corrected) <- "SCT"
    scData <- scData_corrected
  } else {
    scData = scData_noCorrected
  }
  scData@reductions$tsne_noCorrected <- Reductions(scData_noCorrected, "tsne")
  scData@reductions$umap_noCorrected <- Reductions(scData_noCorrected, "umap")
  scData@meta.data$ident_noCorrected <- Idents(scData_noCorrected)
  scData <- PrepSCTFindMarkers(scData)
  
  # get annotation information
  anno <- getSeuratMarkersAndAnnotate(scData, param)
  
  return(list(scData=scData, 
              markers=anno$markers,
              enrichRout=anno$enrichRout, 
              pathwayActivity=anno$pathwayActivity, 
              TFActivity=anno$TFActivity,
              cells.AUC=anno$cells.AUC,
              singler.results=anno$singler.results,
              aziResults=anno$aziResults))
}
