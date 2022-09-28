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
                                        SCT.regress=ezFrame(Type="character", 
                                                            DefaultValue="none", 
                                                            Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."),
                                        DE.method=ezFrame(Type="charVector", 
                                                          DefaultValue="wilcox", 
                                                          Description="Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"),
                                        DE.regress=ezFrame(Type="charVector", 
                                                           DefaultValue="Batch", 
                                                           Description="Variables to regress out if the test LR is chosen"),
                                        maxSamplesSupported=ezFrame(Type="numeric", 
                                                              DefaultValue=5, 
                                                              Description="Maximum number of samples to compare"))
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

  library(BiocParallel)
  
  BPPARAM <- MulticoreParam(workers = param$cores)
  register(BPPARAM)
  

  if(input$getLength() > param$maxSamplesSupported){
    stop(paste("It only works for", param$maxSamplesSupported, "at most"))
  }
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  #the individual sce objects can be in hdf5 format (for new reports) or in rds format (for old reports)
  filePath <- file.path("/srv/gstore/projects", input$getColumn("SC Cluster Report"), 'scData.rds')
  filePath_course <- file.path("/srv/GT/analysis/course_sushi/public/projects", input$getColumn("SC Cluster Report"), 'scData.rds')
  
  if(!file.exists(filePath[1])) 
    filePath <- filePath_course
  
  scDataList <- lapply(filePath, readRDS)
  names(scDataList) <- names(input$getColumn("Static Report"))
  
  pvalue_allMarkers <- 0.05
  nrSamples <- length(scDataList)
  
  if(ezIsSpecified(param$chosenClusters)){
    for(eachSample in names(param$chosenClusters)){
      chosenCells <- names(Idents(scDataList[[eachSample]]))[Idents(scDataList[[eachSample]]) %in% param$chosenClusters[[eachSample]]]
      scDataList[[eachSample]] <- scDataList[[eachSample]][, chosenCells]
       scDataList[[eachSample]] <-
        SubsetData(scDataList[[eachSample]],
                   ident.use=param$chosenClusters[[eachSample]])
    }
  }
  
  scData_noCorrected <- cellClustNoCorrection(scDataList, param)
  scData = scData_noCorrected
  if (param$batchCorrection) {
    scData_corrected = cellClustWithCorrection(scDataList, param)
    #in order to compute the markers we switch again to the original assay
    DefaultAssay(scData_corrected) <- "SCT"
    scData <- scData_corrected
  }
  scData@reductions$tsne_noCorrected <- Reductions(scData_noCorrected, "tsne")
  scData@reductions$umap_noCorrected <- Reductions(scData_noCorrected, "umap")
  scData@meta.data$ident_noCorrected <- Idents(scData_noCorrected)
  scData <- PrepSCTFindMarkers(scData)
  
  #positive cluster markers
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  
  #we do cell type identification using AUCell and SingleR
  cells_AUC <- NULL
  singler.results <- NULL
  #cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if(species == "Human" | species == "Mouse") {
    cells.AUC = cellsLabelsWithAUC(GetAssayData(scData, "counts"), species, param$tissue, BPPARAM=BPPARAM)
    singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "counts"), Idents(scData), species)
    saveRDS(cells.AUC, file="cells.AUC.rds")
    saveRDS(singler.results, file="singler.results.rds")
  }
  
  geneMeans <- geneMeansCluster(scData)
  
  saveRDS(param, file="param.rds")
  saveRDS(output, file="output.rds")
  saveRDS(scData, file = "scData.rds")
  
  # Save some results in external files
  dataFiles <- saveExternalFiles(list(pos_markers = posMarkers, gene_means = as_tibble(as.data.frame(geneMeans), rownames = "gene_name")))
  
  makeRmdReport(dataFiles=dataFiles, rmdFile = "ScSeuratCombine.Rmd", reportTitle = param$name) 
  return("Success")
  
}






