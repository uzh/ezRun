###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCMultipleSamples <-
  setRefClass("EzAppSCMultipleSamples",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCMultipleSamples
                  name <<- "EzAppSCMultipleSamples"
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

ezMethodSCMultipleSamples = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(rlist)
  library(HDF5Array)
  library(scanalysis)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  

  if(input$getLength() > param$maxSamplesSupported){
    stop(paste("It only works for", param$maxSamplesSupported, "at most"))
  }
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("SC Cluster Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  #the individual sce objects can be in hdf5 format (for new reports) or in rds format (for old reports)
  sceURLs <- input$getColumn("Static Report")
  filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(sceURLs)), "sce_h5")
  filePath_course <- file.path(paste0("/srv/GT/analysis/course_sushi/public/projects/", input$getColumn("Report"), "/sce_h5"))
  
  #In case it is in hdf5 format the Seurat object is not stored in the metadata slot, so we have to build it and store it there, since the 
  #clustering functions below work with sce objects and take the seurat object from them.
  if(!file.exists(filePath)) 
    filePath <- filePath_course
  
  sceList <- lapply(filePath,loadHDF5SummarizedExperiment)
  names(sceList) <- names(sceURLs)
  scDataList <- lapply(sceList, function(sce) {
      sce = swapAltExp(sce, "RNA", withColData = FALSE)
      colData(sce)[, grep("SCT", colnames(colData(sce)))] = NULL #remove previous clustering done on SCT assay
      CreateSeuratObject(counts=counts(sce),meta.data=data.frame(colData(sce)))})
  
  
  pvalue_allMarkers <- 0.05
  nrSamples <- length(sceList)
  
  if(ezIsSpecified(param$chosenClusters)){
    for(eachSample in names(param$chosenClusters)){
      chosenCells <- names(Idents(metadata(sceList[[eachSample]])$scData))[Idents(metadata(sceList[[eachSample]])$scData) %in% param$chosenClusters[[eachSample]]]
      sceList[[eachSample]] <- sceList[[eachSample]][, chosenCells]
      metadata(sceList[[eachSample]])$scData <-
        SubsetData(metadata(sceList[[eachSample]])$scData,
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
  
  #positive cluster markers
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  
  #we do cell type identification using AUCell and SingleR
  cells_AUC <- NULL
  singler.results <- NULL
  #cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if(species == "Human" | species == "Mouse") {
    cells_AUC = cellsLabelsWithAUC(scData, species, param$tissue)
    singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "counts"), Idents(scData), species)
  }
  
  #Convert scData to Single Cell experiment Object
  #TODO: remove unnecesary dietseurat call when the bug in Seurat is fixed
  scData_diet = DietSeurat(scData, dimreducs = c("pca", "tsne", "umap", "umap_noCorrected"))
  sce <- scData_diet %>% seurat_to_sce(default_assay = "SCT")
  metadata(sce)$PCA_stdev <- Reductions(scData_diet, "pca")@stdev
  metadata(sce)$cells_AUC <- cells_AUC
  metadata(sce)$singler.results <- singler.results
  metadata(sce)$output <- output
  metadata(sce)$param <- param
  metadata(sce)$param$name <- paste(param$name, paste(input$getNames(), collapse=", "), sep=": ")
  
  geneMeans <- geneMeansCluster(sce)
  #Save some results in external files 
  dataFiles = saveExternalFiles(list(pos_markers=posMarkers, gene_means=as_tibble(as.data.frame(geneMeans), rownames="gene_name")))
  
  library(HDF5Array)
  saveHDF5SummarizedExperiment(sce, dir="sce_h5")
  
  makeRmdReport(dataFiles=dataFiles, rmdFile = "SCMultipleSamples.Rmd", reportTitle = metadata(sce)$param$name) 
  return("Success")
  
}






