###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCMultipleSamplesOneGroup <-
  setRefClass("EzAppSCMultipleSamplesOneGroup",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCMultipleSamplesOneGroup
                  name <<- "EzAppSCMultipleSamplesOneGroup"
                  appDefaults <<- rbind(npcs=ezFrame(Type="numeric", 
                                                     DefaultValue=30, 
                                                     Description="The maximal dimensions to use for reduction"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        batchCorrection=ezFrame(Type="logical", 
                                                                DefaultValue="TRUE",
                                                                Description="Perform batch correction"),
                                        SCT.regress=ezFrame(Type="character", 
                                                            DefaultValue="none", 
                                                            Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."),
                                        DE.method=ezFrame(Type="charVector", 
                                                          DefaultValue="wilcox", 
                                                          Description="Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."),
                                        DE.regress=ezFrame(Type="charVector", 
                                                          DefaultValue="Plate", 
                                                          Description="Variables to regress out if the test LR is chosen"),
                                        chosenClusters=ezFrame(Type="charList",
                                                               DefaultValue="",
                                                               Description="The clusters to choose from each sample.In the format of sample1=cluster1,cluster2;sample2=cluster1,cluster2."),
                                        all2allMarkers=ezFrame(Type="logical",
                                                               DefaultValue=FALSE, 
                                                               Description="Run all against all cluster comparisons?"),
                                        markersToShow=ezFrame(Type="numeric", 
                                                              DefaultValue=10, 
                                                              Description="The markers to show in the heatmap of cluster marker genes"),
                                        maxSamplesSupported=ezFrame(Type="numeric", 
                                                                    DefaultValue=5, 
                                                                    Description="Maximum number of samples to compare"))
                }
              )
  )

ezMethodSCMultipleSamplesOneGroup = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(rlist)
  library(HDF5Array)
  

  if(input$getLength() > param$maxSamplesSupported){
    stop(paste("It only works for", param$maxSamplesSupported, "at most"))
  }
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  #the individual sce objects can be in hdf5 format (for new reports) or in rds format (for old reports)
  sceURLs <- input$getColumn("Static Report")
  filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(sceURLs)), "sce_h5")
  #In case it is in hdf5 format the Seurat object is not stored in the metadata slot, so we have to build it and store it there, since the 
  #clustering functions below work with sce objects and take the seurat object from them.
  if(file.exists(filePath)) {
    sceList <- lapply(filePath,loadHDF5SummarizedExperiment)
    names(sceList) <- names(sceURLs)
    #sceList <- lapply(sceList, function(sce) {metadata(sce)$scData <- CreateSeuratObject(counts=counts(sce),meta.data=data.frame(colData(sce)[,c(2:25, which(colnames(colData(sceList[[1]]))%in% "Condition"))])) 
    sceList <- lapply(sceList, function(sce) {metadata(sce)$scData <- CreateSeuratObject(counts=counts(sce),meta.data=data.frame(colData(sce))) 
    sce})
    #if it is an rds object it has been likely generated from old reports, so we need to update the seurat version before using the clustering functions below.                                             )
  } else {
    filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(sceURLs)), "sce.rds")
    sceList <- lapply(filePath,readRDS)
    names(sceList) <- names(sceURLs)
    sceList = lapply(sceList, update_seuratObjectVersion)
    sceList = lapply(sceList,  add_Condition_oldReports)
  }
  
  
  pvalue_allMarkers <- 0.05
  pvalue_all2allMarkers <- 0.01
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
  
  scData_noCorrected = cellClustNoCorrection(sceList, param)
  scData = scData_noCorrected
  if (param$batchCorrection) {
     scData_corrected = cellClustWithCorrection(sceList, param)
     #in order to compute the markers we switch again to the original assay
     DefaultAssay(scData_corrected) <- "SCT"
     scData = scData_corrected
  }
  scData@reductions$tsne_noCorrected <- Reductions(scData_noCorrected, "tsne")
  scData@meta.data$ident_noCorrected <- Idents(scData_noCorrected) 
    
  #positive cluster markers
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  
  #if all2allmarkers are not calculated it will remain as NULL
  all2allMarkers <- NULL
  
  #perform all pairwise comparisons to obtain markers
  if(doEnrichr(param) && param$all2allMarkers) 
    all2allMarkers <- all2all(scData, pvalue_all2allMarkers, param)
  
  #we do cell type identification using AUCell and SingleR
  cells_AUC <- NULL
  singler.results <- NULL
  #cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if(species == "Human" | species == "Mouse") {
    cells_AUC = cellsLabelsWithAUC(scData, species, tissue)
    singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "counts"), Idents(scData), species)
  }
  
  #Convert scData to Single Cell experiment Object
  # Use the SCT logcounts for visualization instead of the RNA logcounts.
  # TODO: save all the assays (RNA, SCT and integrated) in the sce object using the package keshavmot2/scanalysis. The function from Seurat doesn't save everything.
  DefaultAssay(scData) <- "SCT" 
  sce <- as.SingleCellExperiment(scData)
  metadata(sce)$cells_AUC <- cells_AUC
  metadata(sce)$singler.results <- singler.results
  metadata(sce)$output <- output
  metadata(sce)$param <- param
  metadata(sce)$param$name <- paste(param$name, paste(input$getNames(), collapse=", "), sep=": ")
  
  #Save some results in external files 
  saveExternalFiles(sce, list(pos_markers=posMarkers, all2allMarkers=all2allMarkers))
  # rowData(sce) = rowData(sce)[, c("gene_id", "biotypes", "description")]
  
  saveHDF5SummarizedExperiment(sce, dir="sce_h5")
  
  
  # Copy the style files and templates
   styleFiles <- file.path(system.file("templates", package="ezRun"), c("fgcz.css", "SCMultipleSamplesOneGroup.Rmd", "fgcz_header.html", "banner.png"))
 # styleFiles <- paste0("/home/daymegr/workspaceR/dayme-scripts/sushi_scripts_mod/", c("fgcz.css", "SCMultipleSamplesOneGroup.Rmd", "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCMultipleSamplesOneGroup.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, clean = TRUE, quiet=TRUE)
  rm(sceList)
  rm(scData)
  
  return("Success")
  
}

