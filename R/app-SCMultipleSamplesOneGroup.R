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
                                                                Description="Which batch correction method to use? None or CCA"),
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
                                                                    Description="Maximum number of samples to compare"), 
                                        species=ezFrame(Type="character", DefaultValue="Human", Description="Organism"))
                }
              )
  )

ezMethodSCMultipleSamplesOneGroup = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library("Seurat")
  library(rlist)
  library(tibble)
  library(readr)

  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  if(input$getLength() > param$maxSamplesSupported){
    stop(paste("It only works for", param$maxSamplesSupported, "at most"))
  }
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  sceURLs <- input$getColumn("Static Report")
  sceList <- lapply(file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(sceURLs)), "sce.rds"),readRDS)
  names(sceList) <- names(sceURLs)
  sceList = lapply(sceList, update_seuratObjectVersion)
  
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
     DefaultAssay(scData_corrected) <- "RNA"
     scData_corrected <- ScaleData(scData_corrected, verbose = FALSE)
     scData = scData_corrected
  }
     
  #positive cluster markers
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers)
  
  #perform all pairwise comparisons to obtain markers
  if(doEnrichr(param) && param$all2allMarkers) 
    all2allMarkers <- all2all(scData, pvalue_all2allMarkers, param)
  
  if(param$species == "Human" | param$species == "Mouse") {
    cells_AUC = cellsLabelsWithAUC(scData, param)
    singler.results <- cellsLabelsWithSingleR(scData, param)
  }
  
  if(param$batchCorrection) {
    scData@reductions$tsne_noCorrected <- Reductions(scData_noCorrected, "tsne")
    scData@meta.data$ident_noCorrected <- Idents(scData_noCorrected)
  } 
  
  #Convert scData to Single Cell experiment Object
  sce <- as.SingleCellExperiment(scData)
  metadata(sce)$cells_AUC <- cells_AUC
  metadata(sce)$singler.results <- singler.results
  metadata(sce)$output <- output
  metadata(sce)$param <- param
  metadata(sce)$param$name <- paste(param$name, paste(input$getNames(), collapse=", "), sep=": ")
  
  #Save some results in external files 
  saveExternalFiles(sce, list(pos_markers=posMarkers, all2allMarkers=all2allMarkers))
  # rowData(sce) = rowData(sce)[, c("gene_id", "biotypes", "description")]
  
  library(HDF5Array)
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

