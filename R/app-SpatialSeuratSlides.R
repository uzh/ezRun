###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSpatialSeuratSlides <-
  setRefClass("EzAppSpatialSeuratSlides",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSpatialSeuratSlides
                  name <<- "EzAppSpatialSeuratSlides"
                  appDefaults <<- rbind(npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=30, 
                                                    Description="The maximal dimensions to use for reduction"),
                                        nfeatures = ezFrame(
                                            Type = "numeric",
                                            DefaultValue = 3000,
                                            Description = "number of variable genes for SCT"
                                        ),
                                        pcGenes = ezFrame(
                                          Type = "charVector",
                                          DefaultValue = "",
                                          Description = "The genes used in supvervised clustering"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        SCT.regress.CellCycle=ezFrame(
                                          Type = "logical", 
                                          DefaultValue = FALSE,
                                          Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
                                        ),
                                        batchCorrection=ezFrame(Type="logical", 
                                                                DefaultValue="TRUE",
                                                                Description="Perform batch correction."),
                                        integrationMethod=ezFrame(Type="character", 
                                                                  DefaultValue="CCA", 
                                                                  Description="Choose integration method in Seurat"),
                                        DE.method=ezFrame(Type="charVector", 
                                                          DefaultValue="wilcox", 
                                                          Description="Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"),
                                        DE.regress=ezFrame(Type="charVector", 
                                                           DefaultValue="Batch", 
                                                           Description="Variables to regress out if the test LR is chosen"),
                                        min.pct = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = 0.1,
                                          Description = "Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations."
                                        ),
                                        logfc.threshold = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = 0.25,
                                          Description = "Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
                                        ),
                                        pt.size.factor = ezFrame(
                                            Type = "numeric",
                                            DefaultValue = 1,
                                            Description = "pt.size.factor for spatial plots"
                                        ),
                                        maxSamplesSupported=ezFrame(Type="numeric", 
                                                              DefaultValue=5, 
                                                              Description="Maximum number of samples to compare"))
                }
              )
  )

ezMethodSpatialSeuratSlides = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(rlist)
  library(HDF5Array)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  scDataURLs <- input$getColumn("Static Report")
  filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(scDataURLs)), "scData.rds")
  
  scDataList <- lapply(filePath,readRDS)
  names(scDataList) <- names(scDataURLs)
  scDataList <- lapply(scDataList, function(scData) {
  scData@meta.data[, grep("SCT", colnames(scData@meta.data))] = NULL #remove previous clustering done on SCT assay
  scData})
  
  pvalue_allMarkers <- 0.05
 
  scData_noCorrected <- cellClustNoCorrection(scDataList, param)
  scData = scData_noCorrected
  if (param$batchCorrection) {
    scData_corrected = cellClustWithCorrection(scDataList, param)
    #in order to compute the markers we switch again to the original assay
    varFeatures <- VariableFeatures(scData_corrected)
    DefaultAssay(scData_corrected) <- "SCT"
    scData <- scData_corrected
    VariableFeatures(scData) <- unique(varFeatures)
  }
  scData@reductions$tsne_noCorrected <- Reductions(scData_noCorrected, "tsne")
  scData@reductions$umap_noCorrected <- Reductions(scData_noCorrected, "umap")
  scData@meta.data$ident_noCorrected <- Idents(scData_noCorrected)
  scData <- PrepSCTFindMarkers(scData)
  
  #positive cluster markers
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  posMarkers[['isSpatialMarker']] = FALSE
  #spatially variable genes
  require(readxl)
  filePath_spatialMarkers <- sub('scData.rds', 'spatialMarkers.xlsx',filePath)
  spatialMarkersList <- lapply(filePath_spatialMarkers,read_xlsx)
  names(spatialMarkersList) <- names(scDataList)
  spatialMarkers <- c()
  for (j in 1:length(spatialMarkersList)){
      spatialMarkersList[[j]][['GeneSymbol']] = spatialMarkersList[[j]]$GeneSymbol
      spatialMarkersList[[j]][['SampleID']] = names(spatialMarkersList)[j]
      spatialMarkersList[[j]] = spatialMarkersList[[j]][!is.na(spatialMarkersList[[j]]$MeanRank),]
      spatialMarkers <- rbind(spatialMarkers, spatialMarkersList[[j]])
  }
  
  spatialPosMarkers <- intersect(posMarkers$gene, unique(spatialMarkers$GeneSymbol))
  posMarkers[which(posMarkers$gene %in% spatialPosMarkers), 'isSpatialMarker'] = TRUE
 
  #Save some results in external files 
  dataFiles = saveExternalFiles(list(pos_markers=posMarkers, spatial_markers=spatialMarkers))
  saveRDS(scData, "scData.rds")
  saveRDS(param, "param.rds")
  
  makeRmdReport(dataFiles=dataFiles, output=output, rmdFile = "SpatialSeuratSlides.Rmd", reportTitle = param$name) 
  return("Success")
  
}






