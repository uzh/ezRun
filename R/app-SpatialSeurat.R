###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSpatialSeurat <-
  setRefClass("EzAppSpatialSeurat",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSpatialSeurat
                  name <<- "EzAppspatialSeurat"
                  appDefaults <<- rbind(npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=20,
                                                    Description="The maximal dimensions to use for reduction"),
                                        pcGenes=ezFrame(Type="charVector", 
                                                        DefaultValue="", 
                                                        Description="The genes used in supvervised clustering"),
                                        SCT.regress=ezFrame(Type="character", 
                                                           DefaultValue="none", 
                                                           Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."),
                                        DE.method=ezFrame(Type="charVector", 
                                                          DefaultValue="wilcoxon", 
                                                          Description="Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.5,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        cellsFraction=ezFrame(Type="numeric", 
                                                                DefaultValue=0.05, 
                                                                Description="A gene will be kept if it is expressed in at least this percentage of cells"),
                                        nUMIs=ezFrame(Type="numeric", 
                                                      DefaultValue=1, 
                                                      Description='A gene will be kept if it has at least nUMIs in the fraction of cells specified before'),
                                        nmad=ezFrame(Type="numeric", 
                                                     DefaultValue=3, 
                                                     Description="Median absolute deviation (MAD) from the median value of each metric across all cells")
                                        )
                }
              )
  )

ezMethodSpatialSeurat <- function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  library(Seurat)
  scData <- load10xSpatialData(input, param)
  scData_list <- filterCellsAndGenes(scData, param) # return sce objects filtered and unfiltered to show the QC metrics later in the rmd
  scData <- scData_list$scData
  scData.unfiltered <- scData_list$scData.unfiltered
  rm(scData_list)
  scData <- seuratClusteringV3(scData, param, assay="Spatial")
  
  pvalue_allMarkers <- 0.05
  
  #positive cluster markers
  clusterMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  
  #spatially variable genes
  scData <- spatialMarkers(scData)
  spatialMarkers <- SpatiallyVariableFeatures(scData, selection.method = "markvariogram")
  
  #Save some results in external files 
  dataFiles = saveExternalFiles(list(cluster_markers=clusterMarkers, spatial_markers=data.frame(spatialMarkers)))
  saveRDS(scData, "scData.rds")
  saveRDS(scData.unfiltered, "scData.unfiltered.rds")
  saveRDS(param, "param.rds")
  
  makeRmdReport(dataFiles=dataFiles, rmdFile = "spatialSeurat.Rmd", reportTitle = param$name) 
  return("Success")
}

filterCellsAndGenes <- function(scData, param) {
  library(scater)
  # Cells filtering
  scData <- PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
  scData <- PercentageFeatureSet(scData, "(?i)^RPS|^RPL", col.name = "percent_ribo")
  
  if (param$nreads == "") {
    qc.lib <- isOutlier(scData$nCount_Spatial, log = TRUE, nmads = param$nmad, type = "lower")
  } else {
    qc.lib <- scData$nCount_Spatial < as.double(param$nreads)
  }
  if (param$ngenes == "") {
    qc.nexprs <- isOutlier(scData$nFeature_Spatial, nmads = param$nmad, log = TRUE, type = "lower")
  } else {
    qc.nexprs <- scData$nFeature_Spatial < as.double(param$ngenes)
  }
  if (param$perc_mito == "") {
    qc.mito <- isOutlier(scData$percent_mito, nmads = param$nmad, type = "higher")
  } else {
    qc.mito <- scData$percent_mito > as.double(param$perc_mito)
  }
  if (param$perc_ribo == "") {
    qc.ribo <- isOutlier(scData$percent_ribo, nmads = param$nmad, type = "higher")
  } else {
    qc.ribo <- scData$percent_ribo > as.double(param$perc_ribo)
  }
  
  discard <- qc.lib | qc.nexprs | qc.mito | qc.ribo
  scData$discard <- discard
  scData$qc.lib <- qc.lib
  scData$qc.nexprs <- qc.nexprs
  scData$qc.mito <- qc.mito
  scData$qc.ribo <- qc.ribo
  scData.unfiltered <- scData
  scData <- scData[, !discard]
  
  # Genes filtering
  num.cells <- param$cellsFraction * ncol(scData) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  is.expressed <- Matrix::rowSums(GetAssayData(scData, "counts") >= param$nUMIs) >= num.cells
  scData[["Spatial"]] <- AddMetaData(object = scData[["Spatial"]], metadata = is.expressed,col.name ='is.expressed')
  return(list(scData.unfiltered = scData.unfiltered, scData = scData))
}

