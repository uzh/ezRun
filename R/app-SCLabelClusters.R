###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCLabelClusters <-
  setRefClass("EzAppSCLabelClusters",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCLabelClusters
                  name <<- "EzAppSCLabelClusters"
                  appDefaults <<- rbind(
                    DE.method = ezFrame(
                      Type = "charVector",
                      DefaultValue = "wilcox",
                      Description = "Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."
                    ),
                    resolution = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.5,
                      Description = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
                    ),
                    controlSeqs = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "control sequences to add"
                    )
                  )
                }
              )
  )

ezMethodSCLabelClusters <- function(input = NA, output = NA, param = NA,
                                htmlFile = "00index.html") {
  library(HDF5Array)
  library(AUCell)
  library(GSEABase)
  library(SingleR)
  library(Seurat)
  library(SingleCellExperiment)
  library(tidyverse)
  require(scDblFinder)
  
  #clusterInfo <- readRDS("clusterInfoFile.rds")
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("SC Celltype Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  sce <- loadHDF5SummarizedExperiment(input$getFullPaths("SC H5"))
  counts(sce) <- as(counts(sce), "sparseMatrix")
  logcounts(sce) <- as(logcounts(sce), "sparseMatrix")
  clusterInfo <- readxl::read_xlsx(file.path(param$dataRoot, param$clusterInfo))
  clusterInfo <- clusterInfo[clusterInfo$Sample %in% input$getNames(), ] 
  clusterInfo <- clusterInfo[!is.na(clusterInfo$ClusterLabel) & clusterInfo$ClusterLabel != "", ] 
  clusterMap <- setNames(clusterInfo$ClusterLabel, clusterInfo$Cluster)
  if (length(clusterMap) > 0){
    sce$cellType <- sce$ident %>% recode(!!!clusterMap)
  } else {
    sce$cellType <- sce$ident
  }
  

  pvalue_allMarkers <- 0.05
  pvalue_all2allMarkers <- 0.01
  

  # positive cell type markers
  scData <- as.Seurat(sce)
  Idents(scData) <- scData$cellType
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  

  assayName <- head(intersect(c("originalexp", "RNA"), Assays(scData)), 1)
  sce <- scData %>% seurat_to_sce(default_assay = assayName)
  #metadata(sce)$PCA_stdev <- Reductions(scData, "PCA")@stdev
  metadata(sce)$output <- output
  metadata(sce)$param <- param
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
                                    paste(input$getNames(), collapse = ", "),
                                    sep = ": "
  )
  geneMeans <- geneMeansCluster(sce)

  # species <- getSpecies(param$refBuild)
  # if (species == "Human" | species == "Mouse") {
  #   cells_AUC <- cellsLabelsWithAUC(scData, species, param$tissue)
  #   singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "counts"), Idents(scData), species)
  # }
  # metadata(sce)$singler.results <- singler.results
  # metadata(sce)$cells_AUC <- cells_AUC
  

  # Save some results in external files
  dataFiles <- saveExternalFiles(list(pos_markers = posMarkers, gene_means = as_tibble(as.data.frame(geneMeans), rownames = "gene_name")))
  # rowData(sce) = rowData(sce)[, c("gene_id", "biotypes", "description")]
  
  saveHDF5SummarizedExperiment(sce, dir = "sce_h5")

  makeRmdReport(dataFiles = dataFiles, rmdFile = "SCLabelClusters.Rmd", reportTitle = metadata(sce)$param$name)
  return("Success")
}
