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
                    npcs = ezFrame(
                      Type = "numeric",
                      DefaultValue = 20,
                      Description = "The maximal dimensions to use for reduction"
                    ),
                    pcGenes = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "The genes used in supvervised clustering"
                    ),
                    SCT.regress = ezFrame(
                      Type = "character",
                      DefaultValue = "none",
                      Description = "Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
                    ),
                    DE.method = ezFrame(
                      Type = "charVector",
                      DefaultValue = "wilcoxon",
                      Description = "Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."
                    ),
                    resolution = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.5,
                      Description = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
                    ),
                    nreads = ezFrame(
                      Type = "numeric",
                      DefaultValue = Inf,
                      Description = "Low quality cells have less than \"nreads\" reads. Only when applying fixed thresholds."
                    ),
                    ngenes = ezFrame(
                      Type = "numeric",
                      DefaultValue = Inf,
                      Description = "Low quality cells have less than \"ngenes\" genes. Only when applying fixed thresholds."
                    ),
                    perc_mito = ezFrame(
                      Type = "numeric",
                      DefaultValue = Inf,
                      Description = "Low quality cells have more than \"perc_mito\" percent of mitochondrial genes. Only when applying fixed thresholds."
                    ),
                    perc_ribo = ezFrame(
                      Type = "numeric",
                      DefaultValue = Inf,
                      Description = "Low quality cells have more than \"perc_ribo\" percent of ribosomal genes. Only when applying fixed thresholds."
                    ),
                    cellsFraction = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.01,
                      Description = "A gene will be kept if it is expressed in at least this percentage of cells"
                    ),
                    nUMIs = ezFrame(
                      Type = "numeric",
                      DefaultValue = 1,
                      Description = "A gene will be kept if it has at least nUMIs in the fraction of cells specified before"
                    ),
                    nmad = ezFrame(
                      Type = "numeric",
                      DefaultValue = 3,
                      Description = "Median absolute deviation (MAD) from the median value of each metric across all cells"
                    ),
                    filterByExpression = ezFrame(
                      Type = "character", DefaultValue = FALSE,
                      Description = "Keep cells according to specific gene expression. i.e. Set > 1 | Pkn3 > 1"
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
  library(scanalysis)
  library(HDF5Array)
  library(AUCell)
  library(GSEABase)
  library(SingleR)
  library(Seurat)
  library(SingleCellExperiment)
  library(tidyverse)
  library(scanalysis)
  require(scDblFinder)
  
  #clusterInfo <- readRDS("clusterInfoFile.rds")
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  pvalue_allMarkers <- 0.05
  sce <- loadHDF5SummarizedExperiment(input$getFullFilename("SC H5"))
  counts(sce) <- as(counts(sce), "sparseMatrix")
  logcounts(sce) <- as(logcounts(sce), "sparseMatrix")
  clusterInfo <- readxl::read_xlsx(param$clusterInfo)
  clusterInfo <- clusterInfo[clusterInfo$Sample %in% input$getNames(), ] 
  clusterInfo <- clusterInfo[!is.na(clusterInfo$ClusterLabel) & clusterInfo$ClusterLabel != "", ] 
  clusterMap <- setNames(clusterInfo$ClusterLabel, clusterInfo$Cluster)
  sce$cellType <- sce$ident %>% recode(!!!clusterMap)
  

  pvalue_allMarkers <- 0.05
  pvalue_all2allMarkers <- 0.01
  

  # positive cell type markers
  scData <- as.Seurat(sce)
  Idents(scData) <- scData$cellType
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  
  cells_AUC <- NULL
  singler.results <- NULL
  # cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if (species == "Human" | species == "Mouse") {
    cells_AUC <- cellsLabelsWithAUC(scData, species, param$tissue)
    singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "counts"), Idents(scData), species)
  }
  
  sce <- scData %>% seurat_to_sce(default_assay = "SCT")
  metadata(sce)$PCA_stdev <- Reductions(scData_diet, "pca")@stdev
  metadata(sce)$cells_AUC <- cells_AUC
  metadata(sce)$singler.results <- singler.results
  metadata(sce)$output <- output
  metadata(sce)$param <- param
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
                                    paste(input$getNames(), collapse = ", "),
                                    sep = ": "
  )
  geneMeans <- geneMeansCluster(sce)
  

  # Save some results in external files
  dataFiles <- saveExternalFiles(list(pos_markers = posMarkers, gene_means = as_tibble(as.data.frame(geneMeans), rownames = "gene_name")))
  # rowData(sce) = rowData(sce)[, c("gene_id", "biotypes", "description")]
  
  saveHDF5SummarizedExperiment(sce, dir = "sce_h5")

  makeRmdReport(dataFiles = dataFiles, clusterInfoFile=clusterInfoFile, rmdFile = "SCLabelClusters.Rmd", reportTitle = metadata(sce)$param$name)
  return("Success")
}
