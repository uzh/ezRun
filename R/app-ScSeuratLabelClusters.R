###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratLabelClusters <-
  setRefClass("EzAppScSeuratLabelClusters",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScSeuratLabelClusters
                  name <<- "EzAppScSeuratLabelClusters"
                  appDefaults <<- rbind(
                    DE.method = ezFrame(
                      Type = "charVector",
                      DefaultValue = "wilcoxon",
                      Description = "Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."
                    ),
                    controlSeqs = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "control sequences to add"
                    ),
                    enrichrDatabase = ezFrame(
                      Type = "charVector", DefaultValue = "", Description="enrichR databases to search"
                    ),
                    computePathwayTFActivity=ezFrame(Type="logical", 
                                                     DefaultValue="TRUE",
                                                     Description="Whether we should compute pathway and TF activities.")
                  )
                }
              )
  )

ezMethodScSeuratLabelClusters <- function(input = NA, output = NA, param = NA,
                                          htmlFile = "00index.html") {
  library(HDF5Array)
  library(stringr)
  library(AUCell)
  library(GSEABase)
  library(SingleR)
  library(Seurat)
  library(scDblFinder)
  library(BiocParallel)
  library(scuttle)
  library(DropletUtils)
  library(enrichR)
  library(decoupleR)
  library(Azimuth)
  library(tidyverse)
  
  if (param$cores > 1){
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    ## scDblFinder fails with many cells and MulticoreParam
    BPPARAM <- SerialParam() 
  }
  register(BPPARAM)
  require(future)
  plan("multicore", workers = param$cores)
  set.seed(38)
  future.seed = TRUE
  options(future.rng.onMisuse="ignore")
  options(future.globals.maxSize = param$ram*1024^3)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("SC Cluster Report")))
  on.exit(setwd(cwd), add = TRUE)
  
  # load cluster annotation file
  clusterAnnoFn <- file.path(param$dataRoot, param$ClusterAnnotationFile)
  if (ezIsSpecified(param$ClusterAnnotationFile)) {
    clusterAnnoFn <- file.path(param$dataRoot, param$ClusterAnnotationFile)
    stopifnot("The cluster annotation file does not exist or is not an .xlsx file!" = 
                file.exists(clusterAnnoFn) && str_ends(clusterAnnoFn, ".xlsx$"))
  } else {
    stop("Must supply cluster annotation file path.", call.=FALSE)
  }
  clusterAnno <- readxl::read_xlsx(clusterAnnoFn) %>% 
    as_tibble() %>%
    dplyr::select(1:3) %>% # remove all other columns
    dplyr::rename(c("Sample"=1, "Cluster"=2, "ClusterLabel"=3)) %>%
    dplyr::filter(Sample == input$getNames())
  labelMap <- as.character(clusterAnno$ClusterLabel)
  names(labelMap) <- as.character(clusterAnno$Cluster)
  
  # load cell data
  allCellsMeta <- readRDS(file.path(input$getFullPaths("SC Cluster Report"), "allCellsMeta.rds"))
  scData <- readRDS(input$getFullPaths("SC Seurat"))
  
  # change labels and store in a variable
  scData$cellType <- unname(labelMap[as.character(scData$seurat_clusters)])
  Idents(scData) <- scData$cellType
  
  # get markers and annotations
  anno <- getSeuratMarkersAndAnnotate(scData, param)
  
  # save markers
  markers <- anno$markers
  writexl::write_xlsx(markers, path="posMarkers.xlsx")
  
  ## generate template for manual cluster annotation -----
  ## we only deal with one sample
  stopifnot(length(input$getNames()) == 1)
  clusterInfos <- ezFrame(Sample=input$getNames(), Cluster=levels(Idents(scData)), ClusterLabel="")
  if (!is.null(anno$singler.results)){
    clusterInfos$SinglerCellType <- anno$singler.results$singler.results.cluster[clusterInfos$Cluster, "pruned.labels"]
  }
  nTopMarkers <- 10
  topMarkers <- markers %>% group_by(cluster) %>%
    slice_max(n = nTopMarkers, order_by = avg_log2FC)
  topMarkerString <- sapply(split(topMarkers$gene, topMarkers$cluster), paste, collapse=", ")
  clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
  clusterInfoFile <- "clusterInfos.xlsx"
  writexl::write_xlsx(clusterInfos, path=clusterInfoFile)
  
  makeRmdReport(param=param, output=output, scData=scData, allCellsMeta=allCellsMeta, 
                enrichRout=anno$enrichRout, cells.AUC=anno$cells.AUC, 
                singler.results=anno$singler.results, aziResults=anno$aziResults,
                pathwayActivity=anno$pathwayActivity, TFActivity=anno$TFActivity,
                rmdFile = "ScSeurat.Rmd", reportTitle = paste0(param$name, ": ",  input$getNames()))
  #remove no longer used objects
  rm(scData)
  gc()
  return("Success")
}
