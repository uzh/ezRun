###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratFilterClusters <-
  setRefClass("EzAppScSeuratFilterClusters",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScSeuratFilterClusters
                  name <<- "EzAppScSeuratFilterClusters"
                  appDefaults <<- rbind(
                    nfeatures = ezFrame(
                          Type = "numeric",
                          DefaultValue = 3000,
                          Description = "number of variable genes for SCT"
                      ),
                    npcs = ezFrame(
                      Type = "numeric",
                      DefaultValue = 20,
                      Description = "The maximal dimensions to use for reduction"
                    ),
                    SCT.regress.CellCycle = ezFrame(
                      Type = "logical", 
                      DefaultValue = FALSE,
                      Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
                    ),
                    DE.method = ezFrame(
                      Type = "charVector",
                      DefaultValue = "wilcoxon",
                      Description = "Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."
                    ),
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
                    resolution = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.6,
                      Description = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
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

ezMethodScSeuratFilterClusters <- function(input = NA, output = NA, param = NA,
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
    dplyr::filter(Sample == input$getNames()) %>%
    dplyr::mutate(ClusterLabel=str_replace(ClusterLabel, "(?i)remove", "REMOVE"),
                  ClusterLabel=str_replace(ClusterLabel, "(?i)keep", "KEEP"))
  stopifnot("No matching sample in `Sample` column of cluster annotation! Check the sample names match." = 
              nrow(clusterAnno) > 0)
  
  labelMap <- as.character(clusterAnno$ClusterLabel)
  names(labelMap) <- as.character(clusterAnno$Cluster)
  
  # load cell data
  allCellsMeta <- readRDS(file.path(input$getFullPaths("SC Cluster Report"), "allCellsMeta.rds"))
  scData <- readRDS(input$getFullPaths("SC Seurat"))
  
  # change labels and store in a variable
  toKeep <- unname(labelMap[as.character(Idents(scData))]) != "REMOVE"
  # We don't need to do anything if all clusters from the sample are kept.
  # Happens if other samples other than this one in the dataset are changed.
  if (!all(toKeep)) {
    scData <- scData[, toKeep]
    DefaultAssay(scData) <- "RNA"
    scData <- DietSeurat(scData, assays="RNA", layers="counts")
    
    ## Get information on which variables to regress out in scaling/SCT
    scData <- seuratStandardSCTPreprocessing(scData, param)
    
    ## defaultAssay is now SCT
    scData <- seuratStandardWorkflow(scData, param, ident.name="seurat_clusters")
  }

  if ("cellType" %in% colnames(scData@meta.data)) {
    # in case where the input dataset is labeled
    Idents(scData) <- scData$cellType
  }
  
  # get markers and annotations
  anno <- getSeuratMarkersAndAnnotate(scData, param, BPPARAM = BPPARAM)
  
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
