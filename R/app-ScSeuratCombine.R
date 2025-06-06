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
                  appDefaults <<- rbind(
                                        nfeatures = ezFrame(
                                        Type = "numeric",
                                        DefaultValue = 3000,
                                        Description = "number of variable genes for SCT"
                                        ),
                                        npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=30, 
                                                    Description="The maximal dimensions to use for reduction"),
                                        pcGenes = ezFrame(
                                          Type = "charVector",
                                          DefaultValue = "",
                                          Description = "The genes used in supvervised clustering"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        integrationMethod=ezFrame(Type="character", 
                                                                  DefaultValue="CCA", 
                                                                  Description="Choose integration method in Seurat (CCA or RPCA)"),
                                        enrichrDatabase=ezFrame(Type = "charVector", 
                                                                DefaultValue = "", 
                                                                Description="enrichR databases to search"),
                                        computePathwayTFActivity=ezFrame(Type="logical", 
                                                                 DefaultValue="TRUE",
                                                                 Description="Whether we should compute pathway and TF activities."),
                                        SCT.regress.CellCycle=ezFrame(
                                          Type = "logical", 
                                          DefaultValue = FALSE,
                                          Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
                                        ),
                                        DE.method=ezFrame(
                                          Type="charVector", 
                                          DefaultValue="wilcox", 
                                          Description="Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"),
                                        DE.regress=ezFrame(
                                          Type="charVector", 
                                          DefaultValue="Batch", 
                                          Description="Variables to regress out if the test LR is chosen"),
                                        min.pct = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = 0.1,
                                          Description = "Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations."
                                        ),
                                        min.diff.pct = ezFrame(
                                            Type = "numeric",
                                            DefaultValue = 0,
                                            Description = "Used in filtering cluster markers"
                                        ),
                                        logfc.threshold = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = 0.25,
                                          Description = "Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
                                        ))
                }
              )
  )

ezMethodScSeuratCombine = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(rlist)
  library(HDF5Array)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(AUCell)
  library(enrichR)
  library(decoupleR)
  library(Azimuth)
  library(BiocParallel)
  require(future)
  
  plan("multicore", workers = param$cores)
  set.seed(38)
  future.seed = TRUE
  options(future.rng.onMisuse="ignore")
  options(future.globals.maxSize = param$ram*1024^3)
    
  BPPARAM <- MulticoreParam(workers = param$cores)
  register(BPPARAM)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  message("Attempting to load Seurat data...")
  filePath <- file.path("/srv/gstore/projects", input$getColumn("SC Seurat"))
  filePath_course <- file.path("/srv/GT/analysis/course_sushi/public/projects", input$getColumn("SC Seurat"))
  if(!file.exists(filePath[1])){
      filePath <- filePath_course
  }
  names(filePath) <- input$getNames()
  
  if (length(filePath) < 2) {
    stop("need at least two samples to combine.")
  }
  
  # Load the data and prepare metadata for integration
  scDataList <- lapply(names(filePath), function(sm) {
    scData <- ezLoadRobj(filePath[sm], nthreads=param$cores)
    aziFilePath <- file.path(dirname(filePath[sm]),'aziResults.rds')
    if(file.exists(aziFilePath)){
        aziResults <- readRDS(aziFilePath)
        scData <- AddMetaData(scData, aziResults)
    }
    scData$Sample <- sm
    # If we have new information in the Condition column, add it to the dataset
    if (all(scData$Condition == "NA" | scData$Condition == "") || 
           (ezIsSpecified(param$overwriteCondition) && as.logical(param$overwriteCondition))) {
      if (any(startsWith(colnames(input$meta), "Condition"))) {
        scData$Condition <- unname(input$getColumn("Condition")[sm])
      } else {
        scData$Condition <- scData$Sample
      }
    }
    # Harmony will complain if the Condition is the same across all samples
    if (param$integrationMethod == "Harmony" && 
        length(unique(input$meta$`Condition`)) == 1) {
      scData$Condition <- scData$Sample
    } else if (ezIsSpecified(param$STACASAnnotationFile)) {
      clusterAnnoFn <- file.path(param$dataRoot, param$STACASAnnotationFile)
      clusterAnno <- readxl::read_xlsx(clusterAnnoFn) %>% 
        as_tibble() %>%
        dplyr::select(1:3) %>% # remove all other columns
        dplyr::rename(c("Sample"=1, "Cluster"=2, "ClusterLabel"=3)) %>%
        dplyr::filter(Sample == sm)
      labelMap <- as.character(clusterAnno$ClusterLabel)
      names(labelMap) <- as.character(clusterAnno$Cluster)
      scData$stacasLabelColumn <- unname(labelMap[as.character(Idents(scData))])
    }
    # Also add the other factors in the input dataset to the objects
    if (ezIsSpecified(param$additionalFactors)) {
      additionalFactors <- str_split(param$additionalFactors, ",", simplify=TRUE)[1,]
      metaFactorNames <- paste0("meta_", additionalFactors) %>% str_replace(., " ", ".")
      names(metaFactorNames) <- additionalFactors
      for (hf in additionalFactors) {
        scData[[metaFactorNames[hf]]] <- unname(input$getColumn(hf)[sm])
      }
    }
    # Rename the cells and add original sample-level clusters back in
    scData <- RenameCells(scData, new.names = paste0(scData$Sample, "-", colnames(scData)))
    if(class(scData$seurat_clusters) == 'character'){
        scData$sample_seurat_clusters <- paste0(scData$Sample, "-", sprintf("%s", scData$seurat_clusters))
    } else {
        scData$sample_seurat_clusters <- paste0(scData$Sample, "-", sprintf("%02d", scData$seurat_clusters))
    }
    return(scData)
  })
  
  # perform all of the analysis
  results <- seuratIntegrateDataAndAnnotate(scDataList, input, output, param, BPPARAM)
  
  # generate ClusterInfos table
  clusterInfos <- ezFrame(Samples=paste(input$getNames(),collapse=','), Cluster=levels(Idents(results$scData)), ClusterLabel="")
  if (!is.null(results$singler.results)){
      clusterInfos$SinglerCellType <- results$singler.results$singler.results.cluster[clusterInfos$Cluster, "pruned.labels"]
  }
  nTopMarkers <- 10
  topMarkers <- results$markers %>% group_by(cluster) %>%
      slice_max(n = nTopMarkers, order_by = avg_log2FC)
  topMarkerString <- sapply(split(topMarkers$gene, topMarkers$cluster), paste, collapse=", ")
  clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
  clusterInfoFile <- "clusterInfos.xlsx"
  writexl::write_xlsx(clusterInfos, path=clusterInfoFile)
  
  # save the markers
  writexl::write_xlsx(results$markers, path="posMarkers.xlsx")
  qs2::qs_save(results$scData, "scData.qs2", nthreads = param$cores)
  
  # Save some results in external files
  reportTitle <- 'SCReport - MultipleSamples based on Seurat'
  makeRmdReport(param=param, output=output, scData=results$scData, 
                enrichRout=results$enrichRout, TFActivity=results$TFActivity, 
                pathwayActivity=results$pathwayActivity, aziResults=results$aziResults,
                cells.AUC=results$cells.AUC, singler.results=results$singler.results,
                rmdFile = "ScSeuratCombine.Rmd", reportTitle = reportTitle) 
  return("Success")
  
}

seuratIntegrateDataAndAnnotate <- function(scDataList, input, output, param, BPPARAM = SerialParam()) {
  pvalue_allMarkers <- param$pvalue_allMarkers
  
  if(ezIsSpecified(param$chosenClusters)){
    for(eachSample in names(param$chosenClusters)){
      chosenCells <- names(Idents(scDataList[[eachSample]]))[Idents(scDataList[[eachSample]]) %in% param$chosenClusters[[eachSample]]]
      scDataList[[eachSample]] <- scDataList[[eachSample]][, chosenCells]
    }
  }
  
  scData_noCorrected <- cellClustNoCorrection(scDataList, param)
  if (param$integrationMethod!='none') {
    scData_corrected = cellClustWithCorrection(scDataList, param)
    #in order to compute the markers we switch again to the original assay
    DefaultAssay(scData_corrected) <- "SCT"
    scData <- scData_corrected
  } else {
    scData = scData_noCorrected
  }
  scData@reductions$tsne_noCorrected <- Reductions(scData_noCorrected, "tsne")
  Key(scData@reductions$tsne_noCorrected) <- 'TSNEnoCorrection_'
  scData@reductions$umap_noCorrected <- Reductions(scData_noCorrected, "umap")
  Key(scData@reductions$umap_noCorrected) <- 'UMAPnoCorrection_'

  # Since Seurat v5, the RNA layers are split by sample. Does not affect SCT assay
  # but clutters up RNA assay should it be used in downstream analysis
  scData <- JoinLayers(scData, assay="RNA")
  
  scData@meta.data$ident_noCorrected <- Idents(scData_noCorrected)
  scData <- PrepSCTFindMarkers(scData)
  
  # get annotation information
  anno <- getSeuratMarkersAndAnnotate(scData, param, BPPARAM = BPPARAM)
  
  return(list(scData=scData, 
              markers=anno$markers,
              enrichRout=anno$enrichRout, 
              pathwayActivity=anno$pathwayActivity, 
              TFActivity=anno$TFActivity,
              cells.AUC=anno$cells.AUC,
              singler.results=anno$singler.results,
              aziResults=anno$aziResults))
}
