###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeurat <-
  setRefClass("EzAppScSeurat",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScSeurat
                  name <<- "EzAppScSeurat"
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
                    pcGenes = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "The genes used in unsupervised clustering"
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
                    min.diff.pct = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0,
                      Description = "Used for filtering cluster markers: The minimum difference of cell fraction of the two tested populations."
                    ),
                    logfc.threshold = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.25,
                      Description = "Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
                    ),
                    pvalue_allMarkers = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.01,
                      Description = "Used for filtering cluster markers: adjusted pValue threshold for marker detection"
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
                    perc_riboprot = ezFrame(
                      Type = "numeric",
                      DefaultValue = Inf,
                      Description = "Low quality cells have more than \"perc_ribo\" percent of ribosomal genes. Only when applying fixed thresholds."
                    ),
                    keepDoublets = ezFrame(
                      Type = "logical",
                      DefaultValue = FALSE,
                      Description = "Whether we should keep cells suspected of being doublets. Set to TRUE only for QC purposes."
                    ),
                    maxEmptyDropPValue = ezFrame(
                      Type = "numeric",
                      DefaultValue = 1,
                      Description = "filter droplets based on DropletUtils::emptyDrops method"
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
                    estimateAmbient = ezFrame(
                      Type = "logical", DefaultValue = TRUE,
                      Description = "estimate contamination with ambient RNA"
                    ),
                    controlSeqs = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "control sequences to add"
                    ),
                    enrichrDatabase = ezFrame(
                      Type = "charVector", DefaultValue = "", Description="enrichR databases to search"
                    ),
                    geneCountModel = ezFrame(
                      Type = "character",
                      DefaultValue = "GeneFull_ExonOverIntron",
                      Description = "(STARsolo Input Only) The gene count model, i.e. Solo features, to use from the previous step"
                    ),
                    computePathwayTFActivity=ezFrame(Type="logical", 
                                                     DefaultValue="TRUE",
                                                     Description="Whether we should compute pathway and TF activities."),
                    excludeGenes = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "file path to txt file with gene symbols to exclude from the analysis")
                  )
                }
              )
  )

ezMethodScSeurat <- function(input = NA, output = NA, param = NA,
                               htmlFile = "00index.html") {
  library(HDF5Array)
  library(AUCell)
  library(GSEABase)
  library(SingleR)
  library(Seurat)
  library(tidyverse)
  library(scDblFinder)
  library(BiocParallel)
  library(scuttle)
  library(DropletUtils)
  library(enrichR)
  library(decoupleR)
  library(Azimuth)
  
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
  
  cmDir <- input$getFullPaths("CountMatrix")
  if (file.exists(file.path(cmDir, param$geneCountModel))){
    cmDir <- file.path(cmDir, param$geneCountModel)
  }
  if(grepl('h5$',cmDir[1])){
      param$cellbender = TRUE
  } else {
      param$cellbender = FALSE 
  }
  if(!param$cellbender){
    cts <- Read10X(cmDir, gene.column = 1)
    featInfo <- ezRead.table(paste0(cmDir, "/features.tsv.gz"), header = FALSE, row.names = NULL)
  } else if(param$cellbender){
    cts <- Read10X_h5(file.path(dirname(cmDir), 'cellbender_filtered_seurat.h5'), use.names = FALSE)
    inputDS <- ezRead.table(file.path(dirname(dirname(cmDir)),'input_dataset.tsv'))
    cellrangerDir <- file.path(param$dataRoot,inputDS[input$getNames(),'CountMatrix [Link]'])
    param[['cellrangerDir']] <- cellrangerDir
    featInfo <- ezRead.table(paste0(cellrangerDir, "/features.tsv.gz"), header = FALSE, row.names = NULL)
  }
  featInfo <- featInfo[,1:3]  # in cases where additional column exist, e.g. CellRangerARC output
  colnames(featInfo) <- c("gene_id", "gene_name", "type")
  featInfo$isMito = grepl( "(?i)^MT-", featInfo$gene_name)
  featInfo$isRiboprot = grepl(  "(?i)^RPS|^RPL", featInfo$gene_name)
  geneAnnoFile <- sub("byTranscript", "byGene", param$ezRef@refAnnotationFile)
  if (file.exists(geneAnnoFile)){
    geneAnno <- ezRead.table(geneAnnoFile)
    if (any(geneAnno$type == "rRNA")){
      featInfo$isRibosomal <- geneAnno[featInfo$gene_id, "type"] == "rRNA"
      if(any(is.na(featInfo[, "isRibosomal"]))){
          featInfo[, "isRibosomal"][which(is.na(featInfo[, "isRibosomal"]))] <- FALSE
      }
    }
  }
  
  ## if we have feature barcodes we keep only the expression matrix
  if (is.list(cts)){
    cts <- cts$`Gene Expression`
    featInfo <- featInfo[  featInfo$type == "Gene Expression", ]
  }
  if(param$cellbender){
    rownames(featInfo) <- featInfo$gene_id
    matchingIds <- intersect(rownames(cts), rownames(featInfo))
    cts <- cts[matchingIds,]
    featInfo <- featInfo[matchingIds,]
  }
  
  ## underscores in genenames will become dashes
  rownames(cts) <- rownames(featInfo) <- gsub("_", "-", uniquifyFeatureNames(ID=featInfo$gene_id, names=featInfo$gene_name)) 
  scData <- CreateSeuratObject(counts = cts[rowSums2(cts >0) >0, ])
  scData$Condition <- unname(input$getColumn("Condition"))
  scData@meta.data$Sample <- input$getNames()
  scData[["RNA"]] <- AddMetaData(object = scData[["RNA"]], metadata = featInfo[rownames(scData), ])
  scData$cellBarcode <- sub(".*_", "", colnames(scData))
  scData <- addCellQcToSeurat(scData, param=param, BPPARAM = BPPARAM, ribosomalGenes = featInfo[rownames(scData), "isRibosomal"])
  
  ## use empty drops to test for ambient       
  if ("UnfilteredCountMatrix" %in% input$colNames) {
    rawDir <- input$getFullPaths("UnfilteredCountMatrix")
    if (file.exists(file.path(rawDir, param$geneCountModel))) {
      rawDir <- file.path(rawDir, param$geneCountModel)
    }
  } else {
    # DEPRECATED; all input datasets should specify the path to the unfiltered count matrix
    rawDir <- sub("filtered_", "raw_", cmDir)
  }
  if (file.exists(rawDir) && rawDir != cmDir){
    if(param$cellbender){
        rawCts <- Read10X_h5(file.path(dirname(cmDir), 'cellbender_raw_seurat.h5'), use.names = FALSE)
    } else {
        rawCts <- Read10X(rawDir, gene.column = 1)
    }
    if (is.list(rawCts)) {
      rawCts <- rawCts$`Gene Expression`
      rawCts <- rawCts[featInfo$gene_id,]
    }
      
    if (("SCDataOrigin" %in% input$colNames) && 
        input$getColumn("SCDataOrigin") == 'BDRhapsody') {
      rawCts <- rawCts[featInfo$gene_id,]
    }
    
    if(param$cellbender){
      rawCts <- rawCts[featInfo$gene_id,]
    }
    
     if(length(setdiff(rownames(rawCts), featInfo$gene_id)) > 0){
      rawCts <- rawCts[featInfo$gene_id,] 
     }
    stopifnot(rownames(rawCts) == featInfo$gene_id)
    emptyStats <- emptyDrops(rawCts[!featInfo$isMito & !featInfo$isRiboprot, ],
                             BPPARAM=BPPARAM, niters=1e5)
    scData$negLog10CellPValue <- - log10(emptyStats[colnames(scData), "PValue"])
    emptyStats <- emptyDrops(rawCts, BPPARAM=BPPARAM, niters=1e5)
    scData$negLog10CellPValue <- pmin(scData$negLog10CellPValue, -log10(emptyStats[colnames(scData), "PValue"]))
    scData@meta.data$negLog10CellPValue[is.na(scData$negLog10CellPValue)] <- 0
    scData$qc.empty <- FALSE
    
    if(param$maxEmptyDropPValue < 1){
        scData$qc.empty[scData$negLog10CellPValue < -log10(param$maxEmptyDropPValue)] <- TRUE
        scData$useCell[scData$qc.empty] <- FALSE
    }
    remove(rawCts)
  }
  allCellsMeta <- scData@meta.data
  
  scData <- subset(scData, cells=rownames(allCellsMeta)[allCellsMeta$useCell]) # %>% head(n=1000))
  
  ## remove lowly expressed genes
  num.cells <- param$cellsFraction * ncol(scData) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  cellsPerGene <- Matrix::rowSums(GetAssayData(scData, layer="counts") >= param$nUMIs)
  is.expressed <- cellsPerGene >= num.cells
  cellsPerGeneFraction <- data.frame(frac = cellsPerGene/ncol(scData), row.names = rownames(cellsPerGene))
  scData <- scData[is.expressed,]
  
  if(ezIsSpecified(param$excludeGenes) && param$excludeGenes!=''){
      genesToExclude <- ezRead.table(param$excludeGenes, header = FALSE, row.names = NULL)
      genesToExclude <- unique(genesToExclude$V1)
      genesToKeep <- setdiff(rownames(scData), genesToExclude)
      scData <- subset(scData, features = genesToKeep)
  }
  
  ## Add Cell Cycle information to Seurat object as metadata columns
  scData <- addCellCycleToSeurat(scData, param$refBuild, BPPARAM)
  scData$CC.Difference <- scData$CellCycleS - scData$CellCycleG2M
  
  ## Get information on which variables to regress out in scaling/SCT
  scData <- seuratStandardSCTPreprocessing(scData, param)
  ## defaultAssay is now SCT
  scData <- seuratStandardWorkflow(scData, param, ident.name="seurat_clusters")

  # estimate ambient first
  if (ezIsSpecified(param$estimateAmbient) && param$estimateAmbient) {
    scData <- addAmbientEstimateToSeurat(scData, rawDir=rawDir, param = param)
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
  if (!is.null(anno$aziResults)){
    for (nm in grep("celltype", colnames(anno$aziResults), value=TRUE)){
      cellCounts <- table(cluster=scData$seurat_clusters, sample=anno$aziResults[[nm]])
      cellPerc <- sweep(cellCounts, 1, rowSums(cellCounts), "/")
      percMat <- as.matrix(cellPerc)
      newLabels <- apply(percMat, 1, function(x){colnames(percMat)[x >0.5]}) %>% unlist()
      clusterInfos[[nm]] <- clusterInfos$Cluster %>% as.character() %>% recode(!!!newLabels)
    }
  }
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
                cellsPerGeneFraction = cellsPerGeneFraction, enrichRout=anno$enrichRout, 
                cells.AUC=anno$cells.AUC, singler.results=anno$singler.results, aziResults=anno$aziResults,
                pathwayActivity=anno$pathwayActivity, TFActivity=anno$TFActivity, cellxgeneResults=anno$cellxgeneResults,
                rmdFile = "ScSeurat.Rmd", reportTitle = paste0(param$name, ": ",  input$getNames()))
  #remove no longer used objects
  rm(scData)
  gc()
  return("Success")
}

addCellQcToSeurat <- function(scData, param=NULL, BPPARAM=NULL, ribosomalGenes=NULL){
  
  library(scater)
  
  # Cells filtering
  scData <- PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
  scData <- PercentageFeatureSet(scData, "(?i)^RPS|^RPL", col.name = "percent_riboprot")
  if (!is.null(ribosomalGenes)){
    scData <- PercentageFeatureSet(scData, features=ribosomalGenes, col.name = "percent_ribosomal")
  }
  if(grepl("Spatial", param$appName)) {
    assay <- "Spatial"
    att_nCounts <- "nCount_Spatial"
    att_nGenes <- "nFeature_Spatial"
  } else {
    att_nCounts <- "nCount_RNA"
    att_nGenes <- "nFeature_RNA"
    assay <- "RNA"
  }
  
  if (!ezIsSpecified(param$nreads)) {
    scData$qc.lib <- isOutlier(scData@meta.data[,att_nCounts], log = TRUE, nmads = param$nmad, type = "lower")
  } else {
    scData$qc.lib <- scData@meta.data[,att_nCounts] < param$nreads
  }
  if (!ezIsSpecified(param$ngenes)) {
    scData$qc.nexprs <- isOutlier(scData@meta.data[,att_nGenes], nmads = param$nmad, log = TRUE, type = "lower")
  } else {
    scData$qc.nexprs <- scData@meta.data[,att_nGenes] < param$ngenes
  }
  if (!ezIsSpecified(param$perc_mito)) {
    scData$qc.mito <- isOutlier(scData@meta.data[,"percent_mito"], nmads = param$nmad, type = "higher")
  } else {
    scData$qc.mito <- scData@meta.data[,"percent_mito"] > param$perc_mito
  }
  if (!ezIsSpecified(param$perc_riboprot )) {
    scData$qc.riboprot <- isOutlier(scData@meta.data[,"percent_riboprot"], nmads = param$nmad, type = "higher")
  } else {
    scData$qc.riboprot <- scData@meta.data[,"percent_riboprot"] > as.numeric(param$perc_riboprot)
  }
  
  scData$useCell <- !(scData$qc.lib | scData$qc.nexprs | scData$qc.mito | scData$qc.riboprot)
  
  set.seed(38)
  doubletsInfo <- scDblFinder(GetAssayData(scData, layer="counts")[ , scData$useCell], returnType = "table", clusters=TRUE, BPPARAM = BPPARAM)
  scData$doubletScore <- doubletsInfo[colnames(scData), "score"]
  scData$doubletClass <- doubletsInfo[colnames(scData), "class"]
  scData$qc.doublet <- scData$doubletClass %in% "doublet"
  if (ezIsSpecified(param$keepDoublets) && param$keepDoublets) {
    futile.logger::flog.info("Keeping doublets...")
  } else {
    scData$useCell <- scData$useCell & scData$doubletClass %in% "singlet"
  }
  return(scData)
}

querySignificantClusterAnnotationEnrichR <- function(genesPerCluster, dbs, overlapGeneCutOff = 3, adjPvalueCutOff = 0.001, reportTopN = 5) {
  enrichRout <- list()
  for (cluster in unique(names(genesPerCluster))) {
    enriched <- enrichr(as.character(genesPerCluster[[cluster]]), dbs)
    
    for (db in names(enriched)) {
      enriched_db <- enriched[[db]]
      if (nrow(enriched_db) > 0 && colnames(enriched_db)[1] == "Term"){
        enriched_db$OverlapGenesN <- sub("/.*", "", enriched_db$Overlap) %>% as.numeric()
        enriched_db$Cluster <- cluster
        enriched_db <- enriched_db %>%
          filter(., Adjusted.P.value < adjPvalueCutOff) %>%
          filter(., OverlapGenesN > overlapGeneCutOff) %>%
          head(reportTopN)
        enrichRout[[cluster]][[db]] <- enriched_db[, c("Term", "Cluster", "Overlap", "OverlapGenesN", "Adjusted.P.value", "Odds.Ratio", "Combined.Score")]
      }
    }
  }
  return(enrichRout)
}


computeTFActivityAnalysis <- function(cells, species){
  species <- tolower(species)
  # Retrieve prior knowledge network.
  network <- decoupleR::get_dorothea(organism = species,
                                     levels = c("A", "B", "C"))
  
  # Run weighted means algorithm.
  activities <- decoupleR::run_wmean(mat = as.matrix(GetAssayData(cells)),
                                     network = network,
                                     .source = "source",
                                     .targe = "target",
                                     .mor = "mor",
                                     times = 100,
                                     minsize = 5)
  
  return(activities)
}


computePathwayActivityAnalysis <- function(cells, species){
  species <- tolower(species)
  # Retrieve prior knowledge network.
  network <- decoupleR::get_progeny(organism = species)
  
  # Run weighted means algorithm.
  activities <- decoupleR::run_wmean(mat = as.matrix(GetAssayData(cells)),
                                     network = network,
                                     .source = "source",
                                     .targe = "target",
                                     .mor = "weight",
                                     times = 100,
                                     minsize = 5)
  
  return(activities)
}
