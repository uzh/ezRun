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
                    perc_riboprot = ezFrame(
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
                    estimateAmbient = ezFrame(
                      Type = "logical", DefaultValue = TRUE,
                      Description = "estimate contamination with ambient RNA"
                    ),
                    controlSeqs = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "control sequences to add"
                    ),
                    CellRangerMulti = ezFrame(
                        Type = "logical", DefaultValue = FALSE,
                        Description = "unfiltered data of CellRanger multi output isn't sample specific and needs subsetting"
                    ),
                    enrichrDatabase = ezFrame(
                      Type = "charVector", DefaultValue = "", Description="enrichR databases to search"
                    ),
                    geneCountModel = ezFrame(
                      Type = "character",
                      DefaultValue = "GeneFull_ExonOverIntron",
                      Description = "(STARsolo Input Only) The gene count model, i.e. Solo features, to use from the previous step"
                    )
                  )
                }
              )
  )

ezMethodScSeurat <- function(input = NA, output = NA, param = NA,
                               htmlFile = "00index.html") {
  #library(scanalysis)
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
  cts <- Read10X(cmDir, gene.column = 1)
  featInfo <- ezRead.table(paste0(cmDir, "/features.tsv.gz"), header = FALSE, row.names = NULL)#, col_names = FALSE)
  colnames(featInfo) <- c("gene_id", "gene_name", "type")
  featInfo$isMito = grepl( "(?i)^MT-", featInfo$gene_name)
  featInfo$isRiboprot = grepl(  "(?i)^RPS|^RPL", featInfo$gene_name)
  geneAnnoFile <- sub("byTranscript", "byGene", param$ezRef@refAnnotationFile)
  if (file.exists(geneAnnoFile)){
    geneAnno <- ezRead.table(geneAnnoFile)
    if (any(geneAnno$type == "rRNA")){
      featInfo$isRibosomal <- geneAnno[featInfo$gene_id, "type"] == "rRNA"
    }
  }
  
  
  ## if we have feature barcodes we keep only the expression matrix
  if (is.list(cts)){
    cts <- cts$`Gene Expression`
    featInfo <- featInfo[  featInfo$type == "Gene Expression", ]
  }
  ## underscores in genenames will become dashes
  rownames(cts) <- rownames(featInfo) <- gsub("_", "-", uniquifyFeatureNames(ID=featInfo$gene_id, names=featInfo$gene_name)) 
  scData <- CreateSeuratObject(counts = cts[rowSums2(cts >0) >0, ])
  scData$Condition <- input$getColumn("Condition")
  scData@meta.data$Sample <- input$getNames()
  scData[["RNA"]] <- AddMetaData(object = scData[["RNA"]], metadata = featInfo[rownames(scData), ])
  scData$cellBarcode <- sub(".*_", "", colnames(scData))
  
  
  scData <- addCellQcToSeurat(scData, param=param, BPPARAM = BPPARAM, ribosomalGenes = featInfo[rownames(scData), "isRibosomal"])
  
  ## use empty drops to test for ambient
  rawDir <- sub("filtered_", "raw_", cmDir)
  if (file.exists(rawDir) && rawDir != cmDir){
    rawCts <- Read10X(rawDir, gene.column = 1)
    if(param$CellRangerMulti){
      rawCts <- rawCts[featInfo$gene_id,]
    }
    stopifnot(rownames(rawCts) == featInfo$gene_id)
    emptyStats <- emptyDrops(rawCts[!featInfo$isMito & !featInfo$isRiboprot, ],
                             BPPARAM=BPPARAM, niters=1e5)
    scData$negLog10CellPValue <- - log10(emptyStats[colnames(scData), "PValue"])
    emptyStats <- emptyDrops(rawCts, BPPARAM=BPPARAM, niters=1e5)
    scData$negLog10CellPValue <- pmin(scData$negLog10CellPValue, -log10(emptyStats[colnames(scData), "PValue"]))
    scData@meta.data$negLog10CellPValue[is.na(scData$negLog10CellPValue)] <- 0
    remove(rawCts)
  }
  allCellsMeta <- scData@meta.data
  
  #allCellsMeta$useCell <- scData$doubletClass %in% "singlet" & !scData$qc.lib & !scData$qc.mito & !scData$qc.nexprs & !scData$qc.riboprot
  scData <- subset(scData, cells=rownames(allCellsMeta)[allCellsMeta$useCell]) # %>% head(n=1000))
  
  ## remove low expressed genes
  num.cells <- param$cellsFraction * ncol(scData) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  is.expressed <- Matrix::rowSums(GetAssayData(scData, "counts") >= param$nUMIs) >= num.cells
  scData <- scData[is.expressed,]

  scData <- addCellCycleToSeurat(scData, param$refBuild, BPPARAM)
  
  vars.to.regress <- NULL
  if(identical("CellCycle", param$SCT.regress)){
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")
  } else {
    vars.to.regress <- NULL
  }
  ## generate normalized slots for the RNA assay
  scData <- NormalizeData(scData, normalization.method = "LogNormalize", scale.factor=10000, verbose=FALSE)
  scData <- FindVariableFeatures(scData, selection.method = "vst", verbose = FALSE, nfeatures=3000)
  scData <- ScaleData(scData, vars.to.regress = vars.to.regress, verbose=FALSE, do.scale=FALSE)
  ## generate the SCT assay
  scData <- SCTransform(scData, vst.flavor=2, vars.to.regress = vars.to.regress, seed.use = 38, verbose = FALSE,
                        return.only.var.genes=FALSE)
  #defaultAssay <- "SCT"
  ## defaultAssay is now SCT
  #scData <- FindVariableFeatures(scData, selection.method = "vst", verbose = FALSE)
  scData <- RunPCA(object=scData, npcs = param$npcs, verbose=FALSE)
  scData <- RunTSNE(object = scData, reduction = "pca", dims = 1:param$npcs)
  scData <- RunUMAP(object=scData, reduction = "pca", dims = 1:param$npcs)
  scData <- FindNeighbors(object = scData, reduction = "pca", dims = 1:param$npcs, verbose=FALSE)
  scData <- FindClusters(object=scData, resolution = seq(from = 0.2, to = 1, by = 0.2), verbose=FALSE)  #calculate clusters for a set of resolutions
  scData$seurat_clusters <- scData@meta.data[,paste0(DefaultAssay(scData), "_snn_res.", param$resolution)]  #but keep as the current clusters the ones obtained with the resolution set by the user
  Idents(scData) <- scData$seurat_clusters

  # positive cluster markers
  ## https://github.com/satijalab/seurat/issues/5321
  ## https://github.com/satijalab/seurat/issues/1501
  markers <- FindAllMarkers(object=scData, test.use = param$DE.method, only.pos=TRUE)
  ## Significant markers
  markers <- markers[ ,c("gene","cluster","pct.1", "pct.2", "avg_log2FC","p_val_adj")]
  #cm <- cm[cm$p_val_adj < 0.05, ]
  markers$cluster <- as.factor(markers$cluster)
  markers$diff_pct = abs(markers$pct.1-markers$pct.2)
  markers <- markers[order(markers$diff_pct, decreasing = TRUE),] ## why would we round here?? %>% mutate_if(is.numeric, round, digits=3)
  writexl::write_xlsx(markers, path="posMarkers.xlsx")
  
  if (ezIsSpecified(param$estimateAmbient) && param$estimateAmbient){
    scData <- addAmbientEstimateToSeurat(scData, rawDir=rawDir, threads = param$cores)
  }
  
  # cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if (species == "Human" | species == "Mouse") {
    genesPerCluster <- split(markers$gene, markers$cluster)
    enrichRout <- querySignificantClusterAnnotationEnrichR(genesPerCluster, param$enrichrDatabase)
    cells.AUC <- cellsLabelsWithAUC(GetAssayData(scData, "counts"), species, param$tissue, BPPARAM = BPPARAM)
    singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "data"), Idents(scData), species, BPPARAM = BPPARAM)
    for (r in names(singler.results)) {
      scData[[paste0(r,"_single")]] <- singler.results[[r]]$single.fine$labels
      scData[[paste0(r,"_cluster")]] <- singler.results[[r]]$cluster.fine$labels[match(Idents(scData), rownames(singler.results[[r]]$cluster.fine))]
    }
    saveRDS(cells.AUC, file="cells.AUC.rds")
    saveRDS(singler.results, file="singler.results.rds")
  } else {
    cells.AUC <- NULL
    singler.results <- NULL
    enrichRout <- NULL
  }
  
  #geneMeans <- geneMeansCluster(scData)
  
  ## generate template for manual cluster annotation -----
  ## we only deal with one sample
  stopifnot(length(input$getNames()) == 1)
  clusterInfos <- ezFrame(Sample=input$getNames(), Cluster=levels(Idents(scData)), ClusterLabel="")
  if (!is.null(singler.results)){
    clusterInfos$SinglerCellType <- singler.results$singler.results.cluster[clusterInfos$Cluster, "pruned.labels"]
  }
  nTopMarkers <- 10
  topMarkers <- markers %>% group_by(cluster) %>%
    slice_max(n = nTopMarkers, order_by = avg_log2FC)
  topMarkerString <- sapply(split(topMarkers$gene, topMarkers$cluster), paste, collapse=", ")
  clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
  clusterInfoFile <- "clusterInfos.xlsx"
  writexl::write_xlsx(clusterInfos, path=clusterInfoFile)
  
  makeRmdReport(param=param, output=output, scData=scData, allCellsMeta=allCellsMeta, enrichRout=enrichRout,
                cells.AUC=cells.AUC, singler.results=singler.results, rmdFile = "ScSeurat.Rmd", reportTitle = paste0(param$name, ": ",  input$getNames()))
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
  doubletsInfo <- scDblFinder(GetAssayData(scData, slot="counts")[ , scData$useCell], returnType = "table", clusters=TRUE, BPPARAM = BPPARAM)
  scData$doubletScore <- doubletsInfo[colnames(scData), "score"]
  scData$doubletClass <- doubletsInfo[colnames(scData), "class"]
  scData$qc.doublet <- scData$doubletClass %in% "doublet"
  scData$useCell <- scData$useCell & scData$doubletClass %in% "singlet"
  return(scData)
}

querySignificantClusterAnnotationEnrichR <- function(genesPerCluster, dbs, overlapGeneCutOff = 3, adjPvalueCutOff = 0.001, reportTopN = 5) {
  enrichRout <- list()
  for (cluster in unique(names(genesPerCluster))) {
    enriched <- enrichr(as.character(genesPerCluster[[cluster]]), dbs)
    
    for (db in names(enriched)) {
      enriched_db <- enriched[[db]]
      if (nrow(enriched_db) > 0 && colnames(enriched_db)[1] == "Term"){
        enriched_db$OverlapGenesN <- as.numeric(sapply(enriched_db$Overlap, function(x) str_split(x, "/")[[1]][1]))
        enriched_db$Cluster <- cluster
        enriched_db <- enriched_db %>%
          filter(., Adjusted.P.value < adjPvalueCutOff) %>%
          filter(., OverlapGenesN > overlapGeneCutOff) %>%
          head(reportTopN)
        enrichRout[[cluster]][[db]] <- enriched_db[, c("Term", "Cluster", "OverlapGenesN", "Adjusted.P.value", "Combined.Score")]
      }
    }
  }
  return(enrichRout)
}


