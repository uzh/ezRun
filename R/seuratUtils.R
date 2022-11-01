###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

seuratPreProcess <- function(sce){
  ## parameters to tune: param$minCellsPerGene;
  ##                     param$maxGenesPerCell; param$minGenesPerCell; param$maxMitoFraction
  ##                     param$minReadsPerCell; 
  ##                     param$pcs, param$pcGenes
  require(Seurat)
  require(scater)
  param <- metadata(sce)$param
  rownames(sce) <- uniquifyFeatureNames(ID=rowData(sce)$gene_id,
                                        names=rowData(sce)$gene_name)
  if(toupper(param$scProtocol) == "SMART-SEQ2"){
    sce <- sce[ ,Matrix::colSums(assays(sce)$counts) > param$minReadsPerCell]
  }
  cell_info <- data.frame(colData(sce),
                          Plate=sub("___.*$", "", colnames(sce)),
                          check.names = FALSE)
  scData <- CreateSeuratObject(raw.data=assays(sce)$counts,
                               min.cells=param$minCellsPerGene,
                               min.genes=1,
                               project=param$name,
                               meta.data=cell_info)
  mito.genes <- rownames(sce)[toupper(as.character(seqnames(rowRanges(sce)))) %in% 
                                toupper(c("chrM", "MT"))]
  mito.genes <- intersect(mito.genes, rownames(scData@data))
  
  perc_mito <- Matrix::colSums(scData@raw.data[mito.genes, ])/Matrix::colSums(scData@raw.data)
  scData <- AddMetaData(object = scData, metadata = perc_mito,
                        col.name = "perc_mito")
  
  scData <- FilterCells(object = scData,
                        subset.names = c("nGene", "perc_mito"),
                        low.thresholds = c(param$minGenesPerCell, -Inf), 
                        high.thresholds = c(param$maxGenesPerCell, param$maxMitoFraction))
  scData <- NormalizeData(object=scData, normalization.method="LogNormalize",
                          scale.factor=getSeuratScalingFactor(param$scProtocol))
  scData <- seuratClustering(scData, param)
  
  metadata(sce)$scData <- scData
  
  return(sce)
}

seuratClustering <- function(scData, param){
  set.seed(38)
  scData <- FindVariableGenes(object = scData, do.plot = FALSE,
                              x.low.cutoff=param$x.low.cutoff,
                              x.high.cutoff=param$x.high.cutoff,
                              y.cutoff=param$y.cutoff)
  scData <- ScaleData(object = scData, do.par=TRUE,
                      vars.to.regress = param$vars.to.regress,
                      num.cores=param$cores)
  if(ezIsSpecified(param$pcGenes)){
    indicesMatch <- match(toupper(param$pcGenes),
                          toupper(rownames(scData@data)))
    if(any(is.na(indicesMatch))){
      stop("The following genes don't exist: ", 
           paste(param$pcGenes[is.na(indicesMatch)], collapse = ","))
    }
    pc.genes <- rownames(scData@data)[indicesMatch]
  }else{
    pc.genes <- scData@var.genes
  }
  scData <- RunPCA(object=scData, pc.genes=pc.genes, pcs.compute=20,
                   do.print=FALSE, pcs.print=1:5,
                   genes.print=5)
  scData <- ProjectPCA(object = scData, do.print = FALSE)
  scData <- JackStraw(object=scData, num.replicate=100) #, display.progress=FALSE,
                      #do.par=TRUE, num.cores=min(4L, param$cores))
  
  scData <- FindClusters(object=scData, reduction.type="pca",
                         dims.use = 1:min(param$pcs, length(pc.genes)-1),
                         resolution = param$resolution, print.output = 0, 
                         save.SNN=TRUE, force.recalc=FALSE)
  scData <- RunTSNE(object=scData, reduction.use = "pca",
                    dims.use=1:min(param$pcs, length(pc.genes)-1), tsne.method="Rtsne",
                    perplexity=ifelse(length(scData@ident) > 200, 30, 10),
                    num_threads=param$cores)
  try(scData <- RunUMAP(object=scData, reduction.use = "pca",
                        dims.use=1:min(param$pcs, length(pc.genes)-1),
                        n_neighbors=ifelse(length(scData@ident) > 200, 30, 10)),
      silent=TRUE)
  return(scData)
}

getSeuratScalingFactor <- function(x){
  x <- toupper(x)
  ans <- switch(x,
                "SMART-SEQ2"=1e5,
                "10X"=1e4)
  return(ans)
}

cellTable <- function(scData){
  toTable <- tibble(Cluster=names(summary(scData@ident)),
                    "# of cells"=summary(scData@ident))
  cellCountsByPlate <- tibble(Plate=scData@meta.data$Plate,
                              Cluster=as.character(scData@ident)) %>%
    group_by(Plate, Cluster) %>% summarise(n()) %>%
    spread(Plate, `n()`, fill=0)
  cellPercByPlate <- select(cellCountsByPlate, -Cluster) %>%
    as.matrix()
  rownames(cellPercByPlate) <- cellCountsByPlate$Cluster
  cellPercByPlate <- sweep(cellPercByPlate, 2, colSums(cellPercByPlate), "/")
  colnames(cellPercByPlate) <- paste0(colnames(cellPercByPlate), "_fraction")
  toTable <- left_join(toTable, cellCountsByPlate, by="Cluster") %>%
    left_join(as_tibble(cellPercByPlate, rownames="Cluster"), by="Cluster")
  ## TODO: add fisher test?
  toTable <- bind_rows(toTable,
                       bind_cols("Cluster"="Total", 
                                 summarise_at(toTable, setdiff(colnames(toTable), "Cluster"),
                                              sum)))
  return(toTable)
}


update_seuratObjectVersion = function(se) {
  if(se@metadata$scData@version < 3)
    metadata(se)$scData = UpdateSeuratObject(metadata(se)$scData)
  return(se)
}

add_Condition_oldReports <- function(sce) {
  scData = metadata(sce)$scData
  scData@meta.data$Condition = sce[,colnames(scData)]$Condition
  metadata(sce)$scData = scData
  sce
}

seuratStandardWorkflow <- function(scData, param){
  require(future)
  plan("multicore", workers = param$cores)
  
  set.seed(38)
  if (identical(param$pcGenes, character(0))) 
     features <- NULL
  else {
     species <- getSpecies(param$refBuild)
     if(species == 'Human'){
         features <- stringr::str_to_upper(param$pcGenes)
     } else {
        features <- stringr::str_to_title(param$pcGenes)
     }
  }
  future.seed = TRUE
  options(future.rng.onMisuse="ignore")
  options(future.globals.maxSize = param$ram*1024^3)
  
  scData <- RunPCA(object=scData, npcs = param$npcs, features=features, verbose=FALSE)
  scData <- RunTSNE(object = scData, reduction = "pca", dims = 1:param$npcs)
  scData <- RunUMAP(object=scData, reduction = "pca", dims = 1:param$npcs)
  scData <- FindNeighbors(object = scData, reduction = "pca", dims = 1:param$npcs, verbose=FALSE)
  scData <- FindClusters(object=scData, resolution = seq(from = 0.1, to = 1, by = 0.1), verbose=FALSE)  #calculate clusters for a set of resolutions
  Idents(scData) <- scData@meta.data[,paste0(DefaultAssay(scData), "_snn_res.", param$resolution)]  #but keep as the current clusters the ones obtained with the resolution set by the user
  scData@meta.data$ident <- Idents(scData)
  return(scData)
}  

buildSeuratObject <- function(sce){
  library(Seurat)
  library(scater)
  param <- metadata(sce)$param
  colData(sce)[, grep("SCT", colnames(colData(sce)))] = NULL  #remove all normalization info done on each object separately
  cell_info <- data.frame(colData(sce),
                          Plate=sub("___.*$", "", colnames(sce)),
                          check.names = FALSE)
  scData <- CreateSeuratObject(counts=assays(sce)$counts,
                               project=param$name,
                               meta.data=cell_info)
  
  #scData[["RNA"]]@meta.features <- cbind.data.frame(scData[["RNA"]]@meta.features, data.frame(rowData(sce)[, c("gene_id", "biotypes", "description")]))
  return(scData)
}

seuratClusteringV3 <- function(scData, param, assay="RNA") {
  vars.to.regress <- NULL
  if(identical("CellCycle", param$SCT.regress))
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")

  scData <- SCTransform(scData, vars.to.regress = vars.to.regress, assay=assay, seed.use = 38, verbose = FALSE)
  scData <- seuratStandardWorkflow(scData, param)
  return(scData)
}

seuratClusteringHTO <- function(scData) {
  DefaultAssay(scData) <- "HTO"
  scData <- ScaleData(scData)
  scData <- RunPCA(scData, features = rownames(scData), reduction.name = "pca_hto", reduction.key = "pca_hto_", 
                           verbose = FALSE)
  # Now, we rerun tSNE using the PCA only on ADT (protein) levels.
  scData <- RunTSNE(scData, reduction = "pca_hto", reduction.key = "htoTSNE_", reduction.name = "tsne_hto", check_duplicates = FALSE)
  scData <- FindNeighbors(scData, reduction="pca_hto", features = rownames(scData), dims=NULL)
  scData <- FindClusters(scData, resolution = 0.2)  
  return(scData)
}

cellClustNoCorrection <- function(scDataList, param) {
   if(param[['name']] == 'SpatialSeuratSlides')
    assay = "Spatial"
  else 
    assay= "RNA"
  vars.to.regress <- NULL
  if(identical("CellCycle", param$SCT.regress))
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")
  #1. Data preprocesing is done on each object separately
  for (i in 1:length(scDataList)) 
    scDataList[[i]] <- SCTransform(scDataList[[i]], vars.to.regress = vars.to.regress,  assay = assay, verbose = TRUE)
  #2. Merge all seurat objects
  scData = Reduce(function(x, y) {merge(x,y, merge.data = TRUE)}, scDataList)
  VariableFeatures(scData) <- unlist(lapply(scDataList, VariableFeatures))
  #3. Dimensionality reduction and clustering
  scData <- seuratStandardWorkflow(scData, param)
  
  return(scData)
}

cellClustWithCorrection <- function (scDataList, param) {
  if(param[['name']] == 'SpatialSeuratSlides')
    assay = "Spatial"
  else 
    assay= "RNA"
  vars.to.regress <- NULL
  if(identical("CellCycle", param$SCT.regress))
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")
  #1. Data preprocesing is done on each object separately
  for (i in 1:length(scDataList)) 
    scDataList[[i]] <- SCTransform(scDataList[[i]], vars.to.regress = vars.to.regress,  assay = assay, verbose = TRUE)
  
  #2. Data integration
  #2.1. # Select the most variable features to use for integration
  integ_features <- SelectIntegrationFeatures(object.list = scDataList, nfeatures = 3000) 
  #2.2. Prepare the SCT list object for integration
  scDataList <- PrepSCTIntegration(object.list = scDataList, anchor.features = integ_features)
  if(param$integrationMethod == 'RPCA'){
    scDataList <- lapply(X = scDataList, FUN = RunPCA, features = integ_features)
  }
  #2.3. Find anchors
  if(param$integrationMethod == 'Classic'){
  integ_anchors <- FindIntegrationAnchors(object.list = scDataList, normalization.method = "SCT", 
                                          anchor.features = integ_features, dims = 1:param$npcs)
  } else if(param$integrationMethod == 'RPCA'){
  integ_anchors <- FindIntegrationAnchors(object.list = scDataList, normalization.method = "SCT", 
                                          anchor.features = integ_features, dims = 1:param$npcs, reduction = "rpca", k.anchor = 20)
  }
  
  #2.4. Integrate datasets
  seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", dims = 1:param$npcs)

  #3. Run the standard workflow for visualization and clustering
  seurat_integrated <- seuratStandardWorkflow(seurat_integrated, param)
  
  return(seurat_integrated)
}

posClusterMarkers <- function(scData, pvalue_allMarkers, param) {
  vars.to.regress = NULL
  if(param$name == "SCOneSample" & param$DE.method == "LR") #For the one sample app, we will regress out cell cycle if the test is LR
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")
  else if(param$name == "SCReportMerging" & param$DE.method == "LR") #for multiple samples we will regress out either the cell cycle, plate or both if the test is LR
    vars.to.regress <- param$DE.regress
  
  markers <- FindAllMarkers(object=scData, test.use = param$DE.method, only.pos=TRUE, latent.vars = vars.to.regress, return.thresh = pvalue_allMarkers)
  ## Significant markers
  cm <- markers[ ,c("gene","cluster","pct.1", "pct.2", "avg_log2FC","p_val_adj")]
  #cm <- cm[cm$p_val_adj < 0.05, ]
  cm$cluster <- as.factor(cm$cluster)
  diff_pct = abs(cm$pct.1-cm$pct.2)
  cm$diff_pct <- diff_pct
  cm <- cm[order(cm$diff_pct, decreasing = TRUE),] %>% mutate_if(is.numeric, round, digits=3)
  rownames(cm) <- NULL
  return(cm)
}

spatialMarkers <- function(scData) { 
  scData <- FindSpatiallyVariableFeatures(scData, features = VariableFeatures(scData)[1:1000], 
                                                  selection.method = "markvariogram")
  spatialMarkers <- SpatiallyVariableFeatures(scData, selection.method = "markvariogram")
  return(spatialMarkers)
  }

all2all <- function(scData, pvalue_all2allMarkers, param) {
  clusterCombs <- combn(levels(Idents(scData)), m=2)
  all2allMarkers <- mcmapply(FindMarkers, as.integer(clusterCombs[1, ]), as.integer(clusterCombs[2, ]),
                             MoreArgs = list(object=scData,only.pos=FALSE),
                             mc.preschedule=FALSE,
                             mc.cores=min(4L, param$cores),
                             SIMPLIFY=FALSE)
  all2allMarkers <- lapply(all2allMarkers, function(x){
    x[x$p_val <= pvalue_all2allMarkers, ]
  })
  names(all2allMarkers) <- apply(clusterCombs, 2, paste, collapse="vs")
  return(all2allMarkers)
}

conservedMarkers <- function(scData, grouping.var="Condition") {
  markers <- list()
  
  if("SCT" %in% Seurat::Assays(scData)) {
    assay <- "SCT"
  } else {
    assay <- "RNA"
  }
  for(eachCluster in levels(Idents(scData))){
    markersEach <- try(FindConservedMarkers(scData, ident.1=eachCluster, 
                                            grouping.var=grouping.var, 
                                            print.bar=FALSE, only.pos=TRUE, 
                                            assay = assay), silent=TRUE)
    if(class(markersEach) != "try-error" && nrow(markersEach) > 0){
      markers[[eachCluster]] <- as_tibble(markersEach, rownames="gene")
    }
  }
  markers <- bind_rows(markers, .id="cluster")
  markers <- markers %>%
    mutate(avg_avg_fc=rowMeans(select(., contains("_avg_log2FC"))))
  return(markers)
}

diffExpressedGenes <- function(scData, param) {
  seurat_clusters <- Idents(scData)
  scData@meta.data$cluster.condition <- paste0(seurat_clusters, "_", scData@meta.data$Condition)
  Idents(scData) <- "cluster.condition"
  conditions <- unique(scData@meta.data$Condition)
  
  vars.to.regress = NULL
  if(param$DE.method == "LR") #regress the plate if the test is LR
    vars.to.regress <- param$DE.regress
  
    diffGenes <- list()
    for(eachCluster in gtools::mixedsort(levels(seurat_clusters))){
      markersEach <- try(FindMarkers(scData, ident.1=paste0(eachCluster, "_",
                                                            param$sampleGroup),
                                     ident.2=paste0(eachCluster, "_", 
                                                    param$refGroup),
                                     test.use = param$DE.method, latent.vars = vars.to.regress))
      ## to skip some groups with few cells
      if(class(markersEach) != "try-error"){
        diffGenes[[eachCluster]] <- as_tibble(markersEach, rownames="gene")
      }
    }
    diffGenes <- bind_rows(diffGenes, .id="cluster")
  
  diff_pct = abs(diffGenes$pct.1-diffGenes$pct.2)
  diffGenes$diff_pct <- diff_pct
  diffGenes <- diffGenes[order(diffGenes$diff_pct, decreasing = TRUE),]
  rownames(diffGenes) <- NULL
 # diffGenes <- diffGenes[diffGenes$p_val_adj < 0.05, ]
  
  return(diffGenes)
}

