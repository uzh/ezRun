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
  scData <- JackStraw(object=scData, num.replicate=100, display.progress=FALSE,
                      do.par=TRUE, num.cores=min(4L, param$cores))
  
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
  require(tidyverse)
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

seuratStandardWorkflow <- function(scData, param){
  set.seed(38)
  if (identical(param$pcGenes, character(0))) 
    param$pcGenes <- NULL
  
  scData <- RunPCA(object=scData, npcs = param$npcs, features=param$pcGenes)
  scData <- RunTSNE(object = scData, reduction = "pca", dims = 1:param$npcs, num_threads=param$cores)
  scData <- RunUMAP(object=scData, reduction = "pca", dims = 1:param$npcs, n_threads=param$cores)
  scData <- FindNeighbors(object = scData, reduction = "pca", dims = 1:param$npcs)
  scData <- FindClusters(object=scData, resolution = seq(from = 0.1, to = 1, by = 0.1))  #calculate clusters for a set of resolutions
  Idents(scData) <- scData@meta.data[,paste0(DefaultAssay(scData), "_snn_res.", param$resolution)]  #but keep as the current clusters the ones obtained with the resolution set by the user
  scData@meta.data$seurat_clusters <- Idents(scData)
  return(scData)
}  
  
cellsProportion <- function(scData){
  require(tidyverse)
  toTable <- tibble(Cluster=names(summary(Idents(scData))),
                    "# of cells"=summary(Idents(scData)))
  cellCountsByPlate <- tibble(Plate=scData@meta.data$Plate,
                              Cluster=as.character(Idents(scData))) %>%
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

buildSeuratObject <- function(sce){
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
  scData <- CreateSeuratObject(counts=assays(sce)$counts,
                               project=param$name,
                               meta.data=cell_info)
  
  #scData[["RNA"]]@meta.features <- cbind.data.frame(scData[["RNA"]]@meta.features, data.frame(rowData(sce)[, c("gene_id", "biotypes", "description")]))
  return(scData)
}

seuratClusteringV3 <- function(scData, param) {
  vars.to.regress <- NULL
  if(identical("CellCycle", param$SCT.regress))
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")
  scData <- SCTransform(scData, vars.to.regress = vars.to.regress, seed.use = 38, verbose = TRUE)
  scData <- seuratStandardWorkflow(scData, param)
  return(scData)
}

cellClustNoCorrection <- function(sceList, param) {
  #Merge all seurat objects
  scData = Reduce(merge, lapply(sceList, function(x){metadata(x)$scData}))
  scData@project.name <- param$name
  # when doing the scaling, normalization and feature slection with SCTransform we will only regress out by cell cycle if specified
  vars.to.regress <- NULL
  if(identical("CellCycle", param$SCT.regress))
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")
  scData <- SCTransform(scData, vars.to.regress = vars.to.regress, seed.use = 38, verbose = TRUE)
  scData <- seuratStandardWorkflow(scData, param)
  return(scData)
}

cellClustWithCorrection <- function (sceList, param) {
  seurat_objects = lapply(sceList, function(se) {metadata(se)$scData})
  
  # when doing the scaling, normalization and feature slection with SCTransform we will only regress out by cell cycle if specified
  vars.to.regress <- NULL
  if(identical("CellCycle", param$SCT.regress))
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")
  #1. Data preprocesing
  for (i in 1:length(seurat_objects)) {
    seurat_objects[[i]] <- SCTransform(seurat_objects[[i]], vars.to.regress = vars.to.regress, verbose = TRUE)
  }
  #2. Data integration
  #2.1. # Select the most variable features to use for integration
  integ_features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000) 
  #2.2. Prepare the SCT list object for integration
  seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = integ_features)
  #2.3. Find anchors
  integ_anchors <- FindIntegrationAnchors(object.list = seurat_objects, normalization.method = "SCT", 
                                          anchor.features = integ_features, dims = 1:param$npcs)
  #2.4. Integrate datasets
  seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", dims = 1:param$npcs)

  # switch to integrated assay for downstream analyses
  DefaultAssay(seurat_integrated) <- "integrated"
  # Run the standard workflow for visualization and clustering
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
  cm <- markers[ ,c("gene","cluster","pct.1", "pct.2", "avg_logFC","p_val_adj")]
  cm <- cm[cm$p_val_adj < 0.05, ]
  rownames(cm) <- NULL
  return(cm)
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

conservedMarkers <- function(scData) {
  markers <- list()
  for(eachCluster in levels(Idents(scData))){
    markersEach <- try(FindConservedMarkers(scData, ident.1=eachCluster, grouping.var="orig.ident", print.bar=FALSE, only.pos=TRUE,silent=TRUE))
    ## to skip some groups with few cells
    if(class(markersEach) != "try-error" && nrow(markersEach) > 0){
      markers[[eachCluster]] <- as_tibble(markersEach, rownames="gene")
    }
  }
  ## some of the cluster have no significant conserved markers
  markers <- markers[sapply(markers, nrow) != 0L] 
  markers <- bind_rows(markers, .id="cluster")
  
  return(markers)
}

diffExpressedGenes <- function(scData) {
  seurat_clusters <- Idents(scData)
  scData@meta.data$cluster.condition <- paste0(seurat_clusters, "_", scData@meta.data$orig.ident)
  Idents(scData) <- "cluster.condition"
  conditions <- unique(scData@meta.data$orig.ident)
  
  diffGenesFns <- c()
  conditionsComb <- combn(conditions, m=2)
  for(i in 1:ncol(conditionsComb)){
    diffGenes <- list()
    for(eachCluster in gtools::mixedsort(levels(seurat_clusters))){
      markersEach <- try(FindMarkers(scData, ident.1=paste0(eachCluster, "_",
                                                            conditionsComb[2,i]),
                                     ident.2=paste0(eachCluster, "_", 
                                                    conditionsComb[1,i]),
                                     print.bar=FALSE), silent=TRUE)
      ## to skip some groups with few cells
      if(class(markersEach) != "try-error"){
        diffGenes[[eachCluster]] <- as_tibble(markersEach, rownames="gene")
      }
    }
    diffGenes <- bind_rows(diffGenes, .id="cluster")
  }
  diffGenes <- diffGenes[diffGenes$p_val_adj < 0.05, ]
  
  return(diffGenes)
}

