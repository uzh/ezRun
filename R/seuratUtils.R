###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

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

seuratStandardWorkflow <- function(scData, param, reduction="pca") {
  scData <- RunPCA(object=scData, npcs = param$npcs, verbose=FALSE)
  if(!('Spatial' %in% as.vector(Seurat::Assays(scData)))){
    scData <- RunTSNE(object = scData, reduction = reduction, dims = 1:param$npcs)
  }
  scData <- RunUMAP(object=scData, reduction = reduction, dims = 1:param$npcs)
  scData <- FindNeighbors(object = scData, reduction = reduction, dims = 1:param$npcs, verbose=FALSE)
  scData <- FindClusters(object=scData, resolution = seq(from = 0.2, to = 1, by = 0.2), verbose=FALSE)  #calculate clusters for a set of resolutions
  Idents(scData) <- scData@meta.data[,paste0(DefaultAssay(scData), "_snn_res.", param$resolution)]  #but keep as the current clusters the ones obtained with the resolution set by the user
  scData$ident <- Idents(scData)
  return(scData)
}  

seuratClusteringV3 <- function(scData, param, assay="RNA") {
  vars.to.regress <- getSeuratVarsToRegress(param)
  
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
  vars.to.regress <- getSeuratVarsToRegress(param)
  #1. Data preprocesing is done on each object separately
  for (i in 1:length(scDataList)) 
    scDataList[[i]] <- SCTransform(scDataList[[i]], vars.to.regress = vars.to.regress,  assay = assay, verbose = TRUE)
  #2. Merge all seurat objects
  scData = merge(x=scDataList[[1]], 
                 y=tail(scDataList, length(scDataList) - 1), 
                 merge.data=TRUE)
  VariableFeatures(scData) <- unlist(lapply(scDataList, VariableFeatures))
  #3. Dimensionality reduction and clustering
  scData <- seuratStandardWorkflow(scData, param)
  
  return(scData)
}

cellClustWithCorrection <- function (scDataList, param) {
  require(harmony)
  
  if(param[['name']] == 'SpatialSeuratSlides')
    assay = "Spatial"
  else 
    assay= "RNA"
  vars.to.regress <- getSeuratVarsToRegress(param)
  #1. Data preprocesing is done on each object separately
  for (i in 1:length(scDataList)) {
    scDataList[[i]] <- SCTransform(scDataList[[i]], vars.to.regress = vars.to.regress,  assay = assay, verbose = TRUE)
  }
  
  #2. Data integration
  #2.1. # Select the most variable features to use for integration
  integ_features <- SelectIntegrationFeatures(object.list = scDataList, nfeatures = 3000)
  
  if (param$integrationMethod %in% c("RPCA", "Classic")) {
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
  } else if (param$integrationMethod == "Harmony") {
    #2.2 Merge normalized samples
    scData <- merge(x = scDataList[[1]],
                    y = scDataList[2:length(scDataList)],
                    merge.data = TRUE)
    #2.3.1 Manually set variable features of merged Seurat object
    VariableFeatures(scData) <- integ_features
    #2.3.2 Calculate PCs using manually set variable features
    scData <- RunPCA(scData, assay = "SCT", npcs = param$npcs)
    
    #2.4 Prep and run Harmony algorithm
    # Find the additional harmony factors if we have any
    if (!is.null(param$harmonyFactors)) {
      harmonyFactors <- 
        c("Condition", colnames(scData@meta.data)[startsWith(colnames(scData@meta.data), "har_")])
    } else {
      harmonyFactors <- "Condition"
    }
    # Calculate Harmony reduction
    scData <- RunHarmony(scData, 
                         group.by.vars = harmonyFactors, 
                         reduction = "pca", assay.use = "SCT",
                         reduction.save = "harmony")
    
    #3. Run the standard workflow for visualization and clustering
    seurat_integrated <- seuratStandardWorkflow(scData, param, reduction="harmony")
  }
  
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



SpatiallyVariableFeatures_workaround <- function(object, assay="SCT", selection.method = "moransi") {
  #' This is work around function to replace SeuratObject::SpatiallyVariableFeatures function.
  #' return ranked list of Spatially Variable Features
  
  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }
  
  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features
  
  # Select columns starting with the provided col_prefix
  moransi_cols <- grep(paste0("^", selection.method), colnames(data), value = TRUE)
  
  # Filter rows where "moransi.spatially.variable" is TRUE
  filtered_data <- data[data[[paste0(selection.method, ".spatially.variable")]], moransi_cols]
  
  # Sort filtered data by "moransi.spatially.variable.rank" column in ascending order
  sorted_data <- filtered_data[order(filtered_data[[paste0(selection.method, ".spatially.variable.rank")]]), ]
  sorted_data <- sorted_data[grep('^NA', rownames(sorted_data), invert = TRUE),]
  sorted_data[['Rank']] = 1:nrow(sorted_data)
  # Return row names of the sorted data frame
  return(sorted_data)
}

spatialMarkers <- function(scData, selection.method = "markvariogram") { 
  scData <- FindSpatiallyVariableFeatures(scData, features = VariableFeatures(scData), r.metric = 5, 
                                          selection.method = selection.method)
  #spatialMarkers <- SpatiallyVariableFeatures(scData, selection.method = "markvariogram") #deactivated due to an unfixed bug in the current Seurat version 4.3
  spatialMarkers <- SpatiallyVariableFeatures_workaround(scData, assay="SCT", selection.method = selection.method)
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
  markers$avg_avg_fc <- markers %>% dplyr::select(contains("_avg_log2FC")) %>% rowMeans()
  return(markers)
}

diffExpressedGenes <- function(scData, param, grouping.var="Condition") {
  seurat_clusters <- Idents(scData)
  scData@meta.data$cluster.condition <- 
    paste0(seurat_clusters, "_", scData@meta.data[[grouping.var]])
  Idents(scData) <- "cluster.condition"
  conditions <- unique(scData@meta.data[[grouping.var]])
  
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
  
  diffGenes$diff_pct = abs(diffGenes$pct.1-diffGenes$pct.2)
  rownames(diffGenes) <- NULL
  # diffGenes <- diffGenes[diffGenes$p_val_adj < 0.05, ]
  
  return(diffGenes)
}

getSeuratVarsToRegress <- function(param) {
  vars.to.regress <- NULL
  if (ezIsSpecified(param$SCT.regress.CellCycle) && 
      param$SCT.regress.CellCycle) {
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")
  }
  if (!is.null(param$SCT.regress.var)) {
    vars.to.regress <- c(vars.to.regress, param$SCT.regress.var)
  }
  return(vars.to.regress)
}

getSeuratMarkers <- function(scData, param) {
  # positive cluster markers
  ## https://github.com/satijalab/seurat/issues/5321
  ## https://github.com/satijalab/seurat/issues/1501
  markers <- FindAllMarkers(object=scData, test.use = param$DE.method, only.pos=TRUE)
  ## Significant markers
  markers <- markers[ ,c("gene","cluster","pct.1", "pct.2", "avg_log2FC","p_val_adj")]
  markers$cluster <- as.factor(markers$cluster)
  markers$diff_pct = abs(markers$pct.1-markers$pct.2)
  markers <- markers[order(markers$diff_pct, decreasing = TRUE),]
  return(markers)
}

getSeuratMarkersAndAnnotate <- function(scData, param) {
  # function for general annotation of Seurat objects
  markers <- getSeuratMarkers(scData, param)
  
  # cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if (species == "Human" | species == "Mouse") {
    genesPerCluster <- split(markers$gene, markers$cluster)
    enrichRout <- querySignificantClusterAnnotationEnrichR(genesPerCluster, param$enrichrDatabase)
    cells.AUC <- cellsLabelsWithAUC(GetAssayData(scData, "counts"), species, param$tissue, BPPARAM = BPPARAM)
    singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "data"), Idents(scData), param$SingleR, BPPARAM = BPPARAM)
    for (r in names(singler.results)) {
      scData[[paste0(r,"_single")]] <- singler.results[[r]]$single.fine$labels
      scData[[paste0(r,"_cluster")]] <- singler.results[[r]]$cluster.fine$labels[match(Idents(scData), rownames(singler.results[[r]]$cluster.fine))]
    }
  } else {
    cells.AUC <- NULL
    singler.results <- NULL
    enrichRout <- NULL
  }
  
  ## SCpubr advanced plots, can currently only be computed for human and mouse
  if (ezIsSpecified(param$computePathwayTFActivity) && 
      as.logical(param$computePathwayTFActivity) &&
      (species == "Human" | species == "Mouse")) {
    pathwayActivity <- computePathwayActivityAnalysis(cells = scData, species = species)
    TFActivity <- computeTFActivityAnalysis(cells = scData, species = species)
  } else {
    pathwayActivity <- NULL
    TFActivity <- NULL
    print("Skipping pathway and TF activity")
  }
  
  return(list(markers=markers, cells.AUC=cells.AUC, singler.results=singler.results,
              enrichRout=enrichRout, pathwayActivity=pathwayActivity, TFActivity=TFActivity))
}
