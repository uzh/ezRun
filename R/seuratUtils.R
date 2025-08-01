###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

seuratStandardSCTPreprocessing <- function(scData, param, assay="RNA", seed=38) {
  DefaultAssay(scData) <- assay
  ## Get information on which variables to regress out in scaling/SCT
  vars.to.regress <- getSeuratVarsToRegress(param)
  ## generate normalized slots for the RNA assay
  scData <- NormalizeData(scData, normalization.method = "LogNormalize", scale.factor=10000, verbose=FALSE)
  species <- getSpecies(param$refBuild)
  if(ezIsSpecified(param$featSelectionMethod) && param$featSelectionMethod == 'STACAS'){
    require(STACAS)
    require(SignatuR)
    if(species == 'Human'){
      my.genes.blocklist <- c(GetSignature(SignatuR$Hs$Blocklists),
                              GetSignature(SignatuR$Hs$Compartments))
      scData <- FindVariableFeatures.STACAS(scData, nfeat = param$nfeatures, genesBlockList = my.genes.blocklist)
    } else if(species == 'Mouse'){
      my.genes.blocklist <- c(GetSignature(SignatuR$Mm$Blocklists),
                              GetSignature(SignatuR$Mm$Compartments))
      scData <- FindVariableFeatures.STACAS(scData, nfeat = param$nfeatures, genesBlockList = my.genes.blocklist)
    } else {
      message('Selection method STACAS not supported for this species! Use default method instead.')
      scData <- FindVariableFeatures(scData, selection.method = "vst", verbose = FALSE, nfeatures=param$nfeatures)  
    }
  } else {
    scData <- FindVariableFeatures(scData, selection.method = "vst", verbose = FALSE, nfeatures=param$nfeatures)
  }
  
  scData <- ScaleData(scData, vars.to.regress = vars.to.regress, verbose=FALSE, do.scale=FALSE)
  ## generate the SCT assay
  scData <- SCTransform(scData, vst.flavor="v2", vars.to.regress = vars.to.regress, assay = assay, seed.use = seed, verbose = FALSE,
                        return.only.var.genes=FALSE)
  return(scData)
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
  scData <- RunTSNE(object=scData, reduction.use = "pca", check_duplicates = FALSE,
                    dims.use=1:min(param$pcs, length(pc.genes)-1), tsne.method="Rtsne",
                    perplexity=ifelse(length(scData@ident) > 200, 30, 10),
                    num_threads=param$cores)
  try(scData <- RunUMAP(object=scData, reduction.use = "pca",
                        dims.use=1:min(param$pcs, length(pc.genes)-1),
                        n_neighbors=ifelse(length(scData@ident) > 200, 30, 10)),
      silent=TRUE)
  return(scData)
}

seuratStandardWorkflow <- function(scData, param, reduction="pca", ident.name="ident") {
  scData <- RunPCA(object=scData, verbose=FALSE)
  #scData <- RunPCA(object=scData, npcs = param$npcs, verbose=FALSE)
  if(!('Spatial' %in% as.vector(Seurat::Assays(scData)))){
    scData <- RunTSNE(object = scData, check_duplicates = FALSE, reduction = reduction, dims = 1:param$npcs)
  }
  scData <- RunUMAP(object=scData, reduction = reduction, dims = 1:param$npcs)
  scData <- FindNeighbors(object = scData, reduction = reduction, dims = 1:param$npcs, verbose=FALSE)
  myResolutions <- sort(unique(round(c(seq(from = 0.2, to = 1, by = 0.2), param$resolution), 1)))
  scData <- FindClusters(object=scData, resolution = myResolutions, verbose=FALSE)  #calculate clusters for a set of resolutions
  Idents(scData) <- scData@meta.data[,paste0(DefaultAssay(scData), "_snn_res.", param$resolution)]  #but keep as the current clusters the ones obtained with the resolution set by the user
  scData[[ident.name]] <- Idents(scData)
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
                 y=scDataList[-1], 
                 merge.data=TRUE)
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
  vars.to.regress <- getSeuratVarsToRegress(param)
  #1. Data preprocesing is done on each object separately
  for (i in 1:length(scDataList)) {
    scDataList[[i]] <- SCTransform(scDataList[[i]], vars.to.regress = vars.to.regress,  assay = assay, verbose = TRUE)
  }
  
  #2. Data integration
  #2.1. # Select the most variable features to use for integration
  integ_features <- SelectIntegrationFeatures(object.list = scDataList, nfeatures = param$nfeatures)
  
  if (param$integrationMethod %in% c("RPCA", "CCA")) {
    for (i in 1:length(scDataList)){
      scDataList[[i]] <- ScaleData(scDataList[[i]], features = integ_features)
    }
    #2.2. Prepare the SCT list object for integration
    scDataList <- PrepSCTIntegration(object.list = scDataList, anchor.features = integ_features)
    if(param$integrationMethod == 'RPCA'){
      scDataList <- lapply(X = scDataList, FUN = RunPCA, features = integ_features)
    }
    #2.3. Find anchors
    if(param$integrationMethod == 'CCA'){
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
  } else if (param$integrationMethod == "STACAS") {
    require(STACAS)
    #2.2 Merge normalized samples
    scDataList <- PrepSCTIntegration(object.list = scDataList, 
                                     anchor.features = integ_features)
    #2.3 Find anchor tree, either using prior label information or without
    if (ezIsSpecified(param$STACASAnnotationFile)) {
      stacas_anchors <- FindAnchors.STACAS(scDataList, 
                                           anchor.features = integ_features,
                                           cell.labels = "stacasLabelColumn",
                                           dims = 1:param$npcs)
      isSemisupervised <- TRUE
    } else {
      stacas_anchors <- FindAnchors.STACAS(scDataList, 
                                           anchor.features = integ_features,
                                           dims = 1:param$npcs)
      isSemisupervised <- FALSE
    }
    st1 <- SampleTree.STACAS(anchorset = stacas_anchors,
                             obj.names = sapply(scDataList, function(scData) {return(unique(scData$Sample))}))
    #2.4 Integration
    seurat_integrated <- IntegrateData.STACAS(stacas_anchors,
                                              sample.tree = st1,
                                              semisupervised=isSemisupervised,
                                              dims=1:param$npcs)
    seurat_integrated <- ScaleData(seurat_integrated)
    #3. Run the standard workflow for visualization and clustering
    seurat_integrated <- seuratStandardWorkflow(seurat_integrated, param)
  } else if (param$integrationMethod == "Harmony") {
    require(harmony)
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
    if (!is.null(param$additionalFactors)) {
      harmonyFactors <- 
        c("Condition", colnames(scData@meta.data)[startsWith(colnames(scData@meta.data), "meta_")])
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
  scData <- FindSpatiallyVariableFeatures(scData, features = VariableFeatures(scData), r.metric = 5, verbose = TRUE,
                                          selection.method = selection.method)
  spatialMarkers <- SpatiallyVariableFeatures(scData, selection.method = selection.method)
  #spatialMarkers <- SpatiallyVariableFeatures_workaround(scData, assay="SCT", selection.method = selection.method)
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

conservedMarkers <- function(scData, grouping.var="Condition", pseudoBulkMode=FALSE) {
  markers <- list()
  
  if("SCT" %in% Seurat::Assays(scData)) {
    assay <- "SCT"
  } else {
    assay <- "RNA"
  }
  if (pseudoBulkMode) {
    DE.method <- "DESeq2"
  } else {
    DE.method <- "wilcox"
  }
  
  for(eachCluster in levels(Idents(scData))){
    markersEach <- try(FindConservedMarkers(scData, ident.1=eachCluster, 
                                            grouping.var=grouping.var, 
                                            print.bar=FALSE, only.pos=TRUE, 
                                            assay = assay, test.use=DE.method), 
                       silent=TRUE)
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
  DE.method <- param$DE.method
  if (ezIsSpecified(param$replicateGrouping) && param$pseudoBulkMode == "true") {
    DE.method <- "DESeq2"
  }
  if (param$DE.method == "LR") { #regress the plate if the test is LR
    vars.to.regress <- param$DE.regress
  }
  
  diffGenes <- list()
  for(eachCluster in gtools::mixedsort(levels(seurat_clusters))){
    markersEach <- try(FindMarkers(scData, ident.1=paste0(eachCluster, "_",
                                                          param$sampleGroup),
                                   ident.2=paste0(eachCluster, "_", 
                                                  param$refGroup),
                                   test.use = DE.method, latent.vars = vars.to.regress))
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
    vars.to.regress <- c("CC.Difference")
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
  markers <- FindAllMarkers(
    object=scData, test.use = param$DE.method, only.pos=TRUE,
    min.pct=ifelse(ezIsSpecified(param$min.pct), param$min.pct, 0.1),
    logfc.threshold=ifelse(ezIsSpecified(param$logfc.threshold), param$logfc.threshold, 0.25)
  )
  ## Significant markers
  markers <- markers[ ,c("gene","cluster","pct.1", "pct.2", "avg_log2FC","p_val_adj")]
  markers <- markers[markers$p_val_adj < param$pvalue_allMarkers,]
  markers$cluster <- as.factor(markers$cluster)
  markers$diff_pct = abs(markers$pct.1-markers$pct.2)
  markers <- markers[markers$diff_pct >= param$min.diff.pct,]
  markers <- markers[order(markers$diff_pct, decreasing = TRUE),]
  return(markers)
}

getSeuratMarkersAndAnnotate <- function(scData, param, BPPARAM) {
  # function for general annotation of Seurat objects
  
  # Helper function to strip species labels from parameter values
  strip_species_labels <- function(param_values) {
    if (is.character(param_values)) {
      # Remove species labels like " (human)" or " (mouse)"
      return(gsub(" \\(human\\)| \\(mouse\\)", "", param_values))
    }
    return(param_values)
  }
  
  # Clean Azimuth and SingleR parameters to remove species labels
  if (ezIsSpecified(param$Azimuth)) {
    param$Azimuth <- strip_species_labels(param$Azimuth)
  }
  if (ezIsSpecified(param$SingleR)) {
    param$SingleR <- strip_species_labels(param$SingleR)
  }
  
  markers <- getSeuratMarkers(scData, param)
  
  # cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if (species == "Human" | species == "Mouse") {
    genesPerCluster <- split(markers$gene, markers$cluster)
    enrichRout <- querySignificantClusterAnnotationEnrichR(genesPerCluster, param$enrichrDatabase)
    cells.AUC <- cellsLabelsWithAUC(GetAssayData(scData, layer="counts"), species, param$tissue, BPPARAM = BPPARAM)
    singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, layer="data"), Idents(scData), param$SingleR, BPPARAM = BPPARAM)
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
    futile.logger::flog.info("Skipping pathway and TF activity")
  }
  
  
  # run Azimuth
  if (ezIsSpecified(param$Azimuth) && param$Azimuth != "none"){
    environment(MyDietSeurat) <- asNamespace('Seurat')
    assignInNamespace("DietSeurat", MyDietSeurat, ns = "Seurat")
    scDataAzi <- RunAzimuth(scData, param$Azimuth, assay="RNA") ## TODO support ADT
    
    ##Rename annotion levels if neccessary:
    colnames(scDataAzi@meta.data) <- sub('level_', 'l', colnames(scDataAzi@meta.data))
    
    aziNames <- setdiff(colnames(scDataAzi@meta.data), colnames(scData@meta.data))
    aziResults <- data.frame(
      Azimuth.celltype.l1=scDataAzi@meta.data[ , grep("l1$", aziNames, value=TRUE)],
      Azimuth.celltype.l2=scDataAzi@meta.data[ , grep("l2$", aziNames, value=TRUE)],
      Azimuth.celltype.l3=scDataAzi@meta.data[ , grep("l3$", aziNames, value=TRUE)],
      Azimuth.celltype.l4=scDataAzi@meta.data[ , grep("l4$", aziNames, value=TRUE)],
      row.names=colnames(scDataAzi))
    ## TODO: score should also be stored
    remove(scDataAzi)
  } else {
    aziResults <- NULL
  }
  
  
  # run cellxgene_annotation
  if (ezIsSpecified(param$cellxgeneUrl) && ezIsSpecified(param$cellxgeneLabel)){
    cellxgeneResults <- cellxgene_annotation(scData = scData,param = param)
    
  }else {
    cellxgeneResults <- NULL
  }
  
  
  
  return(list(markers=markers, cells.AUC=cells.AUC, singler.results=singler.results,
              enrichRout=enrichRout, pathwayActivity=pathwayActivity, TFActivity=TFActivity,
              aziResults=aziResults, cellxgeneResults=cellxgeneResults))
}

##' @title Merge Seurat clusters
##' @description Given an input mapping, group Seurat clusters of the active Seurat Ident together, i.e. merging them
##' @slot scData the Seurat object
##' @slot clustList a list() object where the names are the new cluster names and their values are all the clusters to be merged. If the list consists of multiple sets of clusters to be merged, ensure the sets do not overlap.
##' @template roxygen-template
##' @examples
##' data("pbmc_small")
##' clustList <- list("foo"=as.character(c(1, 2)))
##' seuratMergeClusters(pbmc_small, clustList)
seuratMergeClusters <- function(scData, clustList) {
  require(stringr)
  assertthat::assert_that(
    length(clustList) <= 1 || length(Reduce(intersect, clustList))==0, msg="Sets of input clusters overlap!"
  )

  clusters <- as.character(unique(Idents(scData)))
  newIdent <- as.character(unname(Idents(scData)))
  for (targetClust in names(clustList)) {
    clustersToChange <- clustList[[targetClust]]
    targetClustRep <- rep(targetClust, length(clustersToChange))
    names(targetClustRep) <- clustersToChange
    newIdent <- ifelse(newIdent %in% clustersToChange, unname(targetClustRep[clustersToChange]), newIdent)
  }
  newIdent <- factor(newIdent, levels=str_sort(unique(newIdent), numeric=TRUE))
  return(newIdent)
}
