checkAndCleanAntibody <- function(object){
    if (is.list(object)){
        object <- object$`Gene Expression`
    }
    return (object)
}

autoEstContTfidfMin <- function(sc, tfidfMin){
    while (tfidfMin > 0){
        scOut  <- try({ autoEstCont(sc, tfidfMin=tfidfMin, forceAccept=T, doPlot=FALSE)})
        if (length(class(scOut)) == 1 && class(scOut) == "try-error"){
            tfidfMin = tfidfMin -  0.3
        } else {
            break
        }
    }
    return(scOut)
}

addAmbientEstimateToSeurat <- function(scData, rawDir=NULL, param=NULL) {
  library(celda)
  library(SoupX)
  
  threads <- param$cores
  sce <- SingleCellExperiment(assays=list(counts=GetAssayData(scData, layer="counts")), colData = scData@meta.data)
  sce$clusters <- Idents(scData)
  ## decontx
  ctsClean <- decontX(sce, z=sce$clusters)
  ctsClean <- decontXcounts(ctsClean)
  contaminationFraction <- (colSums2(counts(sce)) - colSums2(ctsClean)) / colSums2(counts(sce))
  scData <- AddMetaData(scData, metadata = contaminationFraction, col.name = paste0("DecontX_contFrac"))
  ## SoupX
  if (!is.null(rawDir) && file.exists(rawDir)){
    if(param$cellbender){  
        tod <- checkAndCleanAntibody(Seurat::Read10X_h5(file.path(dirname(rawDir),'cellbender_filtered_seurat.h5'), use.names = FALSE))
        
        # Use explicit features path if available
        if(!is.null(param$featuresPath) && file.exists(param$featuresPath)) {
            featuresFile <- param$featuresPath
            featInfo <- ezRead.table(featuresFile, header = FALSE, row.names = NULL)
        } else {
            # If featuresPath not provided, try standard locations
            countMatrixToUse <- if(dirname(param$cellrangerCountFiltDir) != dirname(param$cellrangerCountRawDir)) {
                param$cellrangerCountFiltDir
            } else {
                param$cellrangerCountRawDir
            }
            featuresFile <- file.path(countMatrixToUse, "features.tsv.gz")
            
            # If standard location doesn't work, add sample_filtered_feature_bc_matrix subdirectory
            if(!file.exists(featuresFile) && grepl("CellRangerMulti", countMatrixToUse)) {
                featuresFile <- file.path(countMatrixToUse, "sample_filtered_feature_bc_matrix", "features.tsv.gz")
            }
            
            if(file.exists(featuresFile)) {
                featInfo <- ezRead.table(featuresFile, header = FALSE, row.names = NULL)
            } else {
                warning(paste0("Could not find features.tsv.gz file at: ", featuresFile))
                # Create empty featInfo to allow processing to continue
                featInfo <- data.frame(V1=character(0), V2=character(0), V3=character(0))
            }
    colnames(featInfo) <- c("ensemblID", "name", "type")
    featInfo <- featInfo[featInfo$type=='Gene Expression',]
    rownames(tod) <- gsub("_", "-", uniquifyFeatureNames(ID=featInfo$ensemblID, names=featInfo$name))
    commGenes <- intersect(rownames(tod), rownames(counts(sce)))
    sc <- SoupChannel(tod[commGenes, ], counts(sce)[commGenes, ])
    sc <- setClusters(sc, sce$clusters)
    #try({sc1 <- autoEstCont(sc, tfidfMin=1, forceAccept=T, doPlot=FALSE)})
    sc <- autoEstContTfidfMin(sc, tfidfMin=1)
    if(length(class(sc)) == 2L){
      if (!is.null(sc$fit) && "rho" %in% colnames(sc$metaData)){
        ctsClean <- adjustCounts(sc) ## NOTE: ctsClean might have less genes than sce
        contaminationFraction <- (colSums2(counts(sce)) - colSums2(ctsClean)) / colSums2(counts(sce))
        scData <- AddMetaData(scData, metadata = contaminationFraction, col.name = paste0("SoupX_contFrac"))
      }
    }
  }
  return(scData)
}

