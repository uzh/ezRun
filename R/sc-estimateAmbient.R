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

addAmbientEstimateToSeurat <- function(scData, rawDir=NULL, threads=1){
  library(celda)
  library(SoupX)
  
  sce <- SingleCellExperiment(assays=list(counts=GetAssayData(scData, layer="counts")), colData = scData@meta.data)
  sce$clusters <- Idents(scData)
  ## decontx
  ctsClean <- decontX(sce, z=sce$clusters)
  ctsClean <- decontXcounts(ctsClean)
  contaminationFraction <- (colSums2(counts(sce)) - colSums2(ctsClean)) / colSums2(counts(sce))
  scData <- AddMetaData(scData, metadata = contaminationFraction, col.name = paste0("DecontX_contFrac"))
  ## SoupX
  if (!is.null(rawDir) && file.exists(rawDir)){
    tod <- checkAndCleanAntibody(Seurat::Read10X(rawDir, gene.column = 1))
    featInfo <- ezRead.table(paste0(rawDir, "/features.tsv.gz"), header = FALSE, row.names = NULL)#, col_names = FALSE)
    colnames(featInfo) <- c("ensemblID", "name", "type")
    featInfo <- featInfo[featInfo$type=='Gene Expression',]
    rownames(tod) <- gsub("_", "-", uniquifyFeatureNames(ID=featInfo$ensemblID, names=featInfo$name))
    commGenes <- intersect(rownames(tod), rownames(counts(sce)))
    sc <- SoupChannel(tod[commGenes, ], counts(sce)[commGenes, ])
    sc <- setClusters(sc, sce$clusters)
    #try({sc1 <- autoEstCont(sc, tfidfMin=1, forceAccept=T, doPlot=FALSE)})
    sc <- autoEstContTfidfMin(sc, tfidfMin=1)
    if (!is.null(sc$fit) && "rho" %in% colnames(sc$metaData)){
      ctsClean <- adjustCounts(sc) ## NOTE: ctsClean might have less genes than sce
      contaminationFraction <- (colSums2(counts(sce)) - colSums2(ctsClean)) / colSums2(counts(sce))
      scData <- AddMetaData(scData, metadata = contaminationFraction, col.name = paste0("SoupX_contFrac"))
    }
  }
  return(scData)
}

