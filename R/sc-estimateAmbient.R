library(DropletUtils)
library(data.table)
library(ggpubr)
library(devtools)
library(data.table)
library(scater)
library(cluster)
library(celda)
library(ggplot2)
library(dbscan)
library(SoupX)
library(parallel)
library(doParallel)
library(Seurat)
library(scDblFinder)
library(SingleR)

library(BiocParallel)
BPPARAM=MulticoreParam(workers=8)


getClusterDBSCAN <- function(umapData, dbscan.eps, dbscan.MinPts){
  totalClusters = 1
  while(totalClusters <= 1 & dbscan.eps > 0) {
    clusterLabels = dbscan(umapData, dbscan.eps, minPts = dbscan.MinPts)
    dbscan.eps = dbscan.eps * 0.75
    totalClusters = length(unique(clusterLabels$cluster))
    
    if (dbscan.eps < 0.05){
      stop("ddbscan finds a single cluster!") 
    }
    
  }
  return(clusterLabels$cluster)
}

autoEstContTfidfMin <- function(sc, tfidfMin){
  if (tfidfMin < 0){
    stop("Parameter tfidfMin cannot be less than 0!") 
  }
  tryCatch({
    sc <- autoEstCont(sc, tfidfMin=tfidfMin, forceAccept=T, doPlot=FALSE)
  }, 
  error=function(cond) {
    tfidfMin = tfidfMin -  0.3
    autoEstContTfidfMin(sc, tfidfMin)
  })
  return (sc)
}

checkAndCleanAntibody <- function(object){
  if (is.list(object)){
    object <- object$`Gene Expression`
  }
  return (object)
}

filterLowQualityCells <- function(object, sampleName, minGenes=400, minCounts=4000){
  
  if (class(object) == "dgCMatrix"){
    matrixObject <- object
  }else{
    matrixObject <- counts(object) # get matrix
  }
  
  scData <- CreateSeuratObject(matrixObject)
  scData <- PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
  
  isHighMito <- isOutlier(scData$percent_mito, type="higher")
  isLowGene <- scData$nFeature_RNA < minGenes
  isLowCount <- scData$nCount_RNA < minCounts
  isLowQuality <- isHighMito | isLowGene | isLowCount
  scData <- scData[ , !isLowQuality]
  
  scData <- NormalizeData(scData, normalization.method = "LogNormalize", scale.factor = 10000)
  scData <- ScaleData(scData, features = rownames(scData))
  scData <- FindVariableFeatures(scData, selection.method = "vst")
  scData <- RunPCA(scData, features = VariableFeatures(object = scData))
  scData <- FindNeighbors(object = scData, dims = 1:30)
  plotFeatureCountRNAPlot(scData, 
                          sampleName, 
                          isLowQuality)
  
  
  scData <- FindClusters(scData)
  scDbl <- scDblFinder(GetAssayData(scData, "counts", assay="RNA"), clusters=scData$seurat_clusters, returnType = "table")
  singletCells <- rownames(scDbl)[ scDbl$class == "singlet"]
  return(object[, colnames(object) %in% singletCells])
}



#cleanCellsWrapper <- function(firstObeject, secondObject, sampleName){
#  
#  firstObeject <- checkAndCleanAntibody(firstObeject)
#  filteredObjects <- filterLowQualityCells(firstObeject, sampleName)
#  firstObeject <- filteredObjects$object
#  clusters <- filteredObjects$cluster
#  nFeature_RNA=filteredObjects$nFeature_RNA
#  nCount_RNA=filteredObjects$nCount_RNA
#  
#  secondObject <- checkAndCleanAntibody(secondObject)
#  
#  secondObject <- secondObject[, colnames(firstObeject) %in% colnames(secondObject)]
#  
#  return (list(firstObeject=firstObeject, 
#               secondObject=secondObject, 
#               clusters=clusters,
#               nFeature_RNA=nFeature_RNA,
#               nCount_RNA=nCount_RNA
#  ))
#  
#}


plotFeatureCountRNAPlot <- function(scData, sampleName, isLowQuality){
  
  percLowQuality <- round(sum(isLowQuality) / length(isLowQuality), 2) * 100
  
  featureCountRNAPlot <- (FeaturePlot(scData, features = "nCount_RNA") 
                          + ggtitle(paste0(sampleName, " - nCount_RNA \n % lowQualityCells = ", 
                                           percLowQuality, 
                                           "% \n LowQualityCell =  HighMito | LowGene | LowCount")
                          ))
  
  ggsave(paste0("plots/featureCountRNAPlot/", sampleName, ".png"), featureCountRNAPlot)
}

shrinkToRange <- function(x, theRange) {
  x[x > theRange[2]] <- theRange[2]
  x[x < theRange[1]] <- theRange[1]
  return(x)
}

contaminationPlot <- function(sceClean, sce, method, sampleName){
  
  clusters <- sce$clusters
  nFeature_RNA <- colSums(counts(sce) > 0)
  nCount_RNA <- colSums(counts(sce))
  
  
  tab <- data.table(contaminationFraction=shrinkToRange(sceClean$contaminationFraction, c(0.01,1)),
                    cluster=as.factor(clusters),
                    nFeature_RNA=nFeature_RNA,
                    nCount_RNA=nCount_RNA
  )
  
  contaminationVSnRNAPlot <- (ggplot(melt(tab, id.vars = c("contaminationFraction", "cluster")), aes(value, contaminationFraction, color=cluster))
                              + geom_point()
                              + facet_wrap(~variable, nrow=2, scales = "free_x")
                              + scale_x_log10() 
                              + ylim(0, 1)
                              # + scale_y_log10()
                              + theme_bw()
                              + ggtitle(paste0(method, " - ", sampleName))
  )
  
  ggsave(paste0("plots/contaminationVSnRNAPlot/", sampleName, "_",method, ".png"), contaminationVSnRNAPlot)
  
  contaminationViolin <- (ggplot(tab, aes(cluster, contaminationFraction, color=cluster))
                          + geom_jitter(alpha=0.1, size=0.4)
                          + geom_violin() 
                          + ylim(0, 1)
                          + theme_bw()
                          + ggtitle(paste0(method, " - ", sampleName))
  )
  
  ggsave(paste0("plots/contaminationViolin/", sampleName, "_",method, ".png"), contaminationViolin)
  
}



calculateEntropy <- function(sce, species){
  
  sce <- scater::logNormCounts(sce)
  
  ref <- switch(species, 
                "Human"=celldex::HumanPrimaryCellAtlasData(),
                "Mouse"=celldex::MouseRNAseqData(),
                NULL)
  
  if (!is.null(ref)){
    singler <- SingleR(test = logcounts(sce), ref = assay(ref, "logcounts"), aggr.ref = TRUE,
                       labels = ref$label.main, clusters = NULL,
                       de.method = "classic", BPPARAM = BPPARAM)
    
    # if negative correlation, set to 0.001
    # ignore how negative cluster correlation is
    # otherwise problem due to log in the entropy formula
    singler$scores[singler$scores < 0] <- 0.001
    
    cellEntropy <- apply(singler$scores, 1, function(x){TFBSTools::shannon.entropy(x/sum(x))})
  } else {
    cellEntropy <- NULL
  }
  
  return(cellEntropy)
  
}

plotViolinBeforeAfterDiffViolin <- function(listData, sampleName, clusters, plotType){
  
  tableData <- as.data.table(listData)
  tableData[, cluster:= as.factor(clusters)]
  tableData <- melt(tableData, id.vars = c("orig", "cluster"), value.name = "After", variable.name = "method")
  tableData[, Before:= orig]
  
  tableData[, diff:= Before - After]
  
  plotViolin <- (ggplot(tableData, aes(cluster, diff, color=cluster))
                 + geom_jitter(alpha=0.1, size=0.4)
                 + geom_violin(alpha=0.7)
                 + ylab("Before - After")
                 + theme_classic()
                 + facet_wrap(~method, nrow=2)
                 + ggtitle(paste0(sampleName, " - ", plotType))
  )
  
  ggsave(paste0("plots/", plotType, "/", sampleName, ".png"), plotViolin)
  
}


plotBeforeAfterScatter <- function(listData, sampleName, clusters, plotType){
  
  tableData <- as.data.table(listData)
  tableData[, cluster:= as.factor(clusters)]
  tableData <- melt(tableData, id.vars = c("orig", "cluster"), value.name = "After", variable.name = "method")
  tableData[, Before:= orig]
  
  
  
  plotScatter <- (ggplot(tableData, aes(Before, After)) 
                  + geom_point(aes(color=cluster), size=0.5, alpha=0.5)
                  + geom_smooth(method="lm", alpha=0.5)
                  + theme_classic()
                  + geom_abline()
                  + facet_wrap(~method)
                  + ggtitle(paste0(sampleName, " - ", plotType))
  )
  
  ggsave(paste0("plots/", plotType, "/", sampleName, ".png"), plotScatter)
  
  
}


plotBar <- function(listData, sampleName, plotType){
  
  tableData <- as.data.table(listData)
  tableData <- melt(tableData, id.vars = "orig", value.name = "After", variable.name = "method")
  
  baseline <- mean(tableData[, orig])
  
  plot <- (ggplot(tableData, aes(method, After)) 
           + geom_boxplot(alpha=0.1)
           + geom_violin(alpha=0.7)
           + geom_hline(yintercept = baseline, linetype="dashed", color="red", alpha=0.5)
           + geom_text(aes(0, baseline, label = "baseline", vjust = -1, hjust=-0.5), color="red")
           + ggtitle(paste0(sampleName, " - ", plotType))
  )
  
  ggsave(paste0("plots/", plotType, "/", sampleName, ".png"), plot)
  
  
}

#getOnlyHighlyVariableGenes <- function(sce, nGenes=2000){
#  # from https://jdblischak.github.io/singlecell-qtl/pca-variable.html
#  
#  # CPM
#  log2cpm <- edgeR::cpm(exprs(sce), log = TRUE)
#  
#  # Calculate coefficient of variation.
#  cv <- apply(log2cpm, 1, function(x) sd(x) / mean(x))
#  
#  # should we use hard cutoff? or take top n %
#  topGenNames <- names(cv[order(-cv)][1:nGenes])
#  return  (sce[topGenNames, ])
#  
#}

plotExplainedVariance <- function(sce, method, sampleName){
  percentExplainedVar <- cumsum(attr(reducedDim(sce, "PCA"), "percentVar"))
  
  explainedVariance <- (qplot(1:length(percentExplainedVar), percentExplainedVar) 
                        + ylim(0, 100)
                        + ggtitle(paste0(method, " - ", sampleName))
                        + theme_classic()
  )
  
  ggsave(paste0("plots/explainedVariance/", sampleName, "_",method, ".png"), explainedVariance)
  
  
}

clusterSce <- function(sce, dbscan.eps, dbscan.MinPts, nTopVariableGenes=2000){
  sce <- scater::logNormCounts(sce)
  ## TODO should we really use UMAP? or rather plain PCA????
  ## COMMENT: the number of comments should be > 2 for UMAP, e.g. 5, and something like > 10 for PCA
  #sceHighlyVariableGenes <- getOnlyHighlyVariableGenes(sce, nGenes=2000)
  
  # NB: ncomponents = 200
  sce <- runPCA(sce, ntop=nTopVariableGenes, ncomponents=200)
  percentVar <- attr(reducedDim(sce, "PCA"), "percentVar")
  nComponentsPCA_90perc <- sum(cumsum(percentVar) < 90) # reduce noise
  sce <- runUMAP(sce, ntop=nTopVariableGenes, ncomponents=5, pca=nComponentsPCA_90perc)
  sce$clusters <- getClusterDBSCAN(reducedDim(sce, "UMAP"), dbscan.eps, dbscan.MinPts)
  return(sce)
}

calculateSilhouette <- function(sce){
  dis <- dist(reducedDim(sce), method = "euclidean")
  return(silhouette(as.integer(sce$clusters), dis)[ ,"sil_width"])
}

computeDecontaminated <- function(method, sampleName, sce, rawFeatureDir=NULL, threads=1){
  library(celda)
  
  cleanCounts <- switch(method,
                        "DecontX"={
                          sce <- decontX(sce, z=sce$clusters)
                          decontXcounts(sce)
                        },
                        "SoupX"={
                          if (is.null(rawFeatureDir) || !file.exists(rawFeatureDir)){
                            NULL
                          } else {
                            tod <- checkAndCleanAntibody(Seurat::Read10X(rawFeatureDir, gene.column = 1))
                            featInfo <- ezRead.table(paste0(rawFeatureDir, "/features.tsv.gz"), header = FALSE, row.names = NULL)#, col_names = FALSE)
                            colnames(featInfo) <- c("ensemblID", "name", "type")
                            featInfo <- featInfo[match(rownames(tod), featInfo$ensemblID), ]
                            rownames(tod) <- gsub("_", "-", uniquifyFeatureNames(ID=featInfo$ensemblID, names=featInfo$name))
                            commGenes <- intersect(rownames(tod), rownames(counts(sce)))
                            
                            sc <- SoupChannel(tod[commGenes, ], counts(sce)[commGenes, ])
                            sc <- setClusters(sc, sce$clusters)
                            sc <- autoEstContTfidfMin(sc, tfidfMin=1)
                            out <- adjustCounts(sc)
                            out
                          }
                        }
                        
  )
  return(cleanCounts)
}

getMetrics <- function(index, dbscan.eps=1, dbscan.MinPts=5, threads=1){
  
  tryCatch({
    
    sampleName <-strsplit(filteredFeatureDirList[index], "/")[[1]][7]
    cts <- checkAndCleanAntibody(Seurat::Read10X(filteredFeatureDirList[index]))
    
    ### keep only those genes that are expressed in more than one cell
    #cts <- cts[rowSums2(cts >0) > 1, ]
    
    cts <- filterLowQualityCells(cts, sampleName)
    
    
    sce <- SingleCellExperiment(assays=list(counts=cts))
    sce <- clusterSce(sce, dbscan.eps, dbscan.MinPts)
    plotExplainedVariance(sce, "rawCounts", sampleName)
    if (length(unique(sce$clusters)) == 1){
      stop("Detected unique clusters  == 1!") 
    }
    
    
    silhouetteList <- list()
    entropyList <- list()
    
    
    dataset <- ezRead.table(file.path(filteredFeatureDirList[index], "../../dataset.tsv"))
    species <- getSpecies(dataset$refBuild[1])
    
    
    
    entropyList[["orig"]] <-  calculateEntropy(sce, species)
    silhouetteList[["orig"]] <- calculateSilhouette(sce)
    for (method in c("SoupX", "DecontX")){
      ctsClean <- computeDecontaminated(method=method, 
                                        sampleName, 
                                        sce, 
                                        rawFeatureDir=rawFeatureDirList[index],
                                        threads=threads)
      
      sceClean <- SingleCellExperiment(assays=list(counts=ctsClean))
      sceClean <- clusterSce(sceClean, dbscan.eps, dbscan.MinPts)
      
      sceClean$contaminationFraction <- (colSums(counts(sce)) - colSums(counts(sceClean))) / colSums(counts(sce))
      
      
      contaminationPlot(sceClean,
                        sce,
                        method, sampleName)
      
      silhouetteList[[method]] <- calculateSilhouette(sceClean)
      entropyList[[method]] <-  calculateEntropy(sceClean, species)
      
      saveRDS(sceClean, file = paste0("raw_data/", sampleName, "_", method, "_sce.rds"))
      
      
    }
    
    plotBeforeAfterScatter(entropyList, sampleName, sce$clusters, "entropyScatterPlot")
    plotBeforeAfterScatter(silhouetteList, sampleName, sce$clusters,"silhouetteScatterPlot")
    
    plotViolinBeforeAfterDiffViolin(entropyList, sampleName, sce$clusters, "entropyViolinPlot")
    plotViolinBeforeAfterDiffViolin(silhouetteList, sampleName, sce$clusters,"silhouetteViolinPlot")
    
    
    plotBar(silhouetteList, sampleName, "silhouetteBarPlot")
    
    metrics <- data.table(
      method=c("Original", "DecontX", "SoupX"),
      meanSilWidth=sapply(silhouetteList, mean),
      sampleName=sampleName
    )
    return(metrics)
  },
  error=function(cond) {
    # has to be removed as it tells the worker about the error
    #message(paste0("Please look manually file: ", filteredFeatureDir))
    #message("Here's the original error message:")
    #message(cond)
    # Choose a return value in case of error
    
    fwrite(list(paste(filteredFeatureDirList[index], cond, sep = " - ")), 
           "log/logFailedSamples.csv", 
           append = T)
    
    metrics <- data.table(
      method=c("Original", "DecontX", "SoupX"),
      intercept=NA,
      slope=NA,
      meanSilWidth=NA,
      sampleName=sampleName
    )
    
    return(metrics)
  })
  
}


addAmbientEstimateToSeurat <- function(scData, rawDir=NULL, threads=1){
  library(celda)
  library(SoupX)
  
  sce <- SingleCellExperiment(assays=list(counts=GetAssayData(scData, slot="counts", assay="RNA")), colData = scData@meta.data)
  sce$clusters <- sce$ident
  for (method in c("SoupX", "DecontX")){
    ctsClean <- computeDecontaminated(method=method, 
                                      sampleName=NULL, 
                                      sce, 
                                      rawFeatureDir=rawDir,
                                      threads=threads)
    if (!is.null(ctsClean)){
      contaminationFraction <- (colSums2(counts(sce)) - colSums2(ctsClean)) / colSums2(counts(sce))
      scData <- AddMetaData(scData, metadata = contaminationFraction, col.name = paste0(method, "_contFrac"))
    }
  }
  return(scData)
}

