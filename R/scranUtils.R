###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

scranIntegration <- function(sceList, method=c("None", "MNN")){
  require(scran)
  require(scater)
  require(batchelor)
  require(SingleCellExperiment)
  require(BiocSingular)
  
  method <- match.arg(method)
  
  # Feature selection across batches
  decList <- lapply(sceList, function(x){metadata(x)$dec})
  decList <- lapply(decList, function(x){x[order(x$bio, decreasing = TRUE), ]})
  universe <- Reduce(intersect, lapply(decList, rownames))
  bioMatrix <- do.call(cbind, lapply(decList, function(x){x[universe, "bio"]}))
  mean.bio <- rowMeans(bioMatrix)
  chosen <- universe[mean.bio > 0]
  rescaled <- do.call(batchelor::multiBatchNorm,
                      lapply(sceList, function(x){x[universe, ]}))
  
  # We always do the simple integration. Even just for comparison
  if(method == "None"){
    sce <- SingleCellExperiment(
      assays=list(counts=do.call(cbind, lapply(rescaled, assay, "counts")),
                  logcounts=do.call(cbind, lapply(rescaled, assay, "logcounts"))),
      rowRanges=rowRanges(rescaled[[1]]),
      colData=do.call(rbind, lapply(rescaled, colData)),
      metadata=list(param=metadata(sceList[[1]])$param,
                    log.exprs.offset=metadata(rescaled[[1]])$log.exprs.offset)
    )
    new.trend <- makeTechTrend(x=sce)
    if(toupper(metadata(sce)$param$scProtocol) == "10X"){
      fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.05))
    }else if(toupper(metadata(sce)$param$scProtocol) == "SMART-SEQ2"){
      fit <- trendVar(sce, use.spikes=FALSE, loess.args=list(span=0.3))
    }
    fit$trend <- new.trend # overwrite trend.
    dec <- decomposeVar(fit=fit) # use per-gene variance estimates in 'fit'.
    metadata(sce)$dec <- dec
    set.seed(1000)
    sce <- denoisePCA(sce, technical=new.trend, BSPARAM=IrlbaParam())
    set.seed(1000)
    sce <- switch(metadata(sce)$param$visMethod,
                  "TSNE"=runTSNE(sce, use_dimred="PCA",
                                 perplexity=getPerplexity(ncol(sce))),
                  "UMAP"=runUMAP(sce, use_dimred="PCA"),
                  "DiffusionMap"=runDiffusionMap(sce, use_dimred="PCA"),
                  stop("The available visusalisation methods: TSNE, UMAP, DiffusionMap.")
                  )
  }else if(method == "MNN"){
    set.seed(1000)
    ## k, d may need some optimisation
    sce <- do.call(batchelor::fastMNN,
                   c(lapply(rescaled, function(x){x[chosen, ]}),
                     k=20, d=50, BSPARAM=IrlbaParam(deferred=FALSE))
    )
    assay(sce, "counts") <-
      do.call(cbind, lapply(rescaled, function(x){assay(x[chosen, ], "counts")}))
    assay(sce, "logcounts") <-
      do.call(cbind, lapply(rescaled, function(x){assay(x[chosen, ], "logcounts")}))
    metadata(sce)$param <- metadata(sceList[[1]])$param
    message("Features are selected across samples for batch correction: ",
            length(chosen))
    set.seed(1000)
    sce <- switch(metadata(sce)$param$visMethod,
                      "TSNE"=runTSNE(sce, use_dimred="corrected",
                                     perplexity=getPerplexity(ncol(sce))),
                      "UMAP"=runUMAP(sce, use_dimred="corrected"),
                      "DiffusionMap"=runDiffusionMap(sce, use_dimred="corrected"),
                  stop("The available visusalisation methods: TSNE, UMAP, DiffusionMap.")
                  )
  }
  return(sce)
}

scranPreprocessing <- function(sce, param) {
  vars.to.regress <- NULL
  if(identical("CellCycle", param$vars.regress))
    vars.to.regress <- colData(sce)$CellCycle
  
  # Normalization
  clusters <- quickCluster(sce)
  sce = computeSumFactors(sce, cluster = clusters)
  sce = logNormCounts(sce)
  # Variance modelling
  dec <- modelGeneVar(sce, design=vars.to.regress)
  # Extract variable genes for downtream analysis
  top.hvgs <- getTopHVGs(dec, n=2000)
  # PCA. The minimum number of PCAs selected will be 30 and the maximum 50
  sce <- denoisePCA(sce, dec, subset.row=top.hvgs, min.rank=20, max.rank=50)
  sce
}

scranClustering<- function(sce, param) {
   # Clustering with different k values
  resolution_values = seq(from =10, to = 50, by=5)
  clustering_res = factor(clusterCells(sce, use.dimred="PCA", 
                                             BLUSPARAM=SNNGraphParam(k=5, type="jaccard", cluster.fun="louvain")))
  for (k in resolution_values) {
     clustering <- factor(clusterCells(sce, use.dimred="PCA", 
                             BLUSPARAM=SNNGraphParam(k=k, type="jaccard", cluster.fun="louvain")))
     clustering_res = cbind.data.frame(clustering_res, clustering)
  }
  colnames(clustering_res) = paste0("k.", seq(from =5, to = 50, by=5))
  colData(sce) = cbind(colData(sce), clustering_res)
  ident <- factor(clustering_res[,paste0("k.", param$resolution)])
  sce$ident <- ident
  sce
}




scranPosMarkers <- function(sce) {
  markers_any <- findMarkers(sce, sce$ident, pval.type="any", direction="up", lfc=1)
  markersList <- list()
  for(cluster in names(markers_any)){
     markersPerCluster <- data.frame(gene_name = rownames(markers_any[[cluster]]), markers_any[[cluster]], row.names = NULL)
     markersPerCluster <- markersPerCluster[markersPerCluster$FDR<0.05, ]
     markersList[[cluster]] <- markersPerCluster
  }
  markers_any <- bind_rows(markersList, .id="cluster") %>% mutate_if(is.numeric, round, digits=3)
  markers_any
}
  
  
scranDiffGenes <- function(sce) {
  #Creating pseudo-bulk samples
  summed <- aggregateAcrossCells(sce, id=colData(sce)[,c("ident", "Sample_name")])
  #filter out all sample-label combinations with insufficient cells.
  summed.filt <- summed[,summed$ncells >= 10]
  # #create the design matrix
  targets <- colData(sce)[!duplicated(sce$Sample_name),]
  design <-  model.matrix(~factor(Batch)+factor(Condition), data=targets)
  rownames(design) <- targets$Sample_name
  #DE analysis
  de.results <- pseudoBulkDGE(summed.filt, 
                              label=summed.filt$ident,
                              design=~factor(Batch)+factor(Condition),
                              coef=ncol(design),
                              condition=summed.filt$Condition
  )
  is.de <- decideTestsPerLabel(de.results, threshold=0.05)
  summarizeTestsPerLabel(is.de)
  
  diffGenesList <- list()
  for(cluster in names(de.results)){
    diffGenesPerCluster <- data.frame(gene_name = rownames(de.results[[cluster]]), de.results[[cluster]], row.names = NULL)
    diffGenesPerCluster <- diffGenesPerCluster[which(diffGenesPerCluster$FDR<0.05), ]
    diffGenesList[[cluster]] <- diffGenesPerCluster
  }
  diffGenes <- bind_rows(diffGenesList, .id="Cluster")
  diffGenes[,-which(colnames(diffGenes)=="LR")]
}

