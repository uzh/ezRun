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
  if(param$scProtocol == "smart-Seq2"){
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
                         dims.use = 1:param$pcs,
                         resolution = param$resolution, print.output = 0, 
                         save.SNN=TRUE, force.recalc=FALSE)
  scData <- RunTSNE(object=scData, reduction.use = "pca",
                    dims.use=1:param$pcs, tsne.method="Rtsne",
                    perplexity=ifelse(length(scData@ident) > 200, 30, 10),
                    num_threads=param$cores)
  scData <- RunUMAP(object=scData, reduction.use = "pca",
                    dims.use=1:param$pcs,
                    n_neighbors=ifelse(length(scData@ident) > 200, 30, 10))
  return(scData)
}

getSeuratScalingFactor <- function(x=c("10X", "smart-Seq2")){
  x <- match.arg(x)
  ans <- switch(x,
                "smart-Seq2"=1e5,
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
