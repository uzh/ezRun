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
                               project=param$name,
                               meta.data=cell_info)
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = scData@data),
                     value = TRUE, ignore.case=TRUE)
  perc_mito <- Matrix::colSums(scData@raw.data[mito.genes, ])/Matrix::colSums(scData@raw.data)
  scData <- AddMetaData(object = scData, metadata = perc_mito,
                        col.name = "perc_mito")
  
  perplexityTsne <- switch(param$scProtocol,
                           "smart-Seq2"=10,
                           "10x"=30,
                           stop("Unknown single cell protocol."))
  scData <- FilterCells(object = scData,
                        subset.names = c("nGene", "perc_mito"),
                        low.thresholds = c(param$minGenesPerCell, -Inf), 
                        high.thresholds = c(param$maxGenesPerCell, param$maxMitoFraction))
  scData <- NormalizeData(object=scData, normalization.method="LogNormalize",
                          scale.factor=getSeuratScalingFactor(param$scProtocol))
  scData <- FindVariableGenes(object = scData, do.plot = FALSE,
                              x.low.cutoff=param$x.low.cutoff,
                              x.high.cutoff=param$x.high.cutoff,
                              y.cutoff=param$y.cutoff)
  scData <- ScaleData(object = scData, do.par=TRUE,
                      vars.to.regress = param$vars.to.regress,
                      num.cores=param$cores)
  if(ezIsSpecified(param$pcGenes)){
    indicesMatch <- match(toupper(param$pcGenes), rownames(scData@data))
    if(any(is.na(indicesMatch))){
      warning("The following genes don't exist: ", 
              paste(param$pcGenes[is.na(indicesMatch)], collapse = ","))
    }
    pc.genes <- rownames(scData@data)[which(indicesMatch)]
    metadata(sce)[["pc.genes"]] <- pc.genes
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
                    perplexity=ifelse(length(scData@ident) > 384, 30, 10),
                    num_threads=param$cores)
  metadata(sce)$scData <- scData
  
  return(sce)
}

getSeuratScalingFactor <- function(x=c("10x", "smart-Seq2")){
  x <- match.arg(x)
  ans <- switch(x,
                "smart-Seq2"=1e5,
                "10x"=1e4)
  return(ans)
}
