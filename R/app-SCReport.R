###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCReport <-
  setRefClass("EzAppSCReport",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCReport
                  name <<- "EzAppSCReport"
                  appDefaults <<- rbind(min_genes=ezFrame(Type="numeric", DefaultValue=500, Description="Minimal number of genes for Seurat filtering"),
                                        max_genes=ezFrame(Type="numeric", DefaultValue=3000, Description="Minimal number of genes for Seurat filtering"),
                                        min_counts=ezFrame(Type="numeric", DefaultValue=5e4, Description="Minimal counts of smart-Seq2 for Seurat filtering"),
                                        pcs=ezFrame(Type="numeric", DefaultValue=10, Description="The maximal dimensions to use for reduction"))
                }
              )
  )

ezMethodSCReport = function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  param$scProtocol <- ifelse("STARLog" %in% input$colNames, "smart-Seq2", "10x")
  
  sce <- loadSCCountDataset(input, param)
  scData <- seuratPreProcess(sce)
  metadata(sce)$scData <- scData
  
  ## debug
  saveRDS(sce, file="sce.rds")
  
}

seuratPreProcess <- function(sce){
  ## parameters to tune: param$min_counts; param$max_genes; param$min_genes
  ##                     param$pcs
  
  require(Seurat)
  require(scater)
  sceSeurat <- sce
  param <- metadata(sceSeurat)$param
  
  rownames(sceSeurat) <- uniquifyFeatureNames(ID=rowData(sceSeurat)$gene_id,
                                              names=rowData(sceSeurat)$gene_name)
  countsSeurat <- assays(sceSeurat)$counts
  if(param$scProtocol == "smart-Seq2"){
    countsSeurat <- countsSeurat[ ,Matrix::colSums(countsSeurat) > param$min_counts]
  }
  
  scData <- CreateSeuratObject(raw.data = countsSeurat, min.cells = 5,
                               project = param$name)
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = scData@data),
                     value = TRUE, ignore.case = TRUE)
  percent.mito <- Matrix::colSums(scData@raw.data[mito.genes, ]) / Matrix::colSums(scData@raw.data)
  scData <- AddMetaData(object = scData, metadata = percent.mito,
                        col.name = "percent.mito")
  if(param$scProtocol == "smart-Seq2"){
    scalingFactorSeurat <- 1e5
  }else if(param$scProtocol == "10x"){
    scalingFactorSeurat <- 1e4
  }
  
  scData <- FilterCells(object = scData,
                        subset.names = c("nGene", "percent.mito"),
                        low.thresholds = c(param$min_genes, -Inf), 
                        high.thresholds = c(param$max_genes, 0.25))
  scData <- NormalizeData(object = scData, normalization.method = "LogNormalize",
                          scale.factor = scalingFactorSeurat)
  scData <- FindVariableGenes(object = scData, do.plot = FALSE)
  scData <- ScaleData(object = scData,
                      vars.to.regress = c("nUMI", "percent.mito"))
  scData <- RunPCA(object = scData, pc.genes = scData@var.genes,
                   do.print = FALSE, pcs.print = 1:5,
                   genes.print = 5)
  scData <- ProjectPCA(object = scData, do.print = FALSE)
  scData <- JackStraw(object=scData, num.replicate=100, display.progress=FALSE,
                      do.par=TRUE, num.cores=param$cores)
  
  scData <- FindClusters(object=scData, reduction.type="pca",
                         dims.use = 1:param$pcs,
                         resolution = 0.6, print.output = 0, save.SNN=TRUE,
                         force.recalc=TRUE)
  set.seed(10)
  scData <- RunTSNE(object = scData, dims.use = 1:param$pcs, do.fast = TRUE, 
                    perplexity = 30)
  return(scData)
}

