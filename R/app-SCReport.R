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
                                        min_counts=ezFrame(Type="numeric", DefaultValue=5e4, Description="Minimal counts of smart-Seq2 for Seurat filtering"))
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
  ## debug
  saveRDS(sce, file="sce.rds")
  
  
}

seuratProcess <- function(sce){
  require(Seurat)
  
  scData <- CreateSeuratObject(raw.data = cts, min.cells = 5, project = paste0("10X_",project))
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = scData@data), value = TRUE)
  percent.mito <- Matrix::colSums(scData@raw.data[mito.genes, ])/Matrix::colSums(scData@raw.data)
  scData <- AddMetaData(object = scData, metadata = percent.mito, col.name = "percent.mito")
  VlnPlot(object=scData, features.plot=c("nGene", "nUMI", "percent.mito"),
          nCol=3, x.lab.rot=ifelse(nchar(scData@project.name) > 20, TRUE, FALSE))
  scData <- FilterCells(object = scData, subset.names = c("nGene", "percent.mito"), 
                        low.thresholds = c(min_genes, -Inf), high.thresholds = c(max_genes, 0.25))
  scData <- NormalizeData(object = scData, normalization.method = "LogNormalize", 
                          scale.factor = 10000)
  #mm = as.matrix(scData@raw.data)
  #ezWrite.table(mm, "/srv/GT/analysis/p2497/10x_o4048_read_counts.txt")
  #mm = as.matrix(scData@data)
  #ezWrite.table(mm, "/srv/GT/analysis/p2497/10x_o4048_norm_counts.txt")
  
  scData <- FindVariableGenes(object = scData, mean.function = ExpMean, dispersion.function = LogVMR, 
                              x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  scData <- FindVariableGenes(object = scData)
  length(x = scData@var.genes)
  scData <- ScaleData(object = scData, vars.to.regress = c("nUMI", "percent.mito"))
  scData <- RunPCA(object = scData, pc.genes = scData@var.genes, do.print = TRUE, pcs.print = 1:5, 
                   genes.print = 5)
  PrintPCA(object = scData, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
  VizPCA(object = scData, pcs.use = 1:2)
  PCAPlot(object = scData, dim.1 = 1, dim.2 = 2)
  scData <- ProjectPCA(object = scData, do.print = FALSE)
  PCHeatmap(object = scData, pc.use = 1:16, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
  scData <- JackStraw(object = scData, num.replicate = 100)
  JackStrawPlot(object = scData, PCs = 1:20)
  PCElbowPlot(object = scData)
  scData <- FindClusters(object = scData, reduction.type = "pca", dims.use = 1:14, 
                         resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
  PrintFindClustersParams(object = scData)
  set.seed(10);
  scData <- RunTSNE(object = scData, dims.use = 1:14, do.fast = TRUE, perplexity = 30)
  
  
}

