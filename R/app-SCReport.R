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
                  appDefaults <<- rbind(minGenesPerCell=ezFrame(Type="numeric", DefaultValue=500, Description="Minimal number of genes per cell for Seurat filtering"),
                                        maxGenesPerCell=ezFrame(Type="numeric", DefaultValue=3000, Description="Maximal number of genes per cell for Seurat filtering"),
                                        minReadsPerCell=ezFrame(Type="numeric", DefaultValue=5e4, Description="Minimal reads per cell of smart-Seq2 for Seurat filtering"),
                                        pcs=ezFrame(Type="numeric", DefaultValue=10, Description="The maximal dimensions to use for reduction"),
                                        pcGenes=ezFrame(Type="charVector", DefaultValue="", Description="The genes used in supvervised clustering"),
                                        x.low.cutoff=ezFrame(Type="numeric", DefaultValue=0.1, Description="Bottom cutoff on x-axis for identifying variable genes"),
                                        x.high.cutoff=ezFrame(Type="numeric", DefaultValue=8, Description="Top cutoff on x-axis for identifying variable genes"),
                                        y.cutoff=ezFrame(Type="numeric", DefaultValue=1, Description="Bottom cutoff on y-axis for identifying variable genes"),
                                        resolution=ezFrame(Type="numeric", DefaultValue=0.8, Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        markersToCheck=ezFrame(Type="charList", DefaultValue="", Description="The markers to check"),
                                        runPseudoTime=ezFrame(Type="logical", DefaultValue=FALSE, Description="Run PseudoTime for single cell data?"))
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
  sce <- seuratPreProcess(sce)
  
  ## debug
  saveRDS(sce, file="sce.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCReport.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCReport.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  return("Success")
}

seuratPreProcess <- function(sce){
  ## parameters to tune: param$minReadsPerCell; param$maxGenesPerCell; param$minGenesPerCell
  ##                     param$pcs, param$pcGenes
  
  require(Seurat)
  require(scater)
  param <- metadata(sce)$param
  
  rownames(sce) <- uniquifyFeatureNames(ID=rowData(sce)$gene_id,
                                        names=rowData(sce)$gene_name)
  countsSeurat <- assays(sce)$counts
  if(param$scProtocol == "smart-Seq2"){
    countsSeurat <- countsSeurat[ ,Matrix::colSums(countsSeurat) > param$minReadsPerCell]
  }
  
  isMito <- toupper(as.character(seqnames(rowRanges(sce)))) %in% toupper(c("chrM", "MT"))
  cell_info <- data.frame(row.names=colnames(countsSeurat),
                          lib_size=Matrix::colSums(countsSeurat)/1e6,
                          expr_genes = Matrix::colSums(countsSeurat >= param$minReadsPerGene),
                          perc_mito=Matrix::colSums(countsSeurat[isMito, , drop=FALSE])*100 / Matrix::colSums(countsSeurat))
  cell_info_seurat <- cell_info
  cell_info_seurat$perc_mito <- cell_info_seurat$perc_mito/100
  
  scData <- CreateSeuratObject(raw.data = countsSeurat, min.cells = 5,
                               project = param$name,
                               meta.data=cell_info_seurat)
  if(param$scProtocol == "smart-Seq2"){
    scalingFactorSeurat <- 1e5
  }else if(param$scProtocol == "10x"){
    scalingFactorSeurat <- 1e4
  }
  
  scData <- FilterCells(object = scData,
                        subset.names = c("nGene", "perc_mito"),
                        low.thresholds = c(param$minGenesPerCell, -Inf), 
                        high.thresholds = c(param$maxGenesPerCell, 0.25))
  scData <- NormalizeData(object = scData, normalization.method = "LogNormalize",
                          scale.factor = scalingFactorSeurat)
  scData <- FindVariableGenes(object = scData, do.plot = FALSE,
                              x.low.cutoff=param$x.low.cutoff,
                              x.high.cutoff=param$x.high.cutoff,
                              y.cutoff=param$y.cutoff)
  scData <- ScaleData(object = scData, do.par=TRUE,
                      vars.to.regress = c("nUMI", "perc_mito"),
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
  scData <- RunPCA(object = scData, pc.genes = pc.genes,
                   do.print = FALSE, pcs.print = 1:5,
                   genes.print = 5)
  scData <- ProjectPCA(object = scData, do.print = FALSE)
  scData <- JackStraw(object=scData, num.replicate=100, display.progress=FALSE,
                      do.par=TRUE, num.cores=param$cores)
  
  scData <- FindClusters(object=scData, reduction.type="pca",
                         dims.use = 1:param$pcs,
                         resolution = param$resolution, print.output = 0, 
                         save.SNN=TRUE, force.recalc=TRUE)
  set.seed(10)
  scData <- RunTSNE(object = scData, dims.use = 1:param$pcs, do.fast = TRUE, 
                    perplexity = 30)
  metadata(sce)$scData <- scData
  
  return(sce)
}

