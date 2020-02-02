###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCOneSample <-
  setRefClass("EzAppSCOneSample",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCOneSample
                  name <<- "EzAppSCOneSample"
                  appDefaults <<- rbind(scProtocol=ezFrame(Type="character", DefaultValue="10X", Description="Which single cell protocol?"),
                                        minReadsPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=5e4, 
                                                                Description="Minimal reads per cell of smart-Seq2 for Seurat filtering"),
                                        npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=20,
                                                    Description="The maximal dimensions to use for reduction"),
                                        pcGenes=ezFrame(Type="charVector", 
                                                        DefaultValue="", 
                                                        Description="The genes used in supvervised clustering"),
                                        vars.to.regress=ezFrame(Type="charVector", 
                                                                DefaultValue="cell_cycle,nUMI,perc_mito", 
                                                                Description="Variables to regress out"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.5,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        all2allMarkers=ezFrame(Type="logical", 
                                                               DefaultValue=FALSE, 
                                                               Description="Run all against all cluster comparisons?"),
                                        cellsPercentage=ezFrame(Type="numeric", 
                                                                DefaultValue=0.05, 
                                                                Description="A gene will be kept if it is expressed in at least this percentage of cells"),
                                        nmad=ezFrame(Type="numeric", 
                                                     DefaultValue=3, 
                                                     Description="Median absolute deviation (MAD) from the median value of each metric across all cells"),
                                        species=ezFrame(Type="character", DefaultValue="Human", Description="Organism")
                                        )
                }
              )
  )

ezMethodSCOneSample <- function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  if (identical(param$pcGenes, character(0))) {
     param$pcGenes <- NULL
  }
  sce <- loadSCCountDataset(input, param)
  pvalue_allMarkers <- 0.05
  pvalue_all2allMarkers <- 0.01
  metadata(sce)$output <- output
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
                                    paste(input$getNames(), collapse=", "),
                                    sep=": ")
  
  sce_list <- filterCellsAndGenes(sce, param)
  sce <- sce_list$sce
  sce.unfiltered <- sce_list$sce.unfiltered
  scData <- buildSeuratObject(sce)
  scData <- seuratClusteringV3(scData, param)
  #positive cluster markers
  scData <- posClusterMarkers(scData, pvalue_allMarkers)
  #perform all pairwise comparisons to obtain markers
  if(doEnrichr(param) && param$all2allMarkers) 
    scData = all2all(scData, pvalue_all2allMarkers, param)
  #cell types annotation is only supported for Human and Mouse at the moment
  if(param$species == "Human" | param$species == "Mouse") {
     cells_AUC = cellsLabelsWithAUC(scData, param)
     cellsLabelsWithSingleR(scData, param)
     metadata(sce)$cells_AUC = cells_AUC
  }
  
  scData = saveExternalFiles(scData)
  metadata(sce)$scData = scData
  sce <- findDoublets(sce)
  
  sce_iSEE = as.SingleCellExperiment(scData)
  saveRDS(sce_iSEE, "sce_iSEE.rds")
  saveRDS(sce, "sce.rds")
  saveRDS(sce.unfiltered, "sce.unfiltered.rds")
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCOneSample.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCOneSample.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)

  return("Success")
}

filterCellsAndGenes <- function(sce, param) {
  #browser()
  require(scater)
  require(Matrix)
  library(scDblFinder)
  mito.genes <- grep("^MT-",rowData(sce)$gene_name)
  ribo.genes <- grep("^RP[SL]",rowData(sce)$gene_name)
  
  sce <- addPerCellQC(sce, subsets = list(Mito = mito.genes, Ribo = ribo.genes))
  
  qc.lib <- isOutlier(sce$sum, log=TRUE, nmads=param$nmad, type="lower")
  qc.nexprs <- isOutlier(sce$detected, nmads=param$nmad, log=TRUE, type="lower")
  qc.mito <- isOutlier(sce$subsets_Mito_percent, nmads=param$nmad, type="higher")
  qc.ribo <- isOutlier(sce$subsets_Ribo_percent, nmads=param$nmad, type="higher")
  discard <- qc.lib | qc.nexprs | qc.mito | qc.ribo
  sce$discard <- discard
  sce$qc.lib <- qc.lib
  sce$qc.nexprs <- qc.nexprs
  sce$qc.mito <- qc.mito
  sce$qc.ribo <- qc.ribo
  sce.unfiltered <- sce
  sce <- sce[,!discard]
  
  num.umis <- 1
  num.cells <- param$cellsPercentage*ncol(sce)     #if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  is.expressed <- Matrix::rowSums(counts(sce) >= num.umis) >= num.cells
  sce <- sce[is.expressed,]
  
  rowData(sce.unfiltered)$is.expressed <- is.expressed
  
  return(list(sce.unfiltered=sce.unfiltered, sce = sce))
}

findDoublets <- function(sce) {
  sce <- scDblFinder(sce)
  scData <- metadata(sce)$scData
  scData@meta.data$scDblFinder.score <- colData(sce)$scDblFinder.score
  scData@meta.data$scDblFinder.class <- colData(sce)$scDblFinder.class
  metadata(sce)$scData <- scData
  return(sce)
}

cellsLabelsWithAUC <- function(scData, param) {
  library(AUCell)
  species <- param$species
  if (species == "other")
    return(NULL)
  tissue <- param$tissue
  tissue = unlist(strsplit(tissue, ","))
  all_cell_markers <- read.table("/srv/GT/databases/all_cell_markers.txt", sep = "\t", header = TRUE)
  filtered_cell_markers <- all_cell_markers[all_cell_markers$speciesType == species & all_cell_markers$tissueType %in% tissue, ]
  expressionMatrix <- GetAssayData(scData, slot = "counts")
  geneSets <- createGeneSets(filtered_cell_markers)
  cells_rankings <- AUCell_buildRankings(expressionMatrix, plotStats=FALSE)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  return(cells_AUC)
}

formatGeneSets <- function(name, geneSets) {
  print(geneSets[grep(name, names(geneSets))])
  gs <- unlist(strsplit(unname(unlist(geneSets[grep(name, names(geneSets), fixed=TRUE)])), split=","))
  unique(gsub(" ", "", gs))
} 
  
createGeneSets <- function(filtered_cell_markers) {
 library(GSEABase)
 geneSets <- list()
 cell.names <- as.character(filtered_cell_markers$cellName)
markers <- as.character(filtered_cell_markers$geneSymbol)
geneSets <- mapply(function(cell.names, markers) {markers}, cell.names, markers, SIMPLIFY = FALSE,USE.NAMES = TRUE)
keys = unique(names(geneSets))
geneSets = sapply(keys, formatGeneSets, geneSets = geneSets)
geneSets = sapply(names(geneSets), function(gs.name, geneSets) {GeneSet(geneSets[[gs.name]], setName = gs.name)}, geneSets=geneSets)
geneSets = GeneSetCollection(geneSets)
return(geneSets)
}

cellsLabelsWithSingleR <- function(scData, param) {
library(SingleR) 
singler = CreateSinglerObject(as.matrix(GetAssayData(scData, slot = "counts")), 
                              annot = NULL, project.name = "", min.genes = 0, 
                              technology = "10X", species = param$species, citation = "", 
                              normalize.gene.length = F, variable.genes = "de", fine.tune = T, 
                              do.signatures = F, clusters=Idents(scData), do.main.types = T, 
                              reduce.file.size=T,temp.dir = NULL, numCores = param$cores)

singler$seurat = scData
singler$meta.data$orig.ident = scData@meta.data$Batch
singler$meta.data$xy = scData@reductions$tsne@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = Idents(scData)

#convert the file to use it in the interactive browser
singler.browser = convertSingleR2Browser(singler)
saveRDS(singler.browser, 'singler.browser.rds')
}

saveExternalFiles = function(scData) {
  tr_cnts <- expm1(GetAssayData(scData))
  geneMeans <- rowsum(t(as.matrix(tr_cnts)), group=Idents(scData))
  geneMeans <- sweep(geneMeans, 1, STATS=table(Idents(scData))[rownames(geneMeans)], FUN="/")
  geneMeans <- log1p(t(geneMeans))
  colnames(geneMeans) <- paste("cluster", colnames(geneMeans), sep="_")
  geneMeanPerClusterFn = "gene_means_per_cluster.txt"
  scData@misc$geneMeanPerClusterFn <- geneMeanPerClusterFn
  ezWrite.table(geneMeans, geneMeanPerClusterFn)
  
  geneMeans <- Matrix::rowMeans(tr_cnts)
  geneMeans <- log1p(geneMeans)
  geneMeansFn = "gene_means.txt"
  scData@misc$geneMeansFn <- geneMeansFn
  ezWrite.table(geneMeans, geneMeansFn)
  
  tSNE_data <- as_tibble(scData@reductions$tsne@cell.embeddings,
                         rownames="cells")
  tSNE_data <- dplyr::rename(tSNE_data, X=`tSNE_1`, Y=`tSNE_2`)
  tSNE_data$cluster <- Idents(scData)
  tSNEFn = "tSNE_data.tsv"
  scData@misc$tSNEFn <- tSNEFn
  write_tsv(tSNE_data, path=tSNEFn)
  
  posMarkersFn <- "pos_markers.tsv"
  write_tsv(as_tibble(scData@misc$posMarkers), path=posMarkersFn)
  
  return(scData)
}