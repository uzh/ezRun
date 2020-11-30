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
                                        SCT.regress=ezFrame(Type="character", 
                                                           DefaultValue="none", 
                                                           Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."),
                                        DE.method=ezFrame(Type="charVector", 
                                                          DefaultValue="wilcoxon", 
                                                          Description="Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.5,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        all2allMarkers=ezFrame(Type="logical", 
                                                               DefaultValue=FALSE, 
                                                               Description="Run all against all cluster comparisons?"),
                                        cellsFraction=ezFrame(Type="numeric", 
                                                                DefaultValue=0.05, 
                                                                Description="A gene will be kept if it is expressed in at least this percentage of cells"),
                                        nUMIs=ezFrame(Type="numeric", 
                                                      DefaultValue=1, 
                                                      Description='A gene will be kept if it has at least nUMIs in the fraction of cells specified before'),
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

  sce <- loadSCCountDataset(input, param)
  pvalue_allMarkers <- 0.05
  pvalue_all2allMarkers <- 0.01
  
  #Doublets prediction and removal
  require(scDblFinder)
  sce <- scDblFinder(sce)
  sce <- sce[,which(sce$scDblFinder.class!="doublet")]
  # scData@meta.data$scDblFinder.score <- colData(sce)$scDblFinder.score
  # scData@meta.data$scDblFinder.class <- colData(sce)$scDblFinder.class
  
  #Cells and genes filtering
  sce_list <- filterCellsAndGenes(sce, param)  #return sce objects filtered and unfiltered to show the QC metrics later in the rmd 
  sce <- sce_list$sce
  sce.unfiltered <- sce_list$sce.unfiltered
  
  scData <- buildSeuratObject(sce)   # the Seurat object is built from the filtered sce object
  scData <- seuratClusteringV3(scData, param)
  
  #positive cluster markers
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  #if all2allmarkers are not calculated it will remain as NULL
  all2allMarkers <- NULL
  #perform all pairwise comparisons to obtain markers
  if(doEnrichr(param) && param$all2allMarkers) 
    all2allMarkers = all2all(scData, pvalue_all2allMarkers, param)
  
  cells_AUC <- NULL
  singler.results <- NULL
  #cell types annotation is only supported for Human and Mouse at the moment
  if(param$species == "Human" | param$species == "Mouse") {
     cells_AUC <- cellsLabelsWithAUC(scData, param)
     singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "counts"), Idents(scData), param)
  }
  
  #Convert scData to Single Cell experiment Object
  # Use the SCT logcounts for visualization instead of the RNA logcounts.
  # TODO: save all the assays (RNA, SCT) in the sce object using the package keshavmot2/scanalysis. The function from Seurat doesn't save everything.
  DefaultAssay(scData) <- "SCT" 
  sce <- as.SingleCellExperiment(scData)
  metadata(sce)$PCA_stdev <- Reductions(scData, "pca")@stdev   
  metadata(sce)$cells_AUC <- cells_AUC
  singler.results.cluster <- singler.results$singler.results.cluster
  singler.single.labels <-  singler.results$singler.results.single
  sce$singler.cluster.labels = singler.results.cluster[as.character(sce$ident), "labels"]
  sce$singler.single.labels <- singler.single.labels$labels
  metadata(sce)$singler.results <- singler.results
  metadata(sce)$output <- output
  metadata(sce)$param <- param
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
                                    paste(input$getNames(), collapse=", "),
                                    sep=": ")
  
  
  #Save some results in external files 
  saveExternalFiles(sce, list(pos_markers=posMarkers, all2allMarkers=all2allMarkers))
 # rowData(sce) = rowData(sce)[, c("gene_id", "biotypes", "description")]
  
  library(HDF5Array)
  saveHDF5SummarizedExperiment(sce, dir="sce_h5")
  
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCOneSample.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  while (dev.cur()>1) dev.off()
  rmarkdown::render(input="SCOneSample.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)

  return("Success")
}

filterCellsAndGenes <- function(sce, param) {
  require(scater)
  require(Matrix)
  
  #Cells filtering
  mito.genes <- grep("^MT-",rowData(sce)$gene_name, ignore.case = TRUE)
  
  sce <- addPerCellQC(sce, subsets = list(Mito = mito.genes))
  
  if(!(param$nreads== "") & !(param$ngenes== "") & !(param$perc_mito== "")) {
    qc.lib <- sce$sum < param$nreads
    qc.nexprs <- sce$detected < param$ngenes
    qc.mito <- sce$subsets_Mito_percent > param$perc_mito
  } else {
    qc.lib <- isOutlier(sce$sum, log=TRUE, nmads=param$nmad, type="lower")
    qc.nexprs <- isOutlier(sce$detected, nmads=param$nmad, log=TRUE, type="lower")
    qc.mito <- isOutlier(sce$subsets_Mito_percent, nmads=param$nmad, type="higher")
  }
  discard <- qc.lib | qc.nexprs | qc.mito
  sce$discard <- discard
  sce$qc.lib <- qc.lib
  sce$qc.nexprs <- qc.nexprs
  sce$qc.mito <- qc.mito
  sce.unfiltered <- sce
  sce <- sce[,!discard]
 
  #Genes filtering
  num.cells <- param$cellsFraction*ncol(sce)     #if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  is.expressed <- Matrix::rowSums(counts(sce) >= param$nUMIs) >= num.cells
  sce <- sce[is.expressed,]
  rowData(sce.unfiltered)$is.expressed <- is.expressed
  
  return(list(sce.unfiltered=sce.unfiltered, sce = sce))
}

cellsLabelsWithAUC <- function(scData, param) {
  library(AUCell)
  species <- param$species
  if (species == "other")
    return(NULL)
  tissue <- param$tissue
  tissue = unlist(strsplit(tissue, ","))
  all_cell_markers <- read.table("/srv/GT/databases/scGeneSets/all_cell_markers.txt", sep = "\t", header = TRUE)
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

cellsLabelsWithSingleR <- function(counts, current_clusters, param) {
library(SingleR)
if(grepl("Homo_sapiens", param$refBuild)){
    reference <- HumanPrimaryCellAtlasData()
    singler.results.single <- SingleR(test = counts, ref = reference, 
                               labels = reference$label.main, method="single", de.method = "wilcox")
    singler.results.cluster <- SingleR(test = counts, ref = reference, 
                                      labels = reference$label.main, method="cluster", clusters=current_clusters, de.method = "wilcox")
  }else {
    reference <- MouseRNAseqData()
    singler.results.single <- SingleR(test = counts, ref = reference, labels = reference$label.main)
    singler.results.cluster <- SingleR(test = counts, ref = reference, labels = reference$label.main, method="cluster", clusters=current_clusters)
  }
  
return(list(singler.results.single=singler.results.single, singler.results.cluster=singler.results.cluster))
}

