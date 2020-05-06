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
                                        cellsPercentage=ezFrame(Type="numeric", 
                                                                DefaultValue=0.05, 
                                                                Description="A gene will be kept if it is expressed in at least this percentage of cells"),
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
  
  #cell types annotation is only supported for Human and Mouse at the moment
  if(param$species == "Human" | param$species == "Mouse") {
     cells_AUC <- cellsLabelsWithAUC(scData, param)
     singler.results <- cellsLabelsWithSingleR(scData, param)
  }
  
  #Convert scData to Single Cell experiment Object
  sce <- as.SingleCellExperiment(scData)
  metadata(sce)$PCA_stdev <- Reductions(scData, "pca")@stdev   
  metadata(sce)$cells_AUC <- cells_AUC
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
  #browser()
  require(scater)
  require(Matrix)
  
  #Cells filtering
  mito.genes <- grep("^MT-",rowData(sce)$gene_name)
  sce <- addPerCellQC(sce, subsets = list(Mito = mito.genes))
  sce.unfiltered <- filt.lenient(sce)
  sce <- sce[,!sce.unfiltered$discard]
 
  #Genes filtering
  num.umis <- 1
  num.cells <- param$cellsPercentage*ncol(sce)     #if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  is.expressed <- Matrix::rowSums(counts(sce) >= num.umis) >= num.cells
  sce <- sce[is.expressed,]
  rowData(sce.unfiltered)$is.expressed <- is.expressed
  
  return(list(sce.unfiltered=sce.unfiltered, sce = sce))
}

add_meta <- function(ds){
  ds$total_features <- ds$detected
  ds$log10_total_features <- log10(ds$detected)
  ds$total_counts <- ds$sum
  ds$log10_total_counts <- log10(ds$sum+1)
  ds$featcount_ratio <- ds$log10_total_counts/ds$log10_total_features
  ds$featcount_dist <- getFeatCountDist(ds)
  ds$pct_counts_top_50_features <- ds$percent_top_50
  ds
}

getFeatCountDist <- function(df, do.plot=FALSE, linear=TRUE){
  if(is(df,"SingleCellExperiment")) df <- as.data.frame(colData(df))
  if(linear){
    mod <- lm(df$log10_total_features~df$log10_total_counts)
  }else{
    mod <- loess(df$log10_total_features~df$log10_total_counts)
  }
  pred <- predict(mod, newdata=data.frame(log10_total_counts=df$log10_total_counts))
  df$diff <- df$log10_total_features - pred
  df$diff
}

filt.lenient <- function(x){  
  if(!("featcount_dist" %in% colnames(colData(x)))) x <- add_meta(x)
  filters <- c( "log10_total_counts:both:5",
                "log10_total_features:both:5",
                "pct_counts_top_50_features:both:5",
                "featcount_dist:both:5")
  out <- lapply(strsplit(filters,":"), FUN=function(f) {
    #  browser()
    which(isOutlier(x[[f[1]]], log=FALSE,nmads=as.numeric(f[3]),type=f[2]))
    # x[[paste0("qc.", f[1])]] <- isOutlier(x[[f[1]]], log=FALSE,
    #  nmads=as.numeric(f[3]),
    #  type=f[2])
    # which(x[[paste0("qc.", f[1])]])
  })
  x$qc.total_counts <- FALSE
  x$qc.total_counts[out[[1]]] <- TRUE
  x$qc.total_features <- FALSE
  x$qc.total_features[out[[2]]] <- TRUE
  x$qc.pct_counts_top_50_features <- FALSE
  x$qc.pct_counts_top_50_features[out[[3]]] <- TRUE
  x$qc.featcount_dist <- FALSE
  x$qc.featcount_dist[out[[4]]] <- TRUE
  
  mtout <- isOutlier(x$subsets_Mito_percent, nmads=3, type="lower" ) | 
    (isOutlier(x$subsets_Mito_percent, nmads=3, type="higher" ) & x$subsets_Mito_percent > 0.08)
  out <- c(out, list(mt=which(mtout)))
  out <- table(unlist(out))
  out <- as.numeric(names(out)[which(out>=2)])
  x$discard <- FALSE
  x$discard[out] <- TRUE
  return(x)
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
if(grepl("Homo_sapiens", param$refBuild)){
    hpca.se <- HumanPrimaryCellAtlasData()
    bp.se <- BlueprintEncodeData()
    singler.results.single <- SingleR(test = GetAssayData(scData), ref = list(BP=bp.se, HPCA=hpca.se), 
                               labels = list(bp.se$label.main, hpca.se$label.main), method="single", de.method = "wilcox")
    singler.results.cluster <- SingleR(test = GetAssayData(scData), ref = list(BP=bp.se, HPCA=hpca.se), 
                                      labels = list(bp.se$label.main, hpca.se$label.main), method="cluster", clusters=Idents(scData), de.method = "wilcox")
  }else {
    hpca.se <- MouseRNAseqData()
    singler.results.single <- SingleR(test = GetAssayData(scData), ref = hpca.se, labels = hpca.se$label.main)
    singler.results.cluster <- SingleR(test = GetAssayData(scData), ref = hpca.se, labels = hpca.se$label.main, method="cluster", clusters=Idents(scData))
  }
  
return(list(singler.results.single=singler.results.single, singler.results.cluster=singler.results.cluster))
}

