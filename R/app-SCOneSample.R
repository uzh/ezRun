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
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSCOneSample
        name <<- "EzAppSCOneSample"
        appDefaults <<- rbind(
          scProtocol = ezFrame(Type = "character", DefaultValue = "10X", Description = "Which single cell protocol?"),
          minReadsPerCell = ezFrame(
            Type = "numeric",
            DefaultValue = 5e4,
            Description = "Minimal reads per cell of smart-Seq2 for Seurat filtering"
          ),
          npcs = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "The maximal dimensions to use for reduction"
          ),
          pcGenes = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "The genes used in supvervised clustering"
          ),
          SCT.regress = ezFrame(
            Type = "character",
            DefaultValue = "none",
            Description = "Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
          ),
          DE.method = ezFrame(
            Type = "charVector",
            DefaultValue = "wilcoxon",
            Description = "Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."
          ),
          resolution = ezFrame(
            Type = "numeric",
            DefaultValue = 0.5,
            Description = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
          ),
          all2allMarkers = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Run all against all cluster comparisons?"
          ),
          cellsFraction = ezFrame(
            Type = "numeric",
            DefaultValue = 0.05,
            Description = "A gene will be kept if it is expressed in at least this percentage of cells"
          ),
          nUMIs = ezFrame(
            Type = "numeric",
            DefaultValue = 1,
            Description = "A gene will be kept if it has at least nUMIs in the fraction of cells specified before"
          ),
          nmad = ezFrame(
            Type = "numeric",
            DefaultValue = 3,
            Description = "Median absolute deviation (MAD) from the median value of each metric across all cells"
          ),
          filterByExpression = ezFrame(
            Type = "character", DefaultValue = FALSE,
            Description = "Keep cells according to specific gene expression. i.e. Set > 1 | Pkn3 > 1"
          ),
          controlSeqs = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "control sequences to add"
          )
        )
      }
    )
  )

ezMethodSCOneSample <- function(input = NA, output = NA, param = NA,
                                htmlFile = "00index.html") {
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)

  sce <- load10xData(input, param)

  pvalue_allMarkers <- 0.05
  pvalue_all2allMarkers <- 0.01

  # Doublets prediction and removal
  library(scDblFinder)
  sce <- scDblFinder(sce)
  sce <- sce[, which(sce$scDblFinder.class != "doublet")]
  # scData@meta.data$scDblFinder.score <- colData(sce)$scDblFinder.score
  # scData@meta.data$scDblFinder.class <- colData(sce)$scDblFinder.class

  # Cells and genes filtering
  sce_list <- filterCellsAndGenes(sce, param) # return sce objects filtered and unfiltered to show the QC metrics later in the rmd
  sce <- sce_list$sce
  sce.unfiltered <- sce_list$sce.unfiltered
  rm(sce_list)

  # calculate cellcycle for the filtered sce object
  sce <- addCellCycleToSce(sce, param$refBuild)

  #before converting to a seurat object, replace all - by . to avoid future problems when subsetting the object
  rownames(sce) <- gsub("-", ".", rownames(sce))
  rowData(sce)$Symbol <- gsub("-", ".", rowData(sce)$Symbol)
  
  scData <- buildSeuratObject(sce) # the Seurat object is built from the filtered sce object
  if (param$filterByExpression != "") {
    expression <- param$filterByExpression
    myCommand <- paste("subset(scData,", expression, ")")
    scData <- eval(parse(text = myCommand))
  }
  scData <- seuratClusteringV3(scData, param)

  # positive cluster markers
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  # if all2allmarkers are not calculated it will remain as NULL
  all2allMarkers <- NULL
  # perform all pairwise comparisons to obtain markers
  if (doEnrichr(param) && param$all2allMarkers) {
    all2allMarkers <- all2all(scData, pvalue_all2allMarkers, param)
  }

  cells_AUC <- NULL
  singler.results <- NULL
  # cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if (species == "Human" | species == "Mouse") {
    cells_AUC <- cellsLabelsWithAUC(scData, species, param$tissue)
    singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "counts"), Idents(scData), species)
  }

  # Convert scData to Single Cell experiment Object
  library(scanalysis)
  #TODO: remove unnecesary dietseurat call when the bug in Seurat is fixed
  scData_diet = DietSeurat(scData, dimreducs = c("pca", "tsne", "umap"))
  sce <- scData_diet %>% seurat_to_sce(default_assay = "SCT")
  metadata(sce)$PCA_stdev <- Reductions(scData_diet, "pca")@stdev
  metadata(sce)$cells_AUC <- cells_AUC
  metadata(sce)$singler.results <- singler.results
  metadata(sce)$output <- output
  metadata(sce)$param <- param
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
    paste(input$getNames(), collapse = ", "),
    sep = ": "
  )


  geneMeans <- geneMeansCluster(sce)
  # Save some results in external files
  dataFiles <- saveExternalFiles(list(pos_markers = posMarkers, all2allMarkers = all2allMarkers, gene_means = as_tibble(as.data.frame(geneMeans), rownames = "gene_name")))
  # rowData(sce) = rowData(sce)[, c("gene_id", "biotypes", "description")]

  library(HDF5Array)
  saveHDF5SummarizedExperiment(sce, dir = "sce_h5")
  saveHDF5SummarizedExperiment(sce.unfiltered, dir = "sce.unfiltered_h5")

  makeRmdReport(dataFiles = dataFiles, rmdFile = "SCOneSample.Rmd", reportTitle = metadata(sce)$param$name)
  return("Success")
}

filterCellsAndGenes <- function(sce, param) {
  library(scater)
  library(Matrix)

  # Cells filtering
  mito.genes <- grep("^MT.", rowData(sce)$Symbol, ignore.case = TRUE)
  ribo.genes <- grep("^RPS|^RPL", rownames(sce), ignore.case = TRUE)
  
  sce <- addPerCellQC(sce, subsets = list(Mito = mito.genes, Ribo = ribo.genes))

  if (param$nreads == "") {
    qc.lib <- isOutlier(sce$sum, log = TRUE, nmads = param$nmad, type = "lower")
  } else {
    qc.lib <- sce$sum < as.double(param$nreads)
  }
  if (param$ngenes == "") {
    qc.nexprs <- isOutlier(sce$detected, nmads = param$nmad, log = TRUE, type = "lower")
  } else {
    qc.nexprs <- sce$detected < as.double(param$ngenes)
  }
  if (param$perc_mito == "") {
    qc.mito <- isOutlier(sce$subsets_Mito_percent, nmads = param$nmad, type = "higher")
  } else {
    qc.mito <- sce$subsets_Mito_percent > as.double(param$perc_mito)
  }
  
  if (param$perc_ribo == "") {
    qc.ribo <- isOutlier(sce$subsets_Ribo_percent, nmads = param$nmad, type = "higher")
  } else {
    qc.ribo <- sce$subsets_Ribo_percent > as.double(param$perc_ribo)
  }
  
  discard <- qc.lib | qc.nexprs | qc.mito | qc.ribo
  sce$discard <- discard
  sce$qc.lib <- qc.lib
  sce$qc.nexprs <- qc.nexprs
  sce$qc.mito <- qc.mito
  sce$qc.ribo <- qc.ribo
  sce.unfiltered <- sce
  sce <- sce[, !discard]

  # Genes filtering
  num.cells <- param$cellsFraction * ncol(sce) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  is.expressed <- Matrix::rowSums(counts(sce) >= param$nUMIs) >= num.cells
  sce <- sce[is.expressed, ]
  rowData(sce.unfiltered)$is.expressed <- is.expressed

  return(list(sce.unfiltered = sce.unfiltered, sce = sce))
}
