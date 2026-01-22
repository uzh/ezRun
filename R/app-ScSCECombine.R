###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSCECombine <-
  setRefClass(
    "EzAppScSCECombine",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodScSCECombine
        name <<- "EzAppScSCECombine"
        appDefaults <<- rbind(
          batchCorrection = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "Perform batch correction."
          ),
          resolution = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "A numeric value specifying the number of nearest neighbors to consider during graph construction"
          ),
          block = ezFrame(
            Type = "character",
            DefaultValue = "Batch",
            Description = "Blocking factor for each cell"
          )
        )
      }
    )
  )

ezMethodScSCECombine = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  library(HDF5Array)
  library(scran)
  library(bluster)
  library(harmony)
  library(scater)
  library(edgeR)
  library(Seurat)

  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  reportCwd <- getwd()

  sceURLs <- input$getColumn("Static Report")
  filePath <- file.path(
    "/srv/gstore/projects",
    sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "", dirname(sceURLs)),
    "sce_h5"
  )
  filePath_course <- file.path(
    "/srv/GT/analysis/course_sushi/public/projects",
    dirname(sceURLs),
    "sce_h5"
  )

  if (!file.exists(filePath[1])) {
    filePath <- filePath_course
  }

  if (file.exists(filePath[1])) {
    sceList <- lapply(filePath, loadHDF5SummarizedExperiment)
    names(sceList) <- names(sceURLs)
  }

  common_data <- function(sce) {
    common_genes = Reduce(intersect, lapply(sceList, rownames))
    common_colData = Reduce(
      intersect,
      lapply(lapply(sceList, colData), colnames)
    )
    rowData(sce) = NULL
    reducedDims(sce) = NULL
    sce = sce[common_genes, ]
    colData(sce) = colData(sce)[, common_colData]
    colData(sce)[, grep("k.", colnames(colData(sce)))] = NULL #remove previous clustering done on each sample
    return(sce)
  }

  sceList = lapply(sceList, common_data)

  sce = Reduce(SingleCellExperiment::cbind, sceList)

  set.seed(1000)

  #Pre-processing and clustering
  sce <- scranPreprocessing(sce, param)

  # Clustering
  sce_noCorrected <- scranClustering(sce, param)
  sce_noCorrected <- runUMAP(sce_noCorrected, dimred = "PCA", name = "UMAP")
  reducedDim(sce_noCorrected, "UMAP_NOCORRECTED") = reducedDim(
    sce_noCorrected,
    "UMAP"
  )
  sce_noCorrected$ident_noCorrected <- sce_noCorrected$ident
  # If batch correction is needed, perform integration with Harmony using the PCA embeddings
  if (param$batchCorrection) {
    harmony_PCA <- HarmonyMatrix(
      data_mat = reducedDim(sce, "PCA"),
      meta_data = colData(sce),
      vars_use = "Batch",
      do_pca = FALSE
    )
    reducedDim(sce, "PCA_CORRECTED") = harmony_PCA
    sce <- scranClustering(sce, param)
    sce <- runUMAP(sce, dimred = "PCA_CORRECTED", name = "UMAP")
    sce$ident_noCorrected <- sce_noCorrected$ident
    reducedDim(sce, "PCA_NOCORRECTED") = reducedDim(sce_noCorrected, "PCA")
    reducedDim(sce, "UMAP_NOCORRECTED") = reducedDim(sce_noCorrected, "UMAP")
  } else {
    sce = sce_noCorrected
  }

  #positive cluster markers
  posMarkers <- scranPosMarkers(sce)

  # #differentially expressed genes between clusters and conditions (in case of several conditions)
  # diffGenes <- NULL
  # if(length(unique(sce$Condition))>1) {
  #    diffGenes <- scranDiffGenes(sce)
  #    write_tsv(diffGenes, file="differential_genes.tsv")
  # }
  #we do cell type identification only with SingleR since it implemments  block-processing important for large datasets
  singler.results <- NULL
  species <- getSpecies(param$refBuild)
  if (species == "Human" | species == "Mouse") {
    singler.results <- cellsLabelsWithSingleR(counts(sce), sce$ident, species)
  }

  geneMeans <- geneMeansCluster(sce)

  metadata(sce)$singler.results <- singler.results
  metadata(sce)$output <- output
  metadata(sce)$param <- param
  metadata(sce)$param$name <- paste(
    param$name,
    paste(input$getNames(), collapse = ", "),
    sep = ": "
  )

  # Save some results in external files
  dataFiles <- saveExternalFiles(list(
    pos_markers = posMarkers,
    gene_means = as_tibble(as.data.frame(geneMeans), rownames = "gene_name")
  ))
  saveHDF5SummarizedExperiment(sce, dir = "sce_h5")

  makeRmdReport(
    dataFiles = dataFiles,
    rmdFile = "ScSCECombine.Rmd",
    reportTitle = param$name
  )

  #remove no longer used objects
  rm(sce, sceList)
  gc()
  return("Success")
}
