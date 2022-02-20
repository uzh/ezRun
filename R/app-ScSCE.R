###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSCE <-
  setRefClass("EzAppScSCE",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodScSCE
        name <<- "EzAppScSCE"
        appDefaults <<- rbind(
          vars.regress = ezFrame(
            Type = "character",
            DefaultValue = "none",
            Description = "Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
          ),
          resolution = ezFrame(
            Type = "numeric",
            DefaultValue = 10,
            Description = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
          ),
          nreads = ezFrame(
            Type = "numeric",
            DefaultValue = Inf,
            Description = "Low quality cells have less than \"nreads\" reads. Only when applying fixed thresholds."
          ),
          ngenes = ezFrame(
            Type = "numeric",
            DefaultValue = Inf,
            Description = "Low quality cells have less than \"ngenes\" genes. Only when applying fixed thresholds."
          ),
          perc_mito = ezFrame(
            Type = "numeric",
            DefaultValue = Inf,
            Description = "Low quality cells have more than \"perc_mito\" percent of mitochondrial genes. Only when applying fixed thresholds."
          ),
          perc_ribo = ezFrame(
            Type = "numeric",
            DefaultValue = Inf,
            Description = "Low quality cells have more than \"perc_ribo\" percent of ribosomal genes. Only when applying fixed thresholds."
          ),
          cellsFraction = ezFrame(
            Type = "numeric",
            DefaultValue = 0.01,
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
          controlSeqs = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "control sequences to add"
          )
        )
      }
    )
  )

ezMethodScSCE <- function(input = NA, output = NA, param = NA,
                                htmlFile = "00index.html") {
  library(HDF5Array)
  library(SingleR)
  library(AUCell)
  library(scran)
  library(bluster)
  library(SingleCellExperiment)
  library(tidyverse)
  require(scDblFinder)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("SC Cluster Report")))
  on.exit(setwd(cwd), add = TRUE)

  sce <- load10xData(input, param)

  pvalue_allMarkers <- 0.05
  pvalue_all2allMarkers <- 0.01

  # Doublets prediction and removal
  library(scDblFinder)
  doubletsInfo <- scDblFinder(counts(sce), returnType = "table")
  doublets <- rownames(doubletsInfo)[doubletsInfo$type == "real" & doubletsInfo$class == "doublet"]
  sce <- sce[,setdiff(colnames(sce),doublets)]

  # Cells and genes filtering
  sce_list <- filterCellsAndGenes(sce, param) # return sce objects filtered and unfiltered to show the QC metrics later in the rmd
  sce <- sce_list$sce
  sce.unfiltered <- sce_list$sce.unfiltered
  rm(sce_list)

  # calculate cellcycle for the filtered sce object
  sce <- addCellCycleToSCE(sce, param$refBuild)
  
  #Pre-processing and clustering
  sce <- scranPreprocessing(sce, param)
  sce <- scranClustering(sce, param)
  
  # Run UMAP and TSNE
  sce <- runUMAP(sce, dimred = "PCA")
  sce <- runTSNE(sce, dimred = "PCA") 
  
  
  #positive cluster markers
  posMarkers <- scranPosMarkers(sce)
  
  singler.results <- NULL
  # cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if (species == "Human" | species == "Mouse") {
    cells_AUC <- cellsLabelsWithAUC(counts(sce), species, param$tissue)
    singler.results <- cellsLabelsWithSingleR(counts(sce), sce$ident, species)
  }
 
  geneMeans <- geneMeansCluster(sce)
  
  ## generate template for manual cluster annotation -----
  ## we only deal with one sample
  stopifnot(length(input$getNames()) == 1)
  clusterInfos <- ezFrame(Sample=input$getNames(), Cluster=levels(sce$ident), ClusterLabel="")
  if (!is.null(singler.results)){
    clusterInfos$SinglerCellType <- singler.results$singler.results.cluster[clusterInfos$Cluster, "pruned.labels"]
  }
  nTopMarkers <- 10
  topMarkers <- posMarkers %>% group_by(cluster) %>%
    slice_max(n = nTopMarkers, order_by = summary.logFC)
  topMarkerString <- sapply(split(topMarkers$gene_name, topMarkers$cluster), paste, collapse=", ")
  clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
  clusterInfoFile <- "clusterInfos.xlsx"
  writexl::write_xlsx(clusterInfos, path=clusterInfoFile)
  
  # Save some results in external files
  dataFiles <- saveExternalFiles(list(pos_markers = posMarkers, gene_means = as_tibble(as.data.frame(geneMeans), rownames = "gene_name")))
  
  metadata(sce)$cells_AUC <- cells_AUC
  metadata(sce)$singler.results <- singler.results
  metadata(sce)$output <- output
  metadata(sce)$param <- param
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
                                    paste(input$getNames(), collapse = ", "),
                                    sep = ": "
  )
  saveHDF5SummarizedExperiment(sce, dir = "sce_h5")
  saveHDF5SummarizedExperiment(sce.unfiltered, dir = "sce.unfiltered_h5")
  
  makeRmdReport(dataFiles = dataFiles, clusterInfoFile=clusterInfoFile, rmdFile = "ScSCE.Rmd", reportTitle = param$name)
  #remove no longer used objects
  rm(sce, sce.unfiltered)
  gc()
  return("Success")
}

