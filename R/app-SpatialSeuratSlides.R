###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSpatialSeuratSlides <-
  setRefClass(
    "EzAppSpatialSeuratSlides",
    contains = "EzApp",
    methods = list(
      ## Seurat defaults quoted here were checked against Seurat 5.5.1 formals() (R 4.6)
      ## and are identical in 5.1.0 (Dev/R/4.4.2, which SpatialSeuratSlidesApp.rb loads).
      ## Seurat v5 + SCTransform v2 unconditional; Seurat CCA anchor integration only when
      ## batchCorrection is true and integrationMethod is CCA (both appDefaults, not on the form).
      citation = function(param = list()) {
        c(
          "Hao, Y. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology 42, 293-304 (2024). https://doi.org/10.1038/s41587-023-01767-y",
          "Hafemeister, C. & Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biology 20, 296 (2019). https://doi.org/10.1186/s13059-019-1874-1",
          "Choudhary, S. & Satija, R. Comparison and evaluation of statistical error models for scRNA-seq. Genome Biology 23, 27 (2022). https://doi.org/10.1186/s13059-021-02584-9",
          if (isTRUE(as.logical(param$batchCorrection)) && identical(param$integrationMethod, "CCA")) "Stuart, T. et al. Comprehensive integration of single-cell data. Cell 177, 1888-1902 (2019). https://doi.org/10.1016/j.cell.2019.05.031"
        )
      },
      methods_facts = function(param = list()) {
        batchCorrection <- as.logical(param$batchCorrection)
        c(
          ## ezMethodSpatialSeuratSlides (app-SpatialSeuratSlides.R:105-117)
          "Each slide's Seurat object from its single-slide SpatialSeurat analysis (a SpatialSeurat version from before September 2025) was loaded and its previous SCT-based cluster assignments were removed.",
          ## seuratNormalizeSampleList SCTransform without seed.use; no set.seed in this app (seuratUtils.R:308-335)
          "Each slide was re-normalized separately with SCTransform at Seurat defaults (vst.flavor = v2, 3000 variable features, Seurat's default seed 1448145; no global random seed was set).",
          ## getSeuratVarsToRegress (seuratUtils.R:781-793)
          if (isTRUE(as.logical(param$SCT.regress.CellCycle))) "The cell-cycle score difference (S minus G2M) was regressed out in SCTransform.",
          ## cellClustNoCorrection (seuratUtils.R:354-378); kept as umap_noCorrected / ident_noCorrected (app-SpatialSeuratSlides.R:131-133)
          "An uncorrected view was computed by merging the slides, using the union of the per-slide SCTransform variable genes, and running the PCA, UMAP and clustering steps below.",
          if (isFALSE(batchCorrection)) "Batch correction was not applied (batchCorrection false), so the reported UMAP and clusters are the uncorrected ones.",
          ## cellClustWithCorrection CCA branch (seuratUtils.R:381-446); integrationMethod appDefault "CCA", nfeatures appDefault 3000, neither declared in the Ruby app
          if (isTRUE(batchCorrection) && identical(param$integrationMethod, "CCA")) paste0("Slides were integrated with Seurat CCA anchors (integrationMethod, an app default not set on the form): ", param$nfeatures %||% 3000, " integration features (SelectIntegrationFeatures), PrepSCTIntegration, FindIntegrationAnchors with SCT normalization and dimensions 1 to npcs (Seurat defaults k.anchor 5, k.filter 200, k.score 30), and IntegrateData with dimensions 1 to npcs (k.weight 100); the reported UMAP and clusters come from the integrated assay, and the uncorrected UMAP and clusters were kept alongside."),
          ## seuratStandardWorkflow (seuratUtils.R:170-180)
          "PCA computed 50 components (Seurat default) and the first npcs were used for the neighbour graph and for UMAP (uwot, cosine metric, 30 neighbours, seed 42; RunUMAP defaults); t-SNE was not computed and pcGenes was not used.",
          ## seuratStandardWorkflow (seuratUtils.R:181-235)
          "The shared-nearest-neighbour graph used k = 20 (FindNeighbors default); clusters were found with the Louvain algorithm (FindClusters algorithm 1) at resolutions 0.2, 0.4, 0.6, 0.8 and 1.0 plus the resolution parameter, and the clustering at the resolution parameter is the one reported.",
          ## PrepSCTFindMarkers + posClusterMarkers (app-SpatialSeuratSlides.R:134-137; seuratUtils.R:516-548)
          "Before marker detection the per-slide SCT models were reconciled with PrepSCTFindMarkers, and cluster markers were found with FindAllMarkers on the SCT data using the test in DE.method, positive markers only.",
          ## posClusterMarkers passes neither min.pct nor logfc.threshold; latent.vars only for SCOneSample/SCReportMerging (seuratUtils.R:517-531)
          "The min.pct and logfc.threshold values on the form were not passed to FindAllMarkers, so the Seurat defaults applied (min.pct 0.01, logfc.threshold 0.1), and DE.regress was not applied (no latent variables, also for LR).",
          ## pvalue_allMarkers <- 0.05 (app-SpatialSeuratSlides.R:119) used as return.thresh; no p_val_adj filter (seuratUtils.R:530-547)
          "Markers were reported at an unadjusted p-value below 0.05 (hardcoded), with no filter on the adjusted p-value.",
          ## app-SpatialSeuratSlides.R:139-161
          "Spatially variable genes were not recomputed; genes ranked by both markvariogram and Moran's I in a single-slide analysis (non-missing MeanRank) were reused to flag cluster markers as spatial markers."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSpatialSeuratSlides
        name <<- "EzAppSpatialSeuratSlides"
        appDefaults <<- rbind(
          npcs = ezFrame(
            Type = "numeric",
            DefaultValue = 30,
            Description = "The maximal dimensions to use for reduction"
          ),
          nfeatures = ezFrame(
            Type = "numeric",
            DefaultValue = 3000,
            Description = "number of variable genes for SCT"
          ),
          pcGenes = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "The genes used in supvervised clustering"
          ),
          resolution = ezFrame(
            Type = "numeric",
            DefaultValue = 0.6,
            Description = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
          ),
          SCT.regress.CellCycle = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
          ),
          batchCorrection = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "Perform batch correction."
          ),
          integrationMethod = ezFrame(
            Type = "character",
            DefaultValue = "CCA",
            Description = "Choose integration method in Seurat"
          ),
          DE.method = ezFrame(
            Type = "charVector",
            DefaultValue = "wilcox",
            Description = "Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"
          ),
          DE.regress = ezFrame(
            Type = "charVector",
            DefaultValue = "Batch",
            Description = "Variables to regress out if the test LR is chosen"
          ),
          min.pct = ezFrame(
            Type = "numeric",
            DefaultValue = 0.1,
            Description = "Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations."
          ),
          logfc.threshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.25,
            Description = "Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
          ),
          pt.size.factor = ezFrame(
            Type = "numeric",
            DefaultValue = 1,
            Description = "pt.size.factor for spatial plots"
          ),
          maxSamplesSupported = ezFrame(
            Type = "numeric",
            DefaultValue = 5,
            Description = "Maximum number of samples to compare"
          )
        )
      }
    )
  )

ezMethodSpatialSeuratSlides = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  library(Seurat)
  library(rlist)
  library(HDF5Array)
  library(SummarizedExperiment)
  library(SingleCellExperiment)

  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  reportCwd <- getwd()

  scDataURLs <- input$getColumn("Static Report")
  filePath <- file.path(
    "/srv/gstore/projects",
    sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "", dirname(scDataURLs)),
    "scData.rds"
  )

  scDataList <- lapply(filePath, readRDS)
  names(scDataList) <- names(scDataURLs)
  scDataList <- lapply(scDataList, function(scData) {
    scData@meta.data[, grep("SCT", colnames(scData@meta.data))] = NULL #remove previous clustering done on SCT assay
    scData
  })

  pvalue_allMarkers <- 0.05

  scData_noCorrected <- cellClustNoCorrection(scDataList, param)
  scData = scData_noCorrected
  if (param$batchCorrection) {
    scData_corrected = cellClustWithCorrection(scDataList, param)
    #in order to compute the markers we switch again to the original assay
    varFeatures <- VariableFeatures(scData_corrected)
    DefaultAssay(scData_corrected) <- "SCT"
    scData <- scData_corrected
    VariableFeatures(scData) <- unique(varFeatures)
  }
  scData@reductions$tsne_noCorrected <- Reductions(scData_noCorrected, "tsne")
  scData@reductions$umap_noCorrected <- Reductions(scData_noCorrected, "umap")
  scData@meta.data$ident_noCorrected <- Idents(scData_noCorrected)
  scData <- PrepSCTFindMarkers(scData)

  #positive cluster markers
  posMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  posMarkers[['isSpatialMarker']] = FALSE
  #spatially variable genes
  require(readxl)
  filePath_spatialMarkers <- sub('scData.rds', 'spatialMarkers.xlsx', filePath)
  spatialMarkersList <- lapply(filePath_spatialMarkers, read_xlsx)
  names(spatialMarkersList) <- names(scDataList)
  spatialMarkers <- c()
  for (j in 1:length(spatialMarkersList)) {
    spatialMarkersList[[j]][['GeneSymbol']] = spatialMarkersList[[j]]$GeneSymbol
    spatialMarkersList[[j]][['SampleID']] = names(spatialMarkersList)[j]
    spatialMarkersList[[j]] = spatialMarkersList[[j]][
      !is.na(spatialMarkersList[[j]]$MeanRank),
    ]
    spatialMarkers <- rbind(spatialMarkers, spatialMarkersList[[j]])
  }

  spatialPosMarkers <- intersect(
    posMarkers$gene,
    unique(spatialMarkers$GeneSymbol)
  )
  posMarkers[
    which(posMarkers$gene %in% spatialPosMarkers),
    'isSpatialMarker'
  ] = TRUE

  #Save some results in external files
  dataFiles = saveExternalFiles(list(
    pos_markers = posMarkers,
    spatial_markers = spatialMarkers
  ))
  saveRDS(scData, "scData.rds")
  saveRDS(param, "param.rds")

  makeRmdReport(
    dataFiles = dataFiles,
    output = output,
    rmdFile = "SpatialSeuratSlides.Rmd",
    reportTitle = param$name
  )
  return("Success")
}
