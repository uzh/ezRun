###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratCombine <-
  setRefClass(
    "EzAppScSeuratCombine",
    contains = "EzApp",
    methods = list(
      ## NOTE: unlike ScSeuratApp, this app does NOT declare params for per-sample QC
      ## (doublets/ambient RNA/emptyDrops) or most cell-type annotation tools
      ## (AUCell/SingleR/sc-type/mLLMCelltype/CyteTypeR/Azimuth) -- those ran upstream
      ## in ScSeuratApp already. It only reads a cached cellxgeneResults.rds if
      ## present, never computes it, so STACAS/schard aren't invoked here either.
      citation = function() {
        c(
          "Hao, Y. et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology 42, 293-304 (2024). https://doi.org/10.1038/s41587-023-01767-y",
          "Stuart, T. et al. Comprehensive Integration of Single-Cell Data. Cell 177, 1888-1902 (2019). https://doi.org/10.1016/j.cell.2019.05.031",
          "Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nature Methods 16, 1289-1296 (2019). https://doi.org/10.1038/s41592-019-0619-0",
          "Chen, E.Y. et al. Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. BMC Bioinformatics 14, 128 (2013). https://doi.org/10.1186/1471-2105-14-128",
          "Kuleshov, M.V. et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. Nucleic Acids Research 44(W1), W90-W97 (2016). https://doi.org/10.1093/nar/gkw377",
          "Badia-i-Mompel, P. et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances 2(1), vbac016 (2022). https://doi.org/10.1093/bioadv/vbac016",
          "Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D. & Saez-Rodriguez, J. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Research 29, 1363-1375 (2019). https://doi.org/10.1101/gr.240663.118",
          "Schubert, M. et al. Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications 9, 20 (2018). https://doi.org/10.1038/s41467-017-02391-6"
        )
      },
      ## Seurat defaults quoted here were checked against Seurat 5.5.1 formals() (R 4.6 system lib),
      ## harmony 2.0.5, decoupleR 2.17.0, progeny 1.34.0, AUCell 1.34.0.
      methods_facts = function() {
        c(
          ## ezMethodScSeuratCombine set.seed(38) (app-ScSeuratCombine.R:132); SCTransform seed.use 1448145, RunPCA/RunUMAP 42, FindClusters random.seed 0 are Seurat defaults
          "The R random seed was set to 38 at the start; Seurat functions used their own default seeds (SCTransform 1448145, RunPCA and RunUMAP 42, FindClusters 0).",
          ## Condition handling (app-ScSeuratCombine.R:182-199)
          "When the loaded objects carried no Condition values, Condition was taken from the input dataset's Condition column, or set to the sample name if the dataset has none; when integrationMethod is Harmony and the input dataset has a single Condition value, Condition was replaced by the sample name.",
          ## seuratNormalizeSampleList + seuratScaleMergedLogNorm + getSeuratVarsToRegress (seuratUtils.R:308-352, 781-793)
          "Each sample was normalised separately from its RNA counts: with normalizationMethod SCTransform, SCTransform was re-run per sample (vst.flavor v2, 3000 variable features, Seurat defaults), regressing out the cell-cycle score difference (S minus G2M) when SCT.regress.CellCycle is true; with LogNormalize, counts were log-normalised (scale factor 10000) and variable genes selected per sample (vst), and for Harmony and the uncorrected merge the merged data were scaled with the optional cell-cycle regression, while for CCA and RPCA each sample was scaled on the integration features without regression.",
          ## SelectIntegrationFeatures(nfeatures = param$nfeatures) (seuratUtils.R:395-398); nfeatures is an appDefault, not declared in ScSeuratCombineApp.rb
          "For integration, 3000 features (nfeatures, an app default not on the parameter form) were selected across samples with Seurat SelectIntegrationFeatures.",
          ## cellClustWithCorrection CCA/RPCA branch (seuratUtils.R:400-445); k.anchor 5, k.filter 200, k.score 30, k.weight 100 are Seurat defaults
          "With integrationMethod CCA or RPCA, samples were integrated with Seurat FindIntegrationAnchors and IntegrateData (not IntegrateLayers) over dimensions 1 to npcs; RPCA used per-sample PCA on the integration features and k.anchor = 20, CCA the default k.anchor = 5 (k.filter 200, k.score 30, k.weight 100 were Seurat defaults), and PCA, neighbours and clusters were then computed on the integrated assay.",
          ## cellClustWithCorrection Harmony branch (seuratUtils.R:446-511); harmony 2.0.5 defaults
          "With integrationMethod Harmony, the normalised samples were merged, PCA was computed with npcs components on the integration features, and harmony RunHarmony corrected all npcs components for the harmonyGroupBy column(s) plus any additionalFactors, with harmony defaults (theta 2 per variable, automatic lambda, sigma 0.1, at most 10 iterations); Batch has one level per input sample.",
          ## seuratIntegrateDataAndAnnotate calls cellClustNoCorrection unconditionally (app-ScSeuratCombine.R:310-328; seuratUtils.R:354-379)
          "An uncorrected merge (variable features = union of the per-sample variable features) was always clustered too and kept for comparison plots; with integrationMethod none it is the reported result.",
          ## seuratStandardWorkflow (seuratUtils.R:164-240)
          "The shared-nearest-neighbour graph (k = 20), UMAP (uwot, cosine metric, 30 neighbours, seed 42) and t-SNE used dimensions 1 to npcs of the PCA, integrated PCA or Harmony embedding; clusters were found with the Louvain algorithm at resolutions 0.2, 0.4, 0.6, 0.8 and 1.0 plus the resolution parameter rounded to one decimal, and the clustering at the resolution parameter is the one reported (when the resolution has more than one decimal, the lowest-resolution clustering was reported instead).",
          ## PrepSCTFindMarkers (app-ScSeuratCombine.R:329-332); getSeuratMarkers (seuratUtils.R:795-825); pvalue_allMarkers 0.05 from EZ_PARAM_DEFAULTS.txt, min.diff.pct 0 appDefault
          "Cluster markers were found with Seurat FindAllMarkers on the SCT assay (after PrepSCTFindMarkers) or, with LogNormalize, on the RNA assay, using the test in DE.method and only positive markers; p-values were Bonferroni-adjusted over all genes (Seurat p_val_adj), and markers with p_val_adj below 0.05 (ezRun default pvalue_allMarkers, not on the form) and any detection-fraction difference (min.diff.pct 0, app default) were kept.",
          ## getSeuratMarkers passes no latent.vars (seuratUtils.R:799-809)
          "When DE.method is LR, cluster markers were tested without latent variables; DE.regress was not applied to marker detection.",
          ## getSeuratMarkersAndAnnotate (seuratUtils.R:847-877); app-ScSeurat.R Enrichr query; scTools.R:521-660; AUCell aucMaxRank default 5%
          "For human and mouse data, clusters were annotated with Enrichr on each cluster's significant markers (terms with adjusted p below 0.001 and more than 3 overlapping genes, top 5 per database), with AUCell (top 5% of ranked genes, AUCell default) using CellMarker 2.0 (2023-09-27 release) gene sets of at least 3 genes for the selected tissue, and, when a SingleR reference is set, with SingleR on the fine labels (label.fine) of that celldex reference, per cell and per cluster.",
          ## computeTFActivityAnalysis / computePathwayActivityAnalysis (app-ScSeurat.R); run_wmean seed 42 default; get_progeny top 500, progeny::getModel top 100
          "When computePathwayTFActivity is true (human and mouse), transcription-factor and pathway activities were inferred with decoupleR run_wmean (100 permutations, minsize 5, seed 42) on the normalised data, using DoRothEA regulons of confidence A-C and PROGENy models (top 500 genes per pathway for human, top 100 for mouse)."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodScSeuratCombine
        name <<- "EzAppScSeuratCombine"
        appDefaults <<- rbind(
          nfeatures = ezFrame(
            Type = "numeric",
            DefaultValue = 3000,
            Description = "number of variable genes for SCT"
          ),
          npcs = ezFrame(
            Type = "numeric",
            DefaultValue = 30,
            Description = "The maximal dimensions to use for reduction"
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
          integrationMethod = ezFrame(
            Type = "character",
            DefaultValue = "Harmony",
            Description = "Choose integration method in Seurat (Harmony, CCA, RPCA)"
          ),
          normalizationMethod = ezFrame(
            Type = "charVector",
            DefaultValue = "SCTransform",
            Description = "SCTransform (variance-stabilizing residuals) or LogNormalize (classic log1p of counts per 10k)"
          ),
          harmonyGroupBy = ezFrame(
            Type = "charVector",
            DefaultValue = "Condition",
            Description = "Metadata column(s) Harmony corrects for; Batch is one level per input sample"
          ),
          enrichrDatabase = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "enrichR databases to search"
          ),
          computePathwayTFActivity = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "Whether we should compute pathway and TF activities."
          ),
          SCT.regress.CellCycle = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Choose CellCycle to be regressed out (in SCTransform, or in ScaleData when LogNormalize is used) if it is a bias."
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
          min.diff.pct = ezFrame(
            Type = "numeric",
            DefaultValue = 0,
            Description = "Used in filtering cluster markers"
          ),
          logfc.threshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.25,
            Description = "Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
          )
        )
      }
    )
  )

ezMethodScSeuratCombine = function(
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
  library(AUCell)
  library(decoupleR)
  library(BiocParallel)
  require(future)

  plan("multicore", workers = param$cores)
  set.seed(38)
  future.seed = TRUE
  options(future.rng.onMisuse = "ignore")
  options(future.globals.maxSize = param$ram * 1024^3)

  BPPARAM <- MulticoreParam(workers = param$cores)
  ## Pin BLAS/OpenMP to one thread before forking (MulticoreParam/future) to
  ## avoid the fork-in-multithreaded-process deadlock (e.g. AUCell labeling).
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  register(BPPARAM)

  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  reportCwd <- getwd()
  ezLog("Attempting to load Seurat data...")
  filePath <- file.path("/srv/gstore/projects", input$getColumn("SC Seurat"))
  filePath_course <- file.path(
    "/srv/GT/analysis/course_sushi/public/projects",
    input$getColumn("SC Seurat")
  )
  if (!file.exists(filePath[1])) {
    filePath <- filePath_course
  }
  names(filePath) <- input$getNames()

  if (length(filePath) < 2) {
    stop("need at least two samples to combine.")
  }

  # Load the data and prepare metadata for integration
  scDataList <- lapply(names(filePath), function(sm) {
    scData <- ezLoadRobj(filePath[sm], nthreads = param$cores)
    aziFilePath <- file.path(dirname(filePath[sm]), 'aziResults.rds')
    if (file.exists(aziFilePath)) {
      aziResults <- readRDS(aziFilePath)
      scData <- AddMetaData(scData, aziResults)
    }
    # Load cellxgene results if available
    cellxgeneFilePath <- file.path(
      dirname(filePath[sm]),
      'cellxgeneResults.rds'
    )
    if (file.exists(cellxgeneFilePath)) {
      cellxgeneResults <- readRDS(cellxgeneFilePath)
      scData <- AddMetaData(scData, cellxgeneResults)
    }
    scData$Sample <- sm
    # If we have new information in the Condition column, add it to the dataset
    if (
      all(scData$Condition == "NA" | scData$Condition == "") ||
        (ezIsSpecified(param$overwriteCondition) &&
          as.logical(param$overwriteCondition))
    ) {
      if (any(startsWith(colnames(input$meta), "Condition"))) {
        scData$Condition <- unname(input$getColumn("Condition")[sm])
      } else {
        scData$Condition <- scData$Sample
      }
    }
    # Harmony will complain if the Condition is the same across all samples
    if (
      param$integrationMethod == "Harmony" &&
        length(unique(input$meta$`Condition`)) == 1
    ) {
      scData$Condition <- scData$Sample
    }
    # Also add the other factors in the input dataset to the objects
    if (ezIsSpecified(param$additionalFactors)) {
      additionalFactors <- str_split(
        param$additionalFactors,
        ",",
        simplify = TRUE
      )[1, ]
      metaFactorNames <- paste0("meta_", additionalFactors) %>%
        str_replace(., " ", ".")
      names(metaFactorNames) <- additionalFactors
      for (hf in additionalFactors) {
        scData[[metaFactorNames[hf]]] <- unname(input$getColumn(hf)[sm])
      }
    }
    # Rename the cells and add original sample-level clusters back in
    scData <- RenameCells(
      scData,
      new.names = paste0(scData$Sample, "-", colnames(scData))
    )
    if (class(scData$seurat_clusters) == 'character') {
      scData$sample_seurat_clusters <- paste0(
        scData$Sample,
        "-",
        sprintf("%s", scData$seurat_clusters)
      )
    } else {
      scData$sample_seurat_clusters <- paste0(
        scData$Sample,
        "-",
        sprintf("%02d", scData$seurat_clusters)
      )
    }
    return(scData)
  })

  # perform all of the analysis
  results <- seuratIntegrateDataAndAnnotate(
    scDataList,
    input,
    output,
    param,
    BPPARAM
  )

  # generate ClusterInfos table
  clusterInfos <- ezFrame(
    Samples = paste(input$getNames(), collapse = ','),
    Cluster = levels(Idents(results$scData)),
    ClusterLabel = ""
  )
  if (!is.null(results$singler.results)) {
    clusterInfos$SinglerCellType <- results$singler.results$singler.results.cluster[
      clusterInfos$Cluster,
      "pruned.labels"
    ]
  }
  nTopMarkers <- 10
  topMarkers <- results$markers %>%
    group_by(cluster) %>%
    slice_max(n = nTopMarkers, order_by = avg_log2FC)
  topMarkerString <- sapply(
    split(topMarkers$gene, topMarkers$cluster),
    paste,
    collapse = ", "
  )
  clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
  clusterInfoFile <- "clusterInfos.xlsx"
  writexl::write_xlsx(clusterInfos, path = clusterInfoFile)

  # save the markers
  writexl::write_xlsx(results$markers, path = "posMarkers.xlsx")
  qs2::qs_save(results$scData, "scData.qs2", nthreads = param$cores)

  # Save some results in external files
  reportTitle <- 'SCReport - MultipleSamples based on Seurat'
  makeRmdReport(
    param = param,
    output = output,
    scData = results$scData,
    enrichRout = results$enrichRout,
    TFActivity = results$TFActivity,
    pathwayActivity = results$pathwayActivity,
    aziResults = results$aziResults,
    cellxgeneResults = results$cellxgeneResults,
    cells.AUC = results$cells.AUC,
    singler.results = results$singler.results,
    rmdFile = "ScSeuratCombine.Rmd",
    reportTitle = reportTitle
  )
  return("Success")
}

seuratIntegrateDataAndAnnotate <- function(
  scDataList,
  input,
  output,
  param,
  BPPARAM = SerialParam()
) {
  pvalue_allMarkers <- param$pvalue_allMarkers

  if (ezIsSpecified(param$chosenClusters)) {
    for (eachSample in names(param$chosenClusters)) {
      chosenCells <- names(Idents(scDataList[[eachSample]]))[
        Idents(scDataList[[eachSample]]) %in% param$chosenClusters[[eachSample]]
      ]
      scDataList[[eachSample]] <- scDataList[[eachSample]][, chosenCells]
    }
  }

  scData_noCorrected <- cellClustNoCorrection(scDataList, param)
  if (param$integrationMethod != 'none') {
    scData_corrected = cellClustWithCorrection(scDataList, param)
    #in order to compute the markers we switch again to the original assay
    DefaultAssay(scData_corrected) <- seuratAnalysisAssay(param)
    scData <- scData_corrected
  } else {
    scData = scData_noCorrected
  }
  scData@reductions$tsne_noCorrected <- Reductions(scData_noCorrected, "tsne")
  Key(scData@reductions$tsne_noCorrected) <- 'TSNEnoCorrection_'
  scData@reductions$umap_noCorrected <- Reductions(scData_noCorrected, "umap")
  Key(scData@reductions$umap_noCorrected) <- 'UMAPnoCorrection_'

  # Since Seurat v5, the RNA layers are split by sample. Does not affect SCT assay
  # but clutters up RNA assay should it be used in downstream analysis
  scData <- JoinLayers(scData, assay = "RNA")

  scData@meta.data$ident_noCorrected <- Idents(scData_noCorrected)
  ## Only the SCTransform path has SCT residuals to re-prepare
  if ("SCT" %in% Seurat::Assays(scData)) {
    scData <- PrepSCTFindMarkers(scData)
  }

  # get annotation information
  anno <- getSeuratMarkersAndAnnotate(scData, param, BPPARAM = BPPARAM)

  return(list(
    scData = scData,
    markers = anno$markers,
    enrichRout = anno$enrichRout,
    pathwayActivity = anno$pathwayActivity,
    TFActivity = anno$TFActivity,
    cells.AUC = anno$cells.AUC,
    singler.results = anno$singler.results,
    aziResults = anno$aziResults,
    cellxgeneResults = anno$cellxgeneResults
  ))
}
