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
      methods_facts = function(param = list()) {
        humanMouse <- methodsSpeciesIs(param, c("Human", "Mouse"))
        integ <- param$integrationMethod
        corrected <- any(c("CCA", "RPCA", "Harmony") %in% integ)
        sct <- "SCTransform" %in% param$normalizationMethod
        logNorm <- "LogNormalize" %in% param$normalizationMethod
        ccReg <- isTRUE(as.logical(param$SCT.regress.CellCycle))
        ccClause <- if (ccReg) ", regressing out the cell-cycle score difference (S minus G2M)" else ""
        nfeat <- if (ezIsSpecified(param$nfeatures)) param$nfeatures else 3000
        res <- suppressWarnings(as.numeric(param$resolution[1]))
        embedding <- if ("none" %in% integ) "the PCA of the uncorrected merge"
          else if (any(c("CCA", "RPCA") %in% integ)) "the PCA of the integrated assay"
          else if ("Harmony" %in% integ) "the Harmony embedding"
          else "the PCA, integrated PCA or Harmony embedding"
        c(
          ## ezMethodScSeuratCombine set.seed(38) (app-ScSeuratCombine.R:132); SCTransform seed.use 1448145, RunPCA/RunUMAP 42, FindClusters random.seed 0 are Seurat defaults
          "The R random seed was set to 38 at the start; Seurat functions used their own default seeds (SCTransform 1448145, RunPCA and RunUMAP 42, FindClusters 0).",
          ## Condition handling in ezMethodScSeuratCombine; overwriteCondition is not on the form, so the replacement only happens for empty Condition values
          "When the loaded objects carried no Condition values, Condition was taken from the input dataset's Condition column, or set to the sample name if the dataset has none.",
          if ("Harmony" %in% integ) "When the input dataset has a single Condition value, Condition was replaced by the sample name.",
          ## seuratNormalizeSampleList + getSeuratVarsToRegress (seuratUtils.R:308-335, 781-793)
          if (sct) paste0("Each sample was normalised separately from its RNA counts with SCTransform (vst.flavor v2, 3000 variable features, Seurat defaults)", ccClause, "."),
          if (logNorm) sprintf("Each sample was normalised separately from its RNA counts: counts were log-normalised (scale factor 10000) and %s variable genes (nfeatures, an app default not on the parameter form) were selected per sample (vst).", nfeat),
          ## seuratScaleMergedLogNorm (seuratUtils.R:342-352): uncorrected merge always, Harmony input too; CCA/RPCA ScaleData per sample (seuratUtils.R:400-403)
          if (logNorm) paste0("The merged log-normalised data of the uncorrected merge", if ("Harmony" %in% integ) " and of the Harmony input" else "", " were scaled", ccClause, "."),
          if (logNorm && any(c("CCA", "RPCA") %in% integ)) "For integration, each sample was scaled on the integration features without regression.",
          ## SelectIntegrationFeatures(nfeatures = param$nfeatures) in cellClustWithCorrection only (seuratUtils.R:395-398); nfeatures is an appDefault, not declared in ScSeuratCombineApp.rb
          if (corrected) sprintf("For integration, %s features (nfeatures, an app default not on the parameter form) were selected across samples with Seurat SelectIntegrationFeatures.", nfeat),
          ## cellClustWithCorrection CCA branch (seuratUtils.R:400-445); k.anchor 5, k.filter 200, k.score 30, k.weight 100 are Seurat defaults
          if ("CCA" %in% integ) "Samples were integrated with Seurat CCA: FindIntegrationAnchors and IntegrateData (not IntegrateLayers) over dimensions 1 to npcs, with k.anchor 5, k.filter 200, k.score 30 and k.weight 100 (Seurat defaults); PCA, neighbours and clusters were then computed on the integrated assay.",
          ## cellClustWithCorrection RPCA branch (seuratUtils.R:400-445)
          if ("RPCA" %in% integ) "Samples were integrated with Seurat reciprocal PCA: per-sample PCA on the integration features, then FindIntegrationAnchors (reduction rpca, k.anchor = 20; k.filter 200, k.score 30 and k.weight 100 Seurat defaults) and IntegrateData (not IntegrateLayers) over dimensions 1 to npcs; PCA, neighbours and clusters were then computed on the integrated assay.",
          ## cellClustWithCorrection Harmony branch (seuratUtils.R:446-511); harmony 2.0.5 defaults
          if ("Harmony" %in% integ) sprintf(
            "Samples were integrated with Harmony: the normalised samples were merged, PCA was computed with npcs components on the integration features, and harmony RunHarmony corrected all npcs components for %s%s, with harmony defaults (theta 2 per variable, automatic lambda, sigma 0.1, at most 10 iterations); Batch has one level per input sample.",
            paste(if (ezIsSpecified(param$harmonyGroupBy)) param$harmonyGroupBy else "Condition", collapse = " and "),
            if (ezIsSpecified(param$additionalFactors)) " plus the additionalFactors columns" else ""),
          ## seuratIntegrateDataAndAnnotate calls cellClustNoCorrection unconditionally (seuratUtils.R:354-379)
          "An uncorrected merge (variable features = union of the per-sample variable features) was always clustered too and kept for comparison plots.",
          if ("none" %in% integ) "No integration was applied (integrationMethod none); the uncorrected merge is the reported result.",
          ## seuratStandardWorkflow (seuratUtils.R:164-240)
          paste0("The shared-nearest-neighbour graph (k = 20), UMAP (uwot, cosine metric, 30 neighbours, seed 42) and t-SNE used dimensions 1 to npcs of ", embedding, "; clusters were found with the Louvain algorithm at resolutions 0.2, 0.4, 0.6, 0.8 and 1.0 plus the resolution parameter rounded to one decimal, and the clustering at the resolution parameter is the one reported."),
          ## seuratStandardWorkflow: selectedCol not found for an unrounded resolution -> candidates[1]
          if (isTRUE(res != round(res, 1))) "Because the resolution has more than one decimal, the lowest-resolution clustering was reported instead.",
          ## PrepSCTFindMarkers only when an SCT assay exists (seuratIntegrateDataAndAnnotate); getSeuratMarkers (seuratUtils.R:795-825); pvalue_allMarkers 0.05 from EZ_PARAM_DEFAULTS.txt, min.diff.pct 0 appDefault
          sprintf("Cluster markers were found with Seurat FindAllMarkers %s, using the test in DE.method and only positive markers; p-values were Bonferroni-adjusted over all genes (Seurat p_val_adj), and markers with p_val_adj below %s (pvalue_allMarkers, not on the form) and a detection-fraction difference of at least %s (min.diff.pct, app default not on the form) were kept.",
                  if (sct) "on the SCT assay (after PrepSCTFindMarkers)" else if (logNorm) "on the RNA assay"
                  else "on the SCT assay (after PrepSCTFindMarkers) or, with LogNormalize, on the RNA assay",
                  if (ezIsSpecified(param$pvalue_allMarkers)) param$pvalue_allMarkers else 0.05,
                  if (ezIsSpecified(param$min.diff.pct)) param$min.diff.pct else 0),
          ## getSeuratMarkers passes no latent.vars (seuratUtils.R:799-809)
          if ("LR" %in% param$DE.method) "Cluster markers were tested with the LR test without latent variables; DE.regress was not applied to marker detection.",
          ## getSeuratMarkersAndAnnotate (seuratUtils.R:847-877); querySignificantClusterAnnotationEnrichR skips when no database is set (app-ScSeurat.R)
          if (humanMouse && ezIsSpecified(param$enrichrDatabase)) "Clusters were annotated with Enrichr on each cluster's significant markers (terms with adjusted p below 0.001 and more than 3 overlapping genes, top 5 per database).",
          ## cellsLabelsWithAUC skips when no tissue is set (scTools.R:521-560); AUCell aucMaxRank default 5%
          if (humanMouse && ezIsSpecified(param$tissue)) "Cells and clusters were annotated with AUCell (top 5% of ranked genes, AUCell default) using CellMarker 2.0 (2023-09-27 release) gene sets of at least 3 genes for the selected tissue.",
          ## cellsLabelsWithSingleR skips when SingleR is unset or none (scTools.R:630-660)
          if (humanMouse && ezIsSpecified(param$SingleR) && !("none" %in% param$SingleR)) "Cells and clusters were annotated with SingleR on the fine labels (label.fine) of the celldex reference named in SingleR, per cell and per cluster.",
          ## computeTFActivityAnalysis / computePathwayActivityAnalysis (app-ScSeurat.R); run_wmean seed 42 default; get_progeny top 500, progeny::getModel top 100
          if (humanMouse && isTRUE(as.logical(param$computePathwayTFActivity))) sprintf(
            "Transcription-factor and pathway activities were inferred with decoupleR run_wmean (100 permutations, minsize 5, seed 42) on the normalised data, using DoRothEA regulons of confidence A-C and PROGENy models (top %s genes per pathway).",
            if (methodsSpeciesIs(param, "Human")) "500" else "100")
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
