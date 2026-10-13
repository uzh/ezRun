###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## Per-sample Visium SD analysis: spatially aware QC, SCT, expression
## clusters, BANKSY spatial domains with region annotation, spatially
## variable genes, RCTD (or CARD) deconvolution, domain neighbourhoods,
## pathway/TF activity, distance-constrained ligand-receptor analysis
## (CellChat v2; optional SpatialCellChat, liana), MISTy cell-type
## dependencies and optional BayesSpace domains. The step functions live in
## R/spatialSeuratSteps.R, R/spatialCommunication.R and R/spatialStats.R.

EzAppSpatialSeurat <-
  setRefClass(
    "EzAppSpatialSeurat",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSpatialSeurat
        name <<- "EzAppSpatialSeurat"
        appDefaults <<- rbind(
          nfeatures = ezFrame(
            Type = "numeric",
            DefaultValue = 3000,
            Description = "number of variable genes for SCT"
          ),
          npcs = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "The maximal dimensions to use for reduction"
          ),
          SCT.regress.CellCycle = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Regress out the cell cycle in SCTransform (implies cellCycle)"
          ),
          cellCycle = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Score the cell cycle phase of every spot"
          ),
          enrichrDatabase = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "enrichR databases to search"
          ),
          DE.method = ezFrame(
            Type = "charVector",
            DefaultValue = "wilcox",
            Description = "Test used for cluster and domain markers (wilcox or LR)"
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
          resolution = ezFrame(
            Type = "numeric",
            DefaultValue = 0.6,
            Description = "Resolution of the expression clustering"
          ),
          cellsFraction = ezFrame(
            Type = "numeric",
            DefaultValue = 0,
            Description = "A gene will be kept if it is expressed in at least this percentage of cells"
          ),
          nUMIs = ezFrame(
            Type = "numeric",
            DefaultValue = 1,
            Description = 'A gene will be kept if it has at least nUMIs in the fraction of cells specified before'
          ),
          nmad = ezFrame(
            Type = "numeric",
            DefaultValue = 3,
            Description = "Median absolute deviation (MAD) from the median value of each metric across all cells"
          ),
          nreads = ezFrame(
            Type = "numeric",
            DefaultValue = NA,
            Description = "Fixed threshold: spots with fewer UMIs are removed; NA uses SpotSweeper local outliers"
          ),
          ngenes = ezFrame(
            Type = "numeric",
            DefaultValue = NA,
            Description = "Fixed threshold: spots with fewer genes are removed; NA uses SpotSweeper local outliers"
          ),
          perc_mito = ezFrame(
            Type = "numeric",
            DefaultValue = NA,
            Description = "Fixed threshold: spots with a higher mitochondrial percentage are removed; NA uses SpotSweeper local outliers"
          ),
          perc_ribo = ezFrame(
            Type = "numeric",
            DefaultValue = NA,
            Description = "Fixed threshold: spots with a higher ribosomal protein percentage are removed; NA keeps all"
          ),
          spotClean = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Run spotClean method"
          ),
          pvalue_allMarkers = ezFrame(
            Type = "numeric",
            DefaultValue = 0.01,
            Description = "pValue for marker detection"
          ),
          featSelectionMethod = ezFrame(
            Type = "character",
            DefaultValue = "vst",
            Description = "vst (default) or STACAS (blocklist genes; needs the SignatuR package)"
          ),
          pt.size.factor = ezFrame(
            Type = "numeric",
            DefaultValue = NA,
            Description = "pt.size.factor for spatial plots"
          ),
          banksyLambda = ezFrame(
            Type = "numeric",
            DefaultValue = 0.8,
            Description = "BANKSY weight of the neighbourhood (0 = expression only, 0.8 = domain segmentation)"
          ),
          banksyKgeom = ezFrame(
            Type = "numeric",
            DefaultValue = 18,
            Description = "BANKSY number of spatial neighbours (18 = two hexagonal rings of Visium spots)"
          ),
          banksyResolution = ezFrame(
            Type = "numeric",
            DefaultValue = 0.5,
            Description = "Resolution of the BANKSY domain clustering"
          ),
          banksyNpcs = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "Number of BANKSY principal components"
          ),
          regionMarkerFile = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "TSV (region, gene) used to name the domains; empty uses the built-in mouse brain panel"
          ),
          computeNnSVG = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Rank the top Moran's I genes additionally with nnSVG"
          ),
          nnsvgTopN = ezFrame(
            Type = "numeric",
            DefaultValue = 300,
            Description = "Number of top Moran's I genes tested with nnSVG"
          ),
          svgPatterns = ezFrame(
            Type = "numeric",
            DefaultValue = 6,
            Description = "Number of spatial expression patterns the top SVGs are grouped into"
          ),
          rctdReference = ezFrame(
            Type = "character",
            DefaultValue = "None",
            Description = "RCTD reference relative to /srv/GT/databases/RCTD_References"
          ),
          rctdFile = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Full path to a custom RCTD reference (overrides rctdReference)"
          ),
          rctdUMImin = ezFrame(
            Type = "numeric",
            DefaultValue = 100,
            Description = "Minimal UMIs of a spot to be deconvolved"
          ),
          rctdEngine = ezFrame(
            Type = "character",
            DefaultValue = "rctd-py",
            Description = "RCTD implementation: rctd-py (same model, ~1.5 instead of ~10 min per section) or spacexr (R)"
          ),
          deconvolutionMethod = ezFrame(
            Type = "character",
            DefaultValue = "RCTD",
            Description = "Cell-type deconvolution with the selected reference: RCTD, or CARD (spatially smoothed, ~3 min, 9 GB; needs the MuSiC package)"
          ),
          computePathwayTFActivity = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Score PROGENy pathway and DoRothEA TF activities per spot"
          ),
          computeCellChat = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Run distance-constrained CellChat between domains"
          ),
          cellchatInteractionLength = ezFrame(
            Type = "numeric",
            DefaultValue = 250,
            Description = "Maximal ligand-receptor interaction distance in microns"
          ),
          cellchatNboot = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "Permutations for the CellChat and SpatialCellChat p-values"
          ),
          computeMisty = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "MISTy: explain each cell type's proportion from the same spot, the direct neighbours and the wider neighbourhood (needs the deconvolution)"
          ),
          mistyParaviewRange = ezFrame(
            Type = "numeric",
            DefaultValue = 300,
            Description = "Length scale in microns of the MISTy neighbourhood (paraview)"
          ),
          computeSpatialCellChat = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "SpatialCellChat: communication between neighbouring spots, averaged to cell types with the deconvolution weights (adds ~11 min per section on 8 cores)"
          ),
          computeLiana = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "liana consensus ligand-receptor ranking between touching domains (cross-check of CellChat; ~1.5 min)"
          ),
          computeBayesSpace = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "BayesSpace domains with as many domains as BANKSY found, for comparison (~3 min with 10000 iterations)"
          ),
          bayesspaceNrep = ezFrame(
            Type = "numeric",
            DefaultValue = 10000,
            Description = "BayesSpace MCMC iterations"
          ),
          writeCloupe = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Write a Loupe Browser file with clusters and domains (Space Ranger already provides cloupe.cloupe)"
          )
        )
      }
    )
  )

ezMethodSpatialSeurat <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
  library(Seurat)
  library(future)
  library(BiocParallel)
  stopifnot(length(input$getNames()) == 1)
  sampleName <- input$getNames()

  if (param$cores > 1) {
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    BPPARAM <- SerialParam()
  }
  ## Pin BLAS/OpenMP to one thread before forking (MulticoreParam/future) to
  ## avoid the fork-in-multithreaded-process deadlock.
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  register(BPPARAM)
  plan("multicore", workers = param$cores)
  set.seed(38)
  options(future.rng.onMisuse = "ignore")
  options(future.globals.maxSize = param$ram * 1024^3)
  if (ezIsSpecified(param$SCT.regress.CellCycle) && param$SCT.regress.CellCycle) {
    param$cellCycle <- TRUE
  }
  status <- new.env()
  species <- getSpecies(param$refBuild)

  ## 1. load ----
  featInfo <- spatial_feature_info(input)
  loaded <- spatial_load(input, param)
  scData <- loaded$scData
  scDataRaw <- loaded$scDataRaw
  param <- loaded$param
  spotInfo <- loaded$spotInfo
  remove(loaded)

  ## 2. QC ----
  qc <- spatial_qc(scData, featInfo, param, spotInfo)
  scData <- qc$scData
  scData.unfiltered <- qc$scData.unfiltered
  cellsPerGeneFraction <- qc$cellsPerGeneFraction
  qcInfo <- qc$qcInfo
  remove(qc)

  ## 3. normalisation and expression clusters ----
  scData <- spatial_expr_clusters(scData, param)
  scData <- spatial_add_graphs(scData)
  exprMarkers <- spatial_group_markers(scData, "expr_cluster", param)

  ## 4. BANKSY domains ----
  banksy <- spatial_banksy_domains(scData, param)
  scData <- banksy$scData
  banksySweep <- banksy$sweep
  remove(banksy)
  domainMarkers <- spatial_group_markers(scData, "banksy_domain", param)

  ## 5. optional steps ----
  domainAnno <- run_optional_step(
    "domainAnnotation",
    spatial_annotate_domains(scData, param),
    status
  )
  if (!is.null(domainAnno)) {
    scData <- domainAnno$scData
    domainAnno$scData <- NULL
  }
  svg <- run_optional_step("spatiallyVariableGenes", spatial_svg(scData, param), status)
  if (!is.null(svg)) {
    scData <- svg$scData
    svg$scData <- NULL
  }
  rctd <- if (is.null(resolve_rctd_reference_path(param))) {
    skip_step("deconvolution", "no RCTD reference selected", status)
  } else if (identical(param$deconvolutionMethod, "CARD")) {
    run_optional_step("deconvolution", spatial_card(scData, featInfo, param), status)
  } else {
    run_optional_step("deconvolution", spatial_rctd(scData, featInfo, param), status)
  }
  if (!is.null(rctd)) {
    scData <- rctd$scData
    rctd$scData <- NULL
  }
  neighbourhood <- run_optional_step(
    "neighbourhood",
    spatial_neighbourhood(scData, rctd$weights),
    status
  )
  activity <- if (ezIsSpecified(param$computePathwayTFActivity) && param$computePathwayTFActivity) {
    run_optional_step("activity", spatial_activity(scData, param), status)
  } else {
    skip_step("activity", "disabled by parameter", status)
  }
  if (!is.null(activity)) {
    scData <- activity$scData
    activity$scData <- NULL
  }
  communication <- if (ezIsSpecified(param$computeCellChat) && param$computeCellChat) {
    run_optional_step("communication", spatial_communication(scData, param, spotInfo, rctd$weights), status)
  } else {
    skip_step("communication", "disabled by parameter", status)
  }
  ## steps on the cell-type proportions
  misty <- run_weights_step(
    "misty",
    param$computeMisty,
    rctd,
    spatial_misty(scData, param, spotInfo, rctd$weights),
    status
  )
  spatialCellChat <- run_weights_step(
    "spatialCellChat",
    param$computeSpatialCellChat,
    rctd,
    spatial_spatialcellchat(scData, param, spotInfo, rctd$weights),
    status
  )
  liana <- if (ezIsSpecified(param$computeLiana) && param$computeLiana) {
    run_optional_step(
      "liana",
      spatial_liana(scData, param, neighbourhood$domainAdjacency$observed),
      status
    )
  } else {
    skip_step("liana", "disabled by parameter", status)
  }
  bayesspace <- if (ezIsSpecified(param$computeBayesSpace) && param$computeBayesSpace) {
    run_optional_step("bayesspace", spatial_bayesspace(scData, param), status)
  } else {
    skip_step("bayesspace", "disabled by parameter", status)
  }
  if (!is.null(bayesspace)) {
    scData <- bayesspace$scData
    bayesspace$scData <- NULL
  }
  enrichRout <- if (species %in% c("Human", "Mouse") && ezIsSpecified(param$enrichrDatabase)) {
    run_optional_step(
      "enrichR",
      {
        top <- exprMarkers %>%
          group_by(cluster) %>%
          slice_max(n = 500, order_by = avg_log2FC)
        querySignificantClusterAnnotationEnrichR(
          split(as.character(top$gene), top$cluster),
          param$enrichrDatabase
        )
      },
      status
    )
  } else {
    skip_step("enrichR", "no database selected or unsupported species", status)
  }
  aziResults <- if (ezIsSpecified(param$Azimuth) && param$Azimuth != "none") {
    run_optional_step("azimuth", spatial_azimuth(scData, param), status)
  } else {
    skip_step("azimuth", "no Azimuth reference selected", status)
  }
  if (!is.null(aziResults)) {
    scData <- AddMetaData(scData, aziResults)
  }

  ## 6. tables ----
  writexl::write_xlsx(exprMarkers, path = "posMarkers.xlsx")
  writexl::write_xlsx(domainMarkers, path = "domainMarkers.xlsx")
  if (!is.null(svg)) {
    writexl::write_xlsx(svg$table, path = "spatiallyVariableGenes.xlsx")
  }
  if (!is.null(domainAnno)) {
    writexl::write_xlsx(
      list(
        labels = domainAnno$domainLabels,
        scores = data.frame(
          domain = rownames(domainAnno$domainScores),
          domainAnno$domainScores,
          check.names = FALSE
        )
      ),
      path = "domainAnnotation.xlsx"
    )
  }
  if (!is.null(rctd)) {
    rctd$weightsFile <- paste0(tolower(rctd$assay), "Weights.xlsx")
    writexl::write_xlsx(
      list(
        weights = data.frame(spot = rownames(rctd$weights), rctd$weights, check.names = FALSE),
        compositionPerDomain = rctd$composition
      ),
      path = rctd$weightsFile
    )
  }
  if (!is.null(communication)) {
    sheets <- list()
    for (level in names(communication)) {
      sheets[[paste0("interactions_", level)]] <- communication[[level]]$interactions
      sheets[[paste0("pathways_", level)]] <- communication[[level]]$pathways
    }
    writexl::write_xlsx(sheets, path = "cellchat.xlsx")
  }
  if (!is.null(misty)) {
    writexl::write_xlsx(
      misty[c("improvements", "contributions", "importances")],
      path = "misty.xlsx"
    )
  }
  if (!is.null(spatialCellChat)) {
    writexl::write_xlsx(
      spatialCellChat[c("interactions", "pathways")],
      path = "spatialCellChat.xlsx"
    )
  }
  if (!is.null(liana)) {
    writexl::write_xlsx(list(liana = liana), path = "liana.xlsx")
  }
  bulkSignalPerDomain <- data.frame(
    GeneSymbol = rownames(scData),
    AggregateExpression(scData, assays = "SCT", group.by = "banksy_domain")$SCT,
    check.names = FALSE
  )
  dataFiles <- saveExternalFiles(list(bulkSignalPerDomain = bulkSignalPerDomain))
  ezLog("Tables written")

  ## 7. objects and report ----
  ezLog("Finalising objects")
  ## the BANKSY matrix is only needed for clustering; keep the reductions
  if ("BANKSY" %in% SeuratObject::Assays(scData)) {
    scData[["BANKSY"]] <- NULL
  }
  names(scData@images)[1] <- paste0(sampleName, "_S1")
  Idents(scData) <- "banksy_domain"
  allCellsMeta <- scData.unfiltered@meta.data
  allCellsMeta$useCell <- !allCellsMeta$discard
  if (param$spotClean) {
    qs2::qs_save(scDataRaw, "scData.raw.qs2", nthreads = 4)
  }
  if (ezIsSpecified(param$writeCloupe) && param$writeCloupe) {
    run_optional_step("cloupe", spatial_write_cloupe(scData, sampleName), status)
  } else {
    skip_step("cloupe", "disabled by parameter", status)
  }
  if (!is.null(communication)) {
    qs2::qs_save(lapply(communication, `[[`, "cellchat"), "cellchat.qs2", nthreads = 4)
    for (level in names(communication)) {
      communication[[level]]$cellchat <- NULL
    }
  }
  results <- list(
    sampleName = sampleName,
    spotInfo = spotInfo,
    qcInfo = qcInfo,
    banksySweep = banksySweep,
    exprMarkers = exprMarkers,
    domainMarkers = domainMarkers,
    domainAnno = domainAnno,
    svg = svg,
    rctd = rctd,
    neighbourhood = neighbourhood,
    activity = activity,
    communication = communication,
    misty = misty,
    spatialCellChat = spatialCellChat,
    liana = liana,
    bayesspace = bayesspace,
    enrichRout = enrichRout,
    status = as.list(status)
  )
  ezLog("Saving objects and rendering the report")
  makeQuartoReport(
    param = param,
    output = output,
    dataFiles = dataFiles,
    results = results,
    scData = scData,
    scData.unfiltered = scData.unfiltered,
    allCellsMeta = allCellsMeta,
    cellsPerGeneFraction = cellsPerGeneFraction,
    qmdFile = "SpatialSeurat.qmd",
    reportTitle = paste0(param$name, ": ", sampleName),
    number = TRUE,
    quiet = FALSE
  )
  return("Success")
}

MyDietSeurat <- function(
  object,
  layers = NULL,
  features = NULL,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE,
  counts = deprecated(),
  data = deprecated(),
  scale.data = deprecated(),
  ...
) {
  CheckDots(...)
  dep.args <- c(counts = counts, data = data, scale.data = scale.data)
  for (lyr in names(x = dep.args)) {
    if (is_present(arg = dep.args[[lyr]])) {
      if (is.null(x = layers)) {
        layers <- unique(
          x = unlist(
            x = lapply(X = Assays(object = object), FUN = function(x) {
              return(Layers(object = object[[x]]))
            })
          )
        )
      }
      deprecate_soft(
        when = "5.0.0",
        what = paste0("DietSeurat(", lyr, " = )"),
        with = "DietSeurat(layers = )"
      )
      layers <- if (isTRUE(x = dep.args[[lyr]])) {
        c(layers, lyr)
      } else {
        Filter(f = function(x) x != lyr, x = layers)
      }
    }
  }
  object <- UpdateSlots(object = object)
  assays <- assays %||% Assays(object = object)
  assays <- intersect(x = assays, y = Assays(object = object))
  if (!length(x = assays)) {
    abort(message = "No assays provided were found in the Seurat object")
  }
  if (!DefaultAssay(object = object) %in% assays) {
    abort(
      message = "The default assay is slated to be removed, please change the default assay"
    )
  }
  layers <- layers %||% assays
  layers <- .PropagateList(x = layers, names = assays)
  for (assay in names(x = layers)) {
    layers[[assay]] <- tryCatch(
      expr = Layers(object = object[[assay]], search = layers[[assay]]),
      error = function(...) {
        return(character(length = 0L))
      }
    )
  }
  layers <- Filter(f = length, x = layers)
  if (!length(x = layers)) {
    abort(message = "None of the requested layers found")
  }
  for (assay in Assays(object = object)) {
    if (!(assay %in% assays)) {
      object[[assay]] <- NULL
      next
    }
    layers.rm <- setdiff(
      x = Layers(object = object[[assay]]),
      y = layers[[assay]]
    )
    if (length(x = layers.rm)) {
      if (
        inherits(x = object[[assay]], what = "Assay") &&
          all(c("counts", "data") %in% layers.rm)
      ) {
        abort(message = "Cannot remove both 'counts' and 'data' from v3 Assays")
      }
      for (lyr in layers.rm) {
        suppressWarnings(
          object <- tryCatch(
            expr = {
              object[[assay]][[lyr]] <- NULL
              object
            },
            error = function(e) {
              if (lyr == "data") {
                object[[assay]][[lyr]] <- sparseMatrix(
                  i = 1,
                  j = 1,
                  x = 1,
                  dims = dim(object[[assay]][[lyr]]),
                  dimnames = dimnames(object[[assay]][[lyr]])
                )
              } else {
                slot(object = object[[assay]], name = lyr) <- new(
                  Class = "dgCMatrix"
                )
              }
              ezLog(
                "Converting layer ",
                lyr,
                " in assay ",
                assay,
                " to empty dgCMatrix"
              )
              object
            }
          )
        )
      }
    }
    if (!is.null(x = features)) {
      features.assay <- intersect(
        x = features,
        y = rownames(x = object[[assay]])
      )
      if (!length(x = features.assay)) {
        warn(
          message = paste0(
            "No features found in assay ",
            sQuote(x = assay),
            ", removing..."
          )
        )
        object[[assay]] <- NULL
        next
      }

      fixModel = FALSE
      if ('SCTModel.list' %in% slotNames(object[[assay]])) {
        myModel <- object[[assay]]@SCTModel.list$refmodel
        fixModel = TRUE
      }
      suppressWarnings(
        object[[assay]] <- subset(
          x = object[[assay]],
          features = features.assay
        )
      )
      if (fixModel) {
        object[[assay]]@SCTModel.list$refmodel <- myModel
      }
    }
  }
  if (!isTRUE(x = misc)) {
    slot(object = object, name = "misc") <- list()
  }
  all.objects <- .FilterObjects(
    object = object,
    classes.keep = c("DimReduc", "Graph")
  )
  objects.to.remove <- all.objects[!all.objects %in% c(dimreducs, graphs)]
  for (ob in objects.to.remove) {
    object[[ob]] <- NULL
  }
  cells.keep <- list()
  for (assay in Assays(object = object)) {
    cells.keep[[assay]] <- colnames(x = object[[assay]])
  }
  cells.keep <- intersect(colnames(x = object), unlist(cells.keep))
  if (length(cells.keep) < ncol(x = object)) {
    object <- subset(object, cells = cells.keep)
  }
  return(object)
}
