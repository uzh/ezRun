###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

seuratStandardSCTPreprocessing <- function(
  scData,
  param,
  assay = "RNA",
  seed = 38
) {
  DefaultAssay(scData) <- assay
  ## Get information on which variables to regress out in scaling/SCT
  vars.to.regress <- getSeuratVarsToRegress(param)
  ## generate normalized slots for the RNA assay
  scData <- NormalizeData(
    scData,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  species <- getSpecies(param$refBuild)
  if (
    ezIsSpecified(param$featSelectionMethod) &&
      param$featSelectionMethod == 'STACAS'
  ) {
    require(STACAS)
    require(SignatuR)
    data(SignatuR)
    if (species == 'Human') {
      my.genes.blocklist <- c(
        GetSignature(SignatuR$Hs$Blocklists),
        GetSignature(SignatuR$Hs$Compartments)
      )
      scData <- FindVariableFeatures.STACAS(
        scData,
        nfeat = param$nfeatures,
        genesBlockList = my.genes.blocklist
      )
    } else if (species == 'Mouse') {
      my.genes.blocklist <- c(
        GetSignature(SignatuR$Mm$Blocklists),
        GetSignature(SignatuR$Mm$Compartments)
      )
      scData <- FindVariableFeatures.STACAS(
        scData,
        nfeat = param$nfeatures,
        genesBlockList = my.genes.blocklist
      )
    } else {
      ezLog(
        'Selection method STACAS not supported for this species! Use default method instead.'
      )
      scData <- FindVariableFeatures(
        scData,
        selection.method = "vst",
        verbose = FALSE,
        nfeatures = param$nfeatures
      )
    }
  } else {
    scData <- FindVariableFeatures(
      scData,
      selection.method = "vst",
      verbose = FALSE,
      nfeatures = param$nfeatures
    )
  }

  scData <- ScaleData(
    scData,
    vars.to.regress = vars.to.regress,
    verbose = FALSE,
    do.scale = FALSE
  )
  ## generate the SCT assay
  scData <- SCTransform(
    scData,
    vst.flavor = "v2",
    vars.to.regress = vars.to.regress,
    assay = assay,
    seed.use = seed,
    verbose = FALSE,
    return.only.var.genes = FALSE
  )
  return(scData)
}

seuratClustering <- function(scData, param) {
  set.seed(38)
  scData <- FindVariableGenes(
    object = scData,
    do.plot = FALSE,
    x.low.cutoff = param$x.low.cutoff,
    x.high.cutoff = param$x.high.cutoff,
    y.cutoff = param$y.cutoff
  )

  scData <- ScaleData(
    object = scData,
    do.par = TRUE,
    vars.to.regress = param$vars.to.regress,
    num.cores = param$cores
  )
  if (ezIsSpecified(param$pcGenes)) {
    indicesMatch <- match(
      toupper(param$pcGenes),
      toupper(rownames(scData@data))
    )
    if (any(is.na(indicesMatch))) {
      stop(
        "The following genes don't exist: ",
        paste(param$pcGenes[is.na(indicesMatch)], collapse = ",")
      )
    }
    pc.genes <- rownames(scData@data)[indicesMatch]
  } else {
    pc.genes <- scData@var.genes
  }
  scData <- RunPCA(
    object = scData,
    pc.genes = pc.genes,
    pcs.compute = 20,
    do.print = FALSE,
    pcs.print = 1:5,
    genes.print = 5
  )
  scData <- ProjectPCA(object = scData, do.print = FALSE)
  scData <- JackStraw(object = scData, num.replicate = 100) #, display.progress=FALSE,
  #do.par=TRUE, num.cores=min(4L, param$cores))

  scData <- FindClusters(
    object = scData,
    reduction.type = "pca",
    dims.use = 1:min(param$pcs, length(pc.genes) - 1),
    resolution = param$resolution,
    print.output = 0,
    save.SNN = TRUE,
    force.recalc = FALSE
  )
  scData <- RunTSNE(
    object = scData,
    reduction.use = "pca",
    check_duplicates = FALSE,
    dims.use = 1:min(param$pcs, length(pc.genes) - 1),
    tsne.method = "Rtsne",
    perplexity = ifelse(length(scData@ident) > 200, 30, 10),
    num_threads = param$cores
  )
  try(
    scData <- RunUMAP(
      object = scData,
      reduction.use = "pca",
      dims.use = 1:min(param$pcs, length(pc.genes) - 1),
      n_neighbors = ifelse(length(scData@ident) > 200, 30, 10)
    ),
    silent = TRUE
  )
  return(scData)
}

seuratStandardWorkflow <- function(
  scData,
  param,
  reduction = "pca",
  ident.name = "ident"
) {
  scData <- RunPCA(object = scData, verbose = FALSE)
  #scData <- RunPCA(object=scData, npcs = param$npcs, verbose=FALSE)
  if (!('Spatial' %in% as.vector(Seurat::Assays(scData)))) {
    scData <- RunTSNE(
      object = scData,
      check_duplicates = FALSE,
      reduction = reduction,
      dims = 1:param$npcs
    )
  }
  scData <- RunUMAP(object = scData, reduction = reduction, dims = 1:param$npcs)
  scData <- FindNeighbors(
    object = scData,
    reduction = reduction,
    dims = 1:param$npcs,
    verbose = FALSE
  )
  # param$resolution can arrive as: a numeric scalar (web UI), a numeric vector
  # (param-set TSV), or a stringified Ruby array like "[0.6, 0.2, 0.4, 0.6, 0.8, 1]"
  # (Ruby-API submissions). Normalize to a numeric vector; the FIRST parsed
  # value is the user's selected resolution (matches the SUSHI dropdown
  # convention where [default, ...other_options] sets the default).
  parsedResolution <- suppressWarnings(as.numeric(unlist(strsplit(
    gsub("\\[|\\]", "", as.character(param$resolution)),
    "[,[:space:]]+"
  ))))
  parsedResolution <- parsedResolution[
    !is.na(parsedResolution) & parsedResolution > 0
  ]
  if (length(parsedResolution) == 0) {
    parsedResolution <- 0.6 # documented default
  }
  userResolution <- parsedResolution[1]

  myResolutions <- sort(unique(round(
    c(seq(from = 0.2, to = 1, by = 0.2), parsedResolution),
    1
  )))
  scData <- FindClusters(
    object = scData,
    resolution = myResolutions,
    verbose = FALSE
  ) #calculate clusters for a set of resolutions
  selectedCol <- paste0(
    DefaultAssay(scData),
    "_snn_res.",
    userResolution
  )
  if (!selectedCol %in% colnames(scData@meta.data)) {
    # Fallback: pick whatever <DefaultAssay>_snn_res.<r> column FindClusters
    # actually wrote. Defensive guard against surprise DefaultAssay values;
    # gives a clear error instead of "undefined columns selected".
    candidates <- grep(
      paste0("^", DefaultAssay(scData), "_snn_res\\."),
      colnames(scData@meta.data),
      value = TRUE
    )
    if (length(candidates) == 0) {
      stop(
        "seuratStandardWorkflow: no <assay>_snn_res.<r> columns produced; ",
        "DefaultAssay = '",
        DefaultAssay(scData),
        "' but no matching graph in meta.data."
      )
    }
    selectedCol <- candidates[1]
  }
  Idents(scData) <- scData@meta.data[, selectedCol] #but keep as the current clusters the ones obtained with the resolution set by the user
  scData[[ident.name]] <- Idents(scData)
  return(scData)
}

seuratClusteringV3 <- function(scData, param, assay = "RNA") {
  vars.to.regress <- getSeuratVarsToRegress(param)

  scData <- SCTransform(
    scData,
    vars.to.regress = vars.to.regress,
    assay = assay,
    seed.use = 38,
    verbose = FALSE
  )
  scData <- seuratStandardWorkflow(scData, param)
  return(scData)
}

seuratClusteringHTO <- function(scData) {
  DefaultAssay(scData) <- "HTO"
  scData <- ScaleData(scData)
  scData <- RunPCA(
    scData,
    features = rownames(scData),
    reduction.name = "pca_hto",
    reduction.key = "pca_hto_",
    verbose = FALSE
  )
  # Now, we rerun tSNE using the PCA only on ADT (protein) levels.
  scData <- RunTSNE(
    scData,
    reduction = "pca_hto",
    reduction.key = "htoTSNE_",
    reduction.name = "tsne_hto",
    check_duplicates = FALSE
  )
  scData <- FindNeighbors(
    scData,
    reduction = "pca_hto",
    features = rownames(scData),
    dims = NULL
  )
  scData <- FindClusters(scData, resolution = 0.2)
  return(scData)
}

##' @title Which normalization the multi-sample apps should use
##' @description \code{param$normalizationMethod} is "SCTransform" (the historical
##' behaviour, and the default when the parameter is absent) or "LogNormalize".
seuratNormalizationMethod <- function(param) {
  if (ezIsSpecified(param$normalizationMethod)) {
    param$normalizationMethod
  } else {
    "SCTransform"
  }
}

##' @title The assay PCA, clustering and markers run on
##' @description SCTransform writes its own "SCT" assay; log-normalization keeps
##' working in the raw assay ("RNA", or "Spatial" for the slides app).
seuratAnalysisAssay <- function(param, rawAssay = "RNA") {
  if (seuratNormalizationMethod(param) == "LogNormalize") rawAssay else "SCT"
}

##' @title Normalize each sample separately
##' @description Per-sample normalization, as both integration paths need. Note
##' log-normalization is a per-CELL operation (counts / that cell's total * 1e4),
##' so doing it per sample and doing it on the merged object give identical
##' values - which is what makes the later \code{JoinLayers} lossless.
##' @return the list of Seurat objects, each normalized with variable features set.
seuratNormalizeSampleList <- function(scDataList, param, assay = "RNA") {
  vars.to.regress <- getSeuratVarsToRegress(param)
  logNorm <- seuratNormalizationMethod(param) == "LogNormalize"
  for (i in seq_along(scDataList)) {
    if (logNorm) {
      scDataList[[i]] <- NormalizeData(
        scDataList[[i]],
        normalization.method = "LogNormalize",
        scale.factor = 10000,
        assay = assay,
        verbose = FALSE
      )
      scDataList[[i]] <- FindVariableFeatures(
        scDataList[[i]],
        selection.method = "vst",
        nfeatures = param$nfeatures,
        verbose = FALSE
      )
    } else {
      scDataList[[i]] <- SCTransform(
        scDataList[[i]],
        vars.to.regress = vars.to.regress,
        assay = assay,
        verbose = TRUE
      )
    }
  }
  return(scDataList)
}

##' @title Scale a merged log-normalized object so RunPCA and DoHeatmap have data
##' @description SCTransform produces its own scale.data; the log-normalization
##' path has to do it explicitly, on the joined layers, after the variable
##' features are known. Cell-cycle regression happens here rather than in SCT.
seuratScaleMergedLogNorm <- function(scData, param, assay, features) {
  DefaultAssay(scData) <- assay
  scData <- JoinLayers(scData, assay = assay)
  VariableFeatures(scData) <- unique(features)
  scData <- ScaleData(
    scData,
    vars.to.regress = getSeuratVarsToRegress(param),
    verbose = FALSE
  )
  return(scData)
}

cellClustNoCorrection <- function(scDataList, param) {
  if (param[['name']] == 'SpatialSeuratSlides') {
    assay = "Spatial"
  } else {
    assay = "RNA"
  }
  #1. Data preprocesing is done on each object separately
  scDataList <- seuratNormalizeSampleList(scDataList, param, assay = assay)
  #2. Merge all seurat objects
  scData = merge(x = scDataList[[1]], y = scDataList[-1], merge.data = TRUE)
  allVariableFeatures <- unlist(lapply(scDataList, VariableFeatures))
  if (seuratNormalizationMethod(param) == "LogNormalize") {
    scData <- seuratScaleMergedLogNorm(
      scData,
      param,
      assay = assay,
      features = allVariableFeatures
    )
  } else {
    VariableFeatures(scData) <- allVariableFeatures
  }
  #3. Dimensionality reduction and clustering
  scData <- seuratStandardWorkflow(scData, param)

  return(scData)
}

cellClustWithCorrection <- function(scDataList, param) {
  if (param[['name']] == 'SpatialSeuratSlides') {
    assay = "Spatial"
  } else {
    assay = "RNA"
  }
  logNorm <- seuratNormalizationMethod(param) == "LogNormalize"
  ## Seurat's own name for the normalization, as the anchor/integration API wants it
  anchorNormMethod <- if (logNorm) "LogNormalize" else "SCT"
  #1. Data preprocesing is done on each object separately
  scDataList <- seuratNormalizeSampleList(scDataList, param, assay = assay)

  #2. Data integration
  #2.1. # Select the most variable features to use for integration
  integ_features <- SelectIntegrationFeatures(
    object.list = scDataList,
    nfeatures = param$nfeatures
  )

  if (param$integrationMethod %in% c("RPCA", "CCA")) {
    for (i in 1:length(scDataList)) {
      scDataList[[i]] <- ScaleData(scDataList[[i]], features = integ_features)
    }
    #2.2. Prepare the SCT list object for integration (SCT residuals only)
    if (!logNorm) {
      scDataList <- PrepSCTIntegration(
        object.list = scDataList,
        anchor.features = integ_features
      )
    }
    if (param$integrationMethod == 'RPCA') {
      scDataList <- lapply(
        X = scDataList,
        FUN = RunPCA,
        features = integ_features
      )
    }
    #2.3. Find anchors
    if (param$integrationMethod == 'CCA') {
      integ_anchors <- FindIntegrationAnchors(
        object.list = scDataList,
        normalization.method = anchorNormMethod,
        anchor.features = integ_features,
        dims = 1:param$npcs
      )
    } else if (param$integrationMethod == 'RPCA') {
      integ_anchors <- FindIntegrationAnchors(
        object.list = scDataList,
        normalization.method = anchorNormMethod,
        anchor.features = integ_features,
        dims = 1:param$npcs,
        reduction = "rpca",
        k.anchor = 20
      )
    }

    #2.4. Integrate datasets
    seurat_integrated <- IntegrateData(
      anchorset = integ_anchors,
      normalization.method = anchorNormMethod,
      dims = 1:param$npcs
    )

    #3. Run the standard workflow for visualization and clustering
    seurat_integrated <- seuratStandardWorkflow(seurat_integrated, param)
  } else if (param$integrationMethod == "Harmony") {
    require(harmony)
    #2.2 Merge normalized samples
    scData <- merge(
      x = scDataList[[1]],
      y = scDataList[2:length(scDataList)],
      merge.data = TRUE
    )
    #2.3.1 Manually set variable features of merged Seurat object
    if (logNorm) {
      scData <- seuratScaleMergedLogNorm(
        scData,
        param,
        assay = assay,
        features = integ_features
      )
    } else {
      VariableFeatures(scData) <- integ_features
    }
    #2.3.2 Calculate PCs using manually set variable features
    scData <- RunPCA(
      scData,
      assay = seuratAnalysisAssay(param, rawAssay = assay),
      npcs = param$npcs
    )

    #2.4 Prep and run Harmony algorithm
    # param$harmonyGroupBy (ScSeuratCombine knob, default "Condition" for backward
    # compatibility) names the metadata column(s) Harmony corrects for; "Batch" is
    # one level per input sample. Any additionalFactors are corrected for as well.
    harmonyFactors <- if (ezIsSpecified(param$harmonyGroupBy)) {
      param$harmonyGroupBy
    } else {
      "Condition"
    }
    if (!is.null(param$additionalFactors)) {
      harmonyFactors <- c(
        harmonyFactors,
        colnames(scData@meta.data)[startsWith(colnames(scData@meta.data), "meta_")]
      )
    }
    if ("Batch" %in% harmonyFactors && !"Batch" %in% colnames(scData@meta.data)) {
      scData$Batch <- scData$Sample # same fallback as the report
    }
    missingFactors <- setdiff(harmonyFactors, colnames(scData@meta.data))
    if (length(missingFactors) > 0) {
      stop("harmonyGroupBy column(s) not in the cell metadata: ",
           paste(missingFactors, collapse = ", "))
    }
    # Calculate Harmony reduction
    # harmony >= 1.0 renamed `reduction` -> `reduction.use` and removed `assay.use`
    # (check_legacy_args() now hard-errors on assay.use); the assay is taken from the reduction.
    scData <- RunHarmony(
      scData,
      group.by.vars = harmonyFactors,
      reduction.use = "pca",
      reduction.save = "harmony"
    )

    #3. Run the standard workflow for visualization and clustering
    seurat_integrated <- seuratStandardWorkflow(
      scData,
      param,
      reduction = "harmony"
    )
  }

  return(seurat_integrated)
}

posClusterMarkers <- function(scData, pvalue_allMarkers, param) {
  vars.to.regress = NULL
  if (param$name == "SCOneSample" & param$DE.method == "LR") {
    #For the one sample app, we will regress out cell cycle if the test is LR
    vars.to.regress <- c("CellCycleS", "CellCycleG2M")
  } else if (param$name == "SCReportMerging" & param$DE.method == "LR") {
    #for multiple samples we will regress out either the cell cycle, plate or both if the test is LR
    vars.to.regress <- param$DE.regress
  }

  markers <- FindAllMarkers(
    object = scData,
    test.use = param$DE.method,
    only.pos = TRUE,
    latent.vars = vars.to.regress,
    return.thresh = pvalue_allMarkers
  )
  ## Significant markers
  cm <- markers[, c(
    "gene",
    "cluster",
    "pct.1",
    "pct.2",
    "avg_log2FC",
    "p_val_adj"
  )]
  cm$cluster <- as.factor(cm$cluster)
  diff_pct = abs(cm$pct.1 - cm$pct.2)
  cm$diff_pct <- diff_pct
  cm <- cm[order(cm$diff_pct, decreasing = TRUE), ] %>%
    mutate_if(is.numeric, round, digits = 3)
  rownames(cm) <- NULL
  return(cm)
}


SpatiallyVariableFeatures_workaround <- function(
  object,
  assay = "SCT",
  selection.method = "moransi"
) {
  #' This is work around function to replace SeuratObject::SpatiallyVariableFeatures function.
  #' return ranked list of Spatially Variable Features

  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }

  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features

  # Select columns starting with the provided col_prefix
  moransi_cols <- grep(
    paste0("^", selection.method),
    colnames(data),
    value = TRUE
  )

  # Filter rows where "moransi.spatially.variable" is TRUE
  filtered_data <- data[
    data[[paste0(selection.method, ".spatially.variable")]],
    moransi_cols
  ]

  # Sort filtered data by "moransi.spatially.variable.rank" column in ascending order
  sorted_data <- filtered_data[
    order(filtered_data[[paste0(
      selection.method,
      ".spatially.variable.rank"
    )]]),
  ]
  sorted_data <- sorted_data[
    grep('^NA', rownames(sorted_data), invert = TRUE),
  ]
  sorted_data[['Rank']] = 1:nrow(sorted_data)
  # Return row names of the sorted data frame
  return(sorted_data)
}

spatialMarkers <- function(scData, selection.method = "markvariogram") {
  scData <- FindSpatiallyVariableFeatures(
    scData,
    features = VariableFeatures(scData),
    r.metric = 5,
    verbose = TRUE,
    selection.method = selection.method
  )
  spatialMarkers <- SpatiallyVariableFeatures(
    scData,
    selection.method = selection.method
  )
  #spatialMarkers <- SpatiallyVariableFeatures_workaround(scData, assay="SCT", selection.method = selection.method)
  return(spatialMarkers)
}

all2all <- function(scData, pvalue_all2allMarkers, param) {
  clusterCombs <- combn(levels(Idents(scData)), m = 2)
  all2allMarkers <- mcmapply(
    FindMarkers,
    as.integer(clusterCombs[1, ]),
    as.integer(clusterCombs[2, ]),
    MoreArgs = list(object = scData, only.pos = FALSE),
    mc.preschedule = FALSE,
    mc.cores = min(4L, param$cores),
    SIMPLIFY = FALSE
  )
  all2allMarkers <- lapply(all2allMarkers, function(x) {
    x[x$p_val <= pvalue_all2allMarkers, ]
  })
  names(all2allMarkers) <- apply(clusterCombs, 2, paste, collapse = "vs")
  return(all2allMarkers)
}

conservedMarkers <- function(
  scData,
  grouping.var = "Condition",
  pseudoBulkMode = FALSE
) {
  markers <- list()

  if ("SCT" %in% Seurat::Assays(scData)) {
    assay <- "SCT"
  } else {
    assay <- "RNA"
  }
  if (pseudoBulkMode) {
    DE.method <- "DESeq2"
  } else {
    DE.method <- "wilcox"
  }

  # Track if we fall back to RNA due to SCT errors
  fallback_to_rna <- FALSE

  for (eachCluster in levels(Idents(scData))) {
    markersEach <- try(
      FindConservedMarkers(
        scData,
        ident.1 = eachCluster,
        grouping.var = grouping.var,
        print.bar = FALSE,
        only.pos = TRUE,
        assay = assay,
        test.use = DE.method
      ),
      silent = TRUE
    )

    # If SCT fails with PrepSCTFindMarkers error, retry with RNA
    if (class(markersEach) == "try-error" && assay == "SCT") {
      error_msg <- as.character(markersEach)
      if (grepl("PrepSCTFindMarkers|unequal library sizes", error_msg)) {
        if (!fallback_to_rna) {
          ezLog(
            "SCT assay failed with PrepSCTFindMarkers error. Falling back to RNA assay."
          )
          fallback_to_rna <- TRUE
        }
        markersEach <- try(
          FindConservedMarkers(
            scData,
            ident.1 = eachCluster,
            grouping.var = grouping.var,
            print.bar = FALSE,
            only.pos = TRUE,
            assay = "RNA",
            test.use = DE.method
          ),
          silent = TRUE
        )
      }
    }

    if (class(markersEach) != "try-error" && nrow(markersEach) > 0) {
      markers[[eachCluster]] <- as_tibble(markersEach, rownames = "gene")
    }
  }
  markers <- bind_rows(markers, .id = "cluster")
  markers$avg_avg_fc <- markers %>%
    dplyr::select(contains("_avg_log2FC")) %>%
    rowMeans()
  return(markers)
}

diffExpressedGenes <- function(scData, param, grouping.var = "Condition") {
  seurat_clusters <- Idents(scData)
  scData@meta.data$cluster.condition <-
    paste0(seurat_clusters, "_", scData@meta.data[[grouping.var]])
  Idents(scData) <- "cluster.condition"
  conditions <- unique(scData@meta.data[[grouping.var]])

  vars.to.regress = NULL
  DE.method <- param$DE.method
  if (
    ezIsSpecified(param$replicateGrouping) && param$pseudoBulkMode == "true"
  ) {
    DE.method <- "DESeq2"
  }
  if (param$DE.method == "LR") {
    #regress the plate if the test is LR
    vars.to.regress <- param$DE.regress
  }

  # Determine assay to use
  if ("SCT" %in% Seurat::Assays(scData)) {
    assay <- "SCT"
  } else {
    assay <- "RNA"
  }

  # Track if we fall back to RNA due to SCT errors
  fallback_to_rna <- FALSE

  diffGenes <- list()
  for (eachCluster in gtools::mixedsort(levels(seurat_clusters))) {
    markersEach <- try(FindMarkers(
      scData,
      ident.1 = paste0(eachCluster, "_", param$sampleGroup),
      ident.2 = paste0(eachCluster, "_", param$refGroup),
      test.use = DE.method,
      latent.vars = vars.to.regress,
      assay = assay
    ))

    # If SCT fails with PrepSCTFindMarkers error, retry with RNA
    if (class(markersEach) == "try-error" && assay == "SCT") {
      error_msg <- as.character(markersEach)
      if (grepl("PrepSCTFindMarkers|unequal library sizes", error_msg)) {
        if (!fallback_to_rna) {
          ezLog(
            "SCT assay failed with PrepSCTFindMarkers error. Falling back to RNA assay for DEG analysis."
          )
          fallback_to_rna <- TRUE
        }
        markersEach <- try(FindMarkers(
          scData,
          ident.1 = paste0(eachCluster, "_", param$sampleGroup),
          ident.2 = paste0(eachCluster, "_", param$refGroup),
          test.use = DE.method,
          latent.vars = vars.to.regress,
          assay = "RNA"
        ))
      }
    }

    ## to skip some groups with few cells
    if (class(markersEach) != "try-error") {
      diffGenes[[eachCluster]] <- as_tibble(markersEach, rownames = "gene")
    }
  }
  diffGenes <- bind_rows(diffGenes, .id = "cluster")

  diffGenes$diff_pct = abs(diffGenes$pct.1 - diffGenes$pct.2)
  rownames(diffGenes) <- NULL
  # diffGenes <- diffGenes[diffGenes$p_val_adj < 0.05, ]

  return(diffGenes)
}

getSeuratVarsToRegress <- function(param) {
  vars.to.regress <- NULL
  if (
    ezIsSpecified(param$SCT.regress.CellCycle) &&
      param$SCT.regress.CellCycle
  ) {
    vars.to.regress <- c("CC.Difference")
  }
  if (!is.null(param$SCT.regress.var)) {
    vars.to.regress <- c(vars.to.regress, param$SCT.regress.var)
  }
  return(vars.to.regress)
}

getSeuratMarkers <- function(scData, param) {
  # positive cluster markers
  ## https://github.com/satijalab/seurat/issues/5321
  ## https://github.com/satijalab/seurat/issues/1501
  markers <- FindAllMarkers(
    object = scData,
    test.use = param$DE.method,
    only.pos = TRUE,
    min.pct = ifelse(ezIsSpecified(param$min.pct), param$min.pct, 0.1),
    logfc.threshold = ifelse(
      ezIsSpecified(param$logfc.threshold),
      param$logfc.threshold,
      0.25
    )
  )
  ## Significant markers
  markers <- markers[, c(
    "gene",
    "cluster",
    "pct.1",
    "pct.2",
    "avg_log2FC",
    "p_val_adj"
  )]
  markers <- markers[markers$p_val_adj < param$pvalue_allMarkers, ]
  markers$cluster <- as.factor(markers$cluster)
  markers$diff_pct = abs(markers$pct.1 - markers$pct.2)
  markers <- markers[markers$diff_pct >= param$min.diff.pct, ]
  markers <- markers[order(markers$diff_pct, decreasing = TRUE), ]
  return(markers)
}

getSeuratMarkersAndAnnotate <- function(scData, param, BPPARAM) {
  # function for general annotation of Seurat objects

  # Helper function to strip species labels from parameter values
  strip_species_labels <- function(param_values) {
    if (is.character(param_values)) {
      # Remove species labels like " (human)" or " (mouse)"
      return(gsub(" \\(human\\)| \\(mouse\\)", "", param_values))
    }
    return(param_values)
  }

  # Clean Azimuth and SingleR parameters to remove species labels
  if (ezIsSpecified(param$Azimuth)) {
    param$Azimuth <- strip_species_labels(param$Azimuth)
  }
  if (ezIsSpecified(param$SingleR)) {
    param$SingleR <- strip_species_labels(param$SingleR)
  }

  markers <- getSeuratMarkers(scData, param)

  # cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if (species == "Human" | species == "Mouse") {
    genesPerCluster <- split(markers$gene, markers$cluster)
    enrichRout <- querySignificantClusterAnnotationEnrichR(
      genesPerCluster,
      param$enrichrDatabase
    )
    cells.AUC <- cellsLabelsWithAUC(
      GetAssayData(scData, layer = "counts"),
      species,
      param$tissue,
      BPPARAM = BPPARAM
    )
    singler.results <- cellsLabelsWithSingleR(
      GetAssayData(scData, layer = "data"),
      Idents(scData),
      param$SingleR,
      BPPARAM = BPPARAM
    )
    for (r in names(singler.results)) {
      scData[[paste0(r, "_single")]] <- singler.results[[r]]$single.fine$labels
      scData[[paste0(r, "_cluster")]] <- singler.results[[
        r
      ]]$cluster.fine$labels[match(
        Idents(scData),
        rownames(singler.results[[r]]$cluster.fine)
      )]
    }
  } else {
    cells.AUC <- NULL
    singler.results <- NULL
    enrichRout <- NULL
  }

  ## SCpubr advanced plots, can currently only be computed for human and mouse
  if (
    ezIsSpecified(param$computePathwayTFActivity) &&
      as.logical(param$computePathwayTFActivity) &&
      (species == "Human" | species == "Mouse")
  ) {
    pathwayActivity <- computePathwayActivityAnalysis(
      cells = scData,
      species = species
    )
    TFActivity <- computeTFActivityAnalysis(cells = scData, species = species)
  } else {
    pathwayActivity <- NULL
    TFActivity <- NULL
    futile.logger::flog.info("Skipping pathway and TF activity")
  }

  # run Azimuth
  if (ezIsSpecified(param$Azimuth) && param$Azimuth != "none") {
    if (!requireNamespace("Azimuth", quietly = TRUE)) {
      stop("Azimuth annotation requested, but package 'Azimuth' is not available.")
    }
    environment(MyDietSeurat) <- asNamespace('Seurat')
    assignInNamespace("DietSeurat", MyDietSeurat, ns = "Seurat")
    scData[["RNA"]] <- JoinLayers(scData[["RNA"]]) # Required for Azimuth compatibility with Seurat v5
    scDataAzi <- Azimuth::RunAzimuth(scData, param$Azimuth, assay = "RNA") ## TODO support ADT

    ##Rename annotion levels if neccessary:
    colnames(scDataAzi@meta.data) <- sub(
      'level_',
      'l',
      colnames(scDataAzi@meta.data)
    )

    aziNames <- setdiff(
      colnames(scDataAzi@meta.data),
      colnames(scData@meta.data)
    )
    aziResults <- data.frame(
      Azimuth.celltype.l1 = scDataAzi@meta.data[, grep(
        "l1$",
        aziNames,
        value = TRUE
      )],
      Azimuth.celltype.l2 = scDataAzi@meta.data[, grep(
        "l2$",
        aziNames,
        value = TRUE
      )],
      Azimuth.celltype.l3 = scDataAzi@meta.data[, grep(
        "l3$",
        aziNames,
        value = TRUE
      )],
      Azimuth.celltype.l4 = scDataAzi@meta.data[, grep(
        "l4$",
        aziNames,
        value = TRUE
      )],
      row.names = colnames(scDataAzi)
    )
    ## TODO: score should also be stored
    remove(scDataAzi)
  } else {
    aziResults <- NULL
  }

  # run cellxgene_annotation
  if (
    ezIsSpecified(param$cellxgeneUrl) && ezIsSpecified(param$cellxgeneLabel)
  ) {
    if(!requireNamespace("qs", quietly = TRUE)) {
      warning("cellxgene annotation requested, but package 'qs' is not available. Skip annotation.")
        cellxgeneResults <- NULL
    } else {
        cellxgeneResults <- cellxgene_annotation(scData = scData, param = param)
    }
  } else {
    cellxgeneResults <- NULL
  }

  return(list(
    markers = markers,
    cells.AUC = cells.AUC,
    singler.results = singler.results,
    enrichRout = enrichRout,
    pathwayActivity = pathwayActivity,
    TFActivity = TFActivity,
    aziResults = aziResults,
    cellxgeneResults = cellxgeneResults
  ))
}

##' @title Merge Seurat clusters
##' @description Given an input mapping, group Seurat clusters of the active Seurat Ident together, i.e. merging them
##' @slot scData the Seurat object
##' @slot clustList a list() object where the names are the new cluster names and their values are all the clusters to be merged. If the list consists of multiple sets of clusters to be merged, ensure the sets do not overlap.
##' @template roxygen-template
##' @examples
##' data("pbmc_small")
##' clustList <- list("foo"=as.character(c(1, 2)))
##' seuratMergeClusters(pbmc_small, clustList)
seuratMergeClusters <- function(scData, clustList) {
  require(stringr)
  assertthat::assert_that(
    length(clustList) <= 1 || length(Reduce(intersect, clustList)) == 0,
    msg = "Sets of input clusters overlap!"
  )

  clusters <- as.character(unique(Idents(scData)))
  newIdent <- as.character(unname(Idents(scData)))
  for (targetClust in names(clustList)) {
    clustersToChange <- clustList[[targetClust]]
    targetClustRep <- rep(targetClust, length(clustersToChange))
    names(targetClustRep) <- clustersToChange
    newIdent <- ifelse(
      newIdent %in% clustersToChange,
      unname(targetClustRep[clustersToChange]),
      newIdent
    )
  }
  newIdent <- factor(
    newIdent,
    levels = str_sort(unique(newIdent), numeric = TRUE)
  )
  return(newIdent)
}

Seurat_to_csv <- function(seurat_obj, export_dir) {
  if (!file.exists(export_dir)) {
    dir.create(file.path(export_dir))
  }

  # Gene Counts Matrix
  counts_mat <- as.matrix(GetAssayData(
    seurat_obj,
    assay = DefaultAssay(seurat_obj)
  )) %>%
    t()
  write.table(
    counts_mat,
    paste(export_dir, "/X.csv", sep = ""),
    row.names = FALSE,
    col.names = FALSE,
    sep = ","
  )

  #Metadata Matrix
  obs_mat <- seurat_obj[[]]
  write.table(obs_mat, paste(export_dir, "/obs.csv", sep = ""), sep = ",")

  # Feature Metadata
  var_mat <- matrix(TRUE, nrow = length(Features(seurat_obj)), ncol = 1)
  rownames(var_mat) <- Features(seurat_obj)
  write.table(var_mat, paste(export_dir, "/var.csv", sep = ""), sep = ",")

  # Reductions and Spatial Coordinates
  X_pca <- seurat_obj@reductions$pca@cell.embeddings
  colnames(X_pca) <- sapply(1:ncol(X_pca), function(x) {
    paste("X_pca", x, sep = "")
  })
  X_umap <- seurat_obj@reductions$umap@cell.embeddings
  colnames(X_umap) <- sapply(1:ncol(X_umap), function(x) {
    paste("X_umap", x, sep = "")
  })

  spatial_list <- lapply(Images(seurat_obj), function(x) {
    seurat_obj@images[[x]]$centroids@coords
  })
  spatial <- do.call("rbind", spatial_list)
  colnames(spatial) <- sapply(1:ncol(spatial), function(x) {
    paste("spatial", x, sep = "")
  })

  obsm_mat <- cbind(X_pca, X_umap, spatial)
  write.table(
    obsm_mat,
    paste(export_dir, "/obsm.csv", sep = ""),
    row.names = FALSE,
    sep = ","
  )
}


csv_to_Seurat <- function(import_dir, assay_name) {
  counts_mat <- read.table(
    paste(import_dir, "/X.csv", sep = ""),
    header = FALSE,
    sep = ",",
    fill = TRUE
  ) %>%
    as.matrix() %>%
    t()
  obs_mat <- read.table(
    paste(import_dir, "/obs.csv", sep = ""),
    header = TRUE,
    sep = ",",
    row.names = 1,
    fill = TRUE
  )
  var_mat <- read.table(
    paste(import_dir, "/var.csv", sep = ""),
    header = TRUE,
    sep = ",",
    row.names = 1,
    fill = TRUE
  )
  obsm_mat <- read.table(
    paste(import_dir, "/obsm.csv", sep = ""),
    header = TRUE,
    sep = ",",
    fill = TRUE
  )

  rownames(counts_mat) <- rownames(var_mat)
  colnames(counts_mat) <- rownames(obs_mat)

  X_pca <- select(obsm_mat, starts_with('X_pca'))
  colnames(X_pca) <- sapply(1:ncol(X_pca), function(x) {
    paste("pca", x, sep = "_")
  })
  rownames(X_pca) <- rownames(obs_mat)

  X_umap <- select(obsm_mat, starts_with('X_umap'))
  colnames(X_umap) <- sapply(1:ncol(X_umap), function(x) {
    paste("umap", x, sep = "_")
  })
  rownames(X_umap) <- rownames(obs_mat)

  spatial <- select(obsm_mat, starts_with("spatial"))
  colnames(spatial) <- sapply(1:ncol(spatial), function(x) {
    paste("spatial", x, sep = "_")
  })
  rownames(spatial) <- rownames(obs_mat)

  seurat_import <- CreateSeuratObject(counts = counts_mat, assay = assay_name)
  seurat_import <- AddMetaData(seurat_import, obs_mat)
  seurat_import[["pca"]] <- CreateDimReducObject(
    embeddings = as.matrix(X_pca),
    key = "pca_",
    assay = DefaultAssay(seurat_import)
  )
  seurat_import[["umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(X_umap),
    key = "umap_",
    assay = DefaultAssay(seurat_import)
  )
  seurat_import[["spatial"]] <- CreateFOV(
    CreateCentroids(spatial[, c(1, 2)]),
    assay = assay_name
  )

  return(seurat_import)
}

Seurat_to_SPE <- function(seurat_obj, coords) {
  spE <- SpatialExperiment(
    assay = GetAssayData(seurat_obj),
    colData = seurat_obj[[]],
    spatialCoords = coords
  )
  return(spE)
}

Seurat_to_SPIATSPE <- function(seurat_obj) {
  ## requires SPIAT package
  spE <- format_image_to_spe(
    format = "general",
    intensity_matrix = as.matrix(GetAssayData(seurat_obj)),
    phenotypes = as.character(Idents(seurat_obj)),
    coord_x = seurat_obj@images$image$centroids@coords[, 1],
    coord_y = seurat_obj@images$image$centroids@coords[, 2]
  )
  return(spE)
}


## does the following
## if the image has 4 channels assume, the 4th channel is alpha and blend in
## if the image has more than 4 channels take only the first 3
## if the blue channel is brighter than the red channel swap these two channels
fix_microscopy_image <- function(
  obj,
  img_name = Images(obj)[1],
  eps = 0.02,
  margin = 0.01
) {
  arr <- obj@images[[img_name]]@image

  # Normalize if stored as 8/16-bit (no-op if already 0..1)
  mx <- suppressWarnings(max(arr, na.rm = TRUE))
  if (!is.finite(mx) || mx == 0) {
    mx <- 1
  }
  if (mx > 1) {
    arr <- arr / mx
  }

  if (length(dim(arr)) == 3) {
    if (dim(arr)[3] == 4) {
      # if 4 channels assume RGBA, do composite over white and drop alpha
      A <- arr[,, 4]
      RGB <- arr[,, 1:3]
      for (k in 1:3) {
        RGB[,, k] <- RGB[,, k] * A + (1 - A)
      }
      arr <- RGB
    }
    if (dim(arr)[3] > 4) {
      ## simply take the first 3
      arr <- arr[,, 1:3, drop = FALSE]
    }
  }

  # mask out near-black background; keep pixels with any channel > eps
  mask <- apply(arr > eps, c(1, 2), any)
  if (!any(mask)) {
    mask <- matrix(TRUE, nrow = dim(arr)[1], ncol = dim(arr)[2])
  }

  mR <- mean(arr[,, 1][mask], na.rm = TRUE)
  mB <- mean(arr[,, 3][mask], na.rm = TRUE)

  if (is.finite(mR) && is.finite(mB) && (mB > mR + margin)) {
    arr <- arr[,, c(3, 2, 1), drop = FALSE]
    chosen <- "BGR->RGB (swapped)"
  } else {
    chosen <- "RGB (kept as-is)"
  }

  obj@images[[img_name]]@image <- arr
  ezLog(sprintf(
    "Channel order decided by means on tissue: R=%.3f, B=%.3f -> %s",
    mR,
    mB,
    chosen
  ))
  obj
}


ezSpatialFeaturePlot <- function(
  obj,
  feature,
  title = feature,
  label = feature,
  pt.size.factor = 2,
  plot_segmentations = TRUE,
  min.cutoff = "q01",
  max.cutoff = "q99",
  legend.pos = "right",
  palette = "turbo",
  font_size = 10, ## that's the default fontsize from plotColData
  ...
) {
  SpatialFeaturePlot(
    object = obj,
    features = feature,
    images = Images(obj),
    plot_segmentations = plot_segmentations,
    pt.size.factor = pt.size.factor,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    ...
  ) +
    labs(title = title, fill = label) +
    scale_fill_viridis_c(option = palette, na.value = "grey95") +
    theme(legend.position = legend.pos) +
    cowplot::theme_cowplot(font_size = font_size)
  # theme_void() +
  # theme(
  #   legend.position   = legend.pos,
  #   plot.title        = element_text(size = title_size,
  #                                    face = "bold",
  #                                    margin = margin(b = 6)),
  #   legend.title      = element_text(size = 8),
  #   legend.text       = element_text(size = 8),
  #   legend.key.height = unit(10, "pt"),
  #   legend.key.width  = unit(6, "pt")
  # )
}

##' @title Marker heatmap that survives large cell counts
##' @description \code{Seurat::DoHeatmap} draws one raster column per cell,
##' twice: the expression body (\code{raster = TRUE}) and the group bar (always).
##' It also pads every group with \code{ceiling(nCells * 0.0025)} blank columns.
##' Cairo refuses a raster wider than 32767 px and then silently stops drawing
##' for the rest of the page. Above that width the report got either an
##' all-white heatmap (body too wide; p39836/o43013, 2026-09-02) or, with
##' \code{raster = FALSE}, tiles without gene names, group bar, cluster labels
##' or legend (group bar too wide; p39836, 2026-09-03: 30000 cells + 46 groups
##' = 33450 columns). This wrapper caps the number of cells drawn so the widest
##' raster stays below the limit; the drawn cells are a random sample.
##' @param object a Seurat object.
##' @param features features to plot.
##' @param maxRasterWidth largest raster width handed to the device, in columns.
##' @param ... passed to \code{Seurat::DoHeatmap}; an explicit \code{cells} is
##'   subsampled if it is still too wide.
##' @return the ggplot returned by \code{Seurat::DoHeatmap}.
ezDoHeatmap <- function(object, features, maxRasterWidth = 31000, ...) {
  args <- list(...)
  cells <- if (is.null(args$cells)) colnames(object) else args$cells
  nGroups <- length(unique(Seurat::Idents(object)[cells]))
  ## width = nCells + nGroups * ceiling(nCells * 0.0025); the ceiling adds at most nGroups
  maxCells <- floor((maxRasterWidth - nGroups) / (1 + 0.0025 * nGroups))
  if (length(cells) > maxCells) {
    args$cells <- sample(cells, maxCells)
  }
  do.call(Seurat::DoHeatmap, c(list(object = object, features = features), args))
}
