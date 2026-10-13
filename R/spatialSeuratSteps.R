###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## Analysis steps of the SpatialSeurat app (Visium SD, one sample per job).
## Every step takes the Seurat object and returns a list with the (updated)
## object and step results. ezMethodSpatialSeurat() chains them; optional
## steps run through run_optional_step() so that one failing add-on does not
## kill a course job.

##' @title Run an optional analysis step
##' @param name step name, used in the log and in the report.
##' @param expr expression returning the step result.
##' @param status environment collecting the step status.
##' @return the value of \code{expr}, or NULL when it failed.
run_optional_step <- function(name, expr, status) {
  t0 <- Sys.time()
  ezLog("Step ", name, " started")
  res <- tryCatch(
    {
      value <- force(expr)
      status[[name]] <- list(ok = TRUE, msg = "")
      value
    },
    error = function(e) {
      ezLog("Step ", name, " FAILED: ", conditionMessage(e))
      status[[name]] <- list(ok = FALSE, msg = conditionMessage(e))
      NULL
    }
  )
  status[[name]]$minutes <- round(
    as.numeric(difftime(Sys.time(), t0, units = "mins")),
    2
  )
  return(res)
}

##' @title Mark a step as skipped by parameter choice
skip_step <- function(name, msg, status) {
  ezLog("Step ", name, " skipped: ", msg)
  status[[name]] <- list(ok = NA, msg = msg, minutes = 0)
  return(NULL)
}

##' @title Feature annotation of a Visium count matrix
##' @return data.frame with gene_id, gene_name, type, isMito, isRiboprot;
##'   row names are the feature names Seurat uses.
spatial_feature_info <- function(input) {
  cmDir <- input$getFullPaths("CountMatrix")
  featInfo <- ezRead.table(
    file.path(cmDir, "features.tsv.gz"),
    header = FALSE,
    row.names = NULL
  )
  colnames(featInfo) <- c("gene_id", "gene_name", "type")
  featInfo$isMito <- grepl("(?i)^MT-", featInfo$gene_name)
  featInfo$isRiboprot <- grepl("(?i)^RPS|^RPL", featInfo$gene_name)
  rownames(featInfo) <- gsub(
    "_",
    "-",
    scater::uniquifyFeatureNames(
      ID = featInfo$gene_id,
      names = featInfo$gene_name
    )
  )
  return(featInfo)
}

##' @title Read the Space Ranger tissue positions
##' @return data.frame with in_tissue, array_row, array_col,
##'   pxl_row_in_fullres, pxl_col_in_fullres; row names are barcodes.
spatial_tissue_positions <- function(spatialDir) {
  newFile <- file.path(spatialDir, "tissue_positions.csv")
  if (file.exists(newFile)) {
    pos <- utils::read.csv(newFile, header = TRUE)
  } else {
    pos <- utils::read.csv(
      file.path(spatialDir, "tissue_positions_list.csv"),
      header = FALSE
    )
    colnames(pos) <- c(
      "barcode",
      "in_tissue",
      "array_row",
      "array_col",
      "pxl_row_in_fullres",
      "pxl_col_in_fullres"
    )
  }
  rownames(pos) <- pos$barcode
  return(pos)
}

##' @title Load one Visium SD sample
##' @return list with scData, scDataRaw, param (imageEnlargementFactor and
##'   pt.size.factor filled in) and spotInfo (scale factors, chemistry).
spatial_load <- function(input, param) {
  res <- load10xSpatialData(input, param)
  scData <- res$scData
  param <- res$param
  spatialDir <- file.path(input$getFullPaths("ResultDir"), "spatial")
  scaleFactors <- jsonlite::fromJSON(file.path(
    spatialDir,
    "scalefactors_json.json"
  ))
  pos <- spatial_tissue_positions(spatialDir)
  barcodes <- sub(".*___", "", colnames(scData))
  scData$array_row <- pos[barcodes, "array_row"]
  scData$array_col <- pos[barcodes, "array_col"]
  scData$pxl_row_fullres <- pos[barcodes, "pxl_row_in_fullres"]
  scData$pxl_col_fullres <- pos[barcodes, "pxl_col_in_fullres"]
  scData$Sample <- input$getNames()

  if (!ezIsSpecified(param$pt.size.factor) || is.na(param$pt.size.factor)) {
    img <- scData@images[[1]]
    spotDiameter <- img@scale.factors$spot * img@scale.factors$hires
    param$pt.size.factor <- min(
      spotDiameter * 0.17,
      2 * param$imageEnlargementFactor
    )
  }
  ## Probe-based chemistry (FFPE / CytAssist) ships a probe_set.csv and has
  ## no mitochondrial probes.
  probeSetFile <- file.path(input$getFullPaths("ResultDir"), "probe_set.csv")
  spotInfo <- list(
    scaleFactors = scaleFactors,
    micronsPerPixel = 55 / scaleFactors$spot_diameter_fullres,
    probeBased = file.exists(probeSetFile),
    nSpotsTotal = nrow(pos),
    nSpotsTissue = sum(pos$in_tissue == 1)
  )
  return(list(
    scData = scData,
    scDataRaw = res$scDataRaw,
    param = param,
    spotInfo = spotInfo
  ))
}

##' @title Spot-level outliers with SpotSweeper
##' @description Local outliers compare each spot to its spatial neighbours,
##'   so biologically low-RNA regions such as fiber tracts are not removed as a
##'   whole, which a global threshold would do.
spatial_spot_sweeper <- function(scData, useMito = TRUE) {
  library(SpatialExperiment)
  library(SpotSweeper)
  coords <- as.matrix(data.frame(
    x = scData$pxl_col_fullres,
    y = scData$pxl_row_fullres,
    row.names = colnames(scData)
  ))
  spe <- SpatialExperiment(
    assays = list(counts = GetAssayData(scData, assay = "Spatial", layer = "counts")),
    colData = scData[[]],
    spatialCoords = coords
  )
  metrics <- c(nCount_Spatial = "lower", nFeature_Spatial = "lower")
  if (useMito) {
    metrics <- c(metrics, percent_mito = "higher")
  }
  for (m in names(metrics)) {
    spe <- localOutliers(
      spe,
      metric = m,
      direction = metrics[[m]],
      log = TRUE,
      samples = "Sample"
    )
  }
  res <- data.frame(row.names = colnames(scData))
  for (m in c("nCount_Spatial", "nFeature_Spatial", "percent_mito")) {
    col <- paste0(m, "_outliers")
    res[[paste0(m, "_SpotSweeper_outliers")]] <- if (col %in% colnames(colData(spe))) {
      as.logical(colData(spe)[[col]])
    } else {
      FALSE
    }
  }
  return(res)
}

##' @title Spatially aware spot QC
##' @return list with scData (filtered), scData.unfiltered,
##'   cellsPerGeneFraction and qcInfo.
spatial_qc <- function(scData, featInfo, param, spotInfo) {
  scData[["Spatial"]] <- AddMetaData(
    object = scData[["Spatial"]],
    metadata = featInfo[rownames(scData), ]
  )
  scData <- PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
  nMitoGenes <- sum(featInfo$isMito & rownames(featInfo) %in% rownames(scData))
  ## Probe panels (FFPE/CytAssist) have no mitochondrial probes: the mito
  ## fraction is then 0 everywhere and cannot flag damaged spots.
  useMito <- nMitoGenes >= 5 && stats::median(scData$percent_mito) > 0.1
  ssStats <- spatial_spot_sweeper(scData, useMito = useMito)
  scData <- AddMetaData(scData, ssStats)
  ## without mito probes the SpotSweeper mito flag is FALSE everywhere
  qcList <- filterCellsAndGenes(scData, param)
  if (ezIsSpecified(param$cellCycle) && param$cellCycle) {
    qcList$scData <- addCellCycleToSeurat(
      qcList$scData,
      param$refBuild,
      bpparam(),
      assay = "Spatial",
      method = "seurat"
    )
  }
  unf <- qcList$scData.unfiltered
  qcInfo <- list(
    useMito = useMito,
    nMitoGenes = nMitoGenes,
    probeBased = spotInfo$probeBased,
    nSpotsTissue = ncol(unf),
    nSpotsKept = ncol(qcList$scData),
    nGenesKept = nrow(qcList$scData),
    flags = colSums(as.matrix(unf[[c("qc.lib", "qc.nexprs", "qc.mito", "qc.ribo")]]))
  )
  qcList$qcInfo <- qcInfo
  return(qcList)
}

##' @title Expression-only clustering
spatial_expr_clusters <- function(scData, param) {
  scData <- seuratStandardSCTPreprocessing(scData, param, assay = "Spatial")
  scData <- seuratStandardWorkflow(scData, param, ident.name = "seurat_clusters")
  scData$expr_cluster <- Idents(scData)
  return(scData)
}

##' @title Spot neighbour graphs stored in the Seurat object
spatial_add_graphs <- function(scData) {
  for (r in 1:2) {
    adj <- spatial_hex_graph(
      scData$array_row,
      scData$array_col,
      colnames(scData),
      rings = r
    )
    scData[[paste0("spatial_hex", r)]] <- SeuratObject::as.Graph(adj)
  }
  return(scData)
}

spatial_get_adjacency <- function(scData, rings = 1) {
  adj <- methods::as(scData[[paste0("spatial_hex", rings)]], "dgCMatrix")
  return(adj)
}

##' @title BANKSY spatial domains
##' @description Augments each spot's expression with the mean expression of
##'   its neighbours (and, optionally, the neighbourhood gradient), weighted
##'   by lambda. lambda = 0 is plain expression clustering, lambda = 0.8 is the
##'   BANKSY domain-segmentation setting.
##' @return list with scData (banksy_domain, banksy_lambda_* columns,
##'   pca.banksy and umap.banksy) and a summary of the lambda sweep.
spatial_banksy_domains <- function(scData, param) {
  lambda <- as.numeric(param$banksyLambda)
  kGeom <- as.integer(param$banksyKgeom)
  resolution <- as.numeric(param$banksyResolution)
  npcs <- as.integer(param$banksyNpcs)
  sweep <- sort(unique(c(0, 0.2, 0.5, 0.8, lambda)))
  adj <- spatial_get_adjacency(scData, rings = 1)
  run_one <- function(obj, lam) {
    obj <- RunBanksy_v5(
      obj,
      lambda = lam,
      assay = "SCT",
      layer = "data",
      features = "variable",
      dimx = "pxl_col_fullres",
      dimy = "pxl_row_fullres",
      k_geom = kGeom,
      assay_name = "BANKSY",
      verbose = FALSE
    )
    obj <- RunPCA(
      obj,
      assay = "BANKSY",
      reduction.name = "pca.banksy",
      reduction.key = "BANKSYPC_",
      features = rownames(obj[["BANKSY"]]),
      npcs = npcs,
      verbose = FALSE
    )
    obj <- FindNeighbors(
      obj,
      reduction = "pca.banksy",
      dims = seq_len(npcs),
      graph.name = c("banksy_nn", "banksy_snn"),
      verbose = FALSE
    )
    obj <- FindClusters(
      obj,
      graph.name = "banksy_snn",
      cluster.name = "banksy_tmp",
      resolution = resolution,
      random.seed = 38,
      verbose = FALSE
    )
    obj
  }
  sweepRes <- list()
  for (lam in sweep) {
    tmp <- run_one(scData, lam)
    labels <- paste0("D", as.integer(as.character(tmp$banksy_tmp)) + 1)
    colName <- paste0("banksy_lambda_", lam)
    scData[[colName]] <- factor(
      labels,
      levels = paste0("D", seq_along(unique(labels)))
    )
    sweepRes[[colName]] <- data.frame(
      lambda = lam,
      nDomains = length(unique(labels)),
      smoothness = mean(spatial_label_smoothness(labels, adj), na.rm = TRUE),
      ariVsExprClusters = adjusted_rand_index(labels, scData$expr_cluster)
    )
    if (lam == lambda) {
      scData[["pca.banksy"]] <- tmp[["pca.banksy"]]
      scData[["banksy_snn"]] <- tmp[["banksy_snn"]]
    }
    remove(tmp)
  }
  scData$banksy_domain <- scData[[paste0("banksy_lambda_", lambda)]][, 1]
  scData <- RunUMAP(
    scData,
    reduction = "pca.banksy",
    dims = seq_len(npcs),
    reduction.name = "umap.banksy",
    reduction.key = "BANKSYUMAP_",
    verbose = FALSE
  )
  ## FindClusters() in Seurat 5.5 overwrites seurat_clusters; restore it
  scData$seurat_clusters <- scData$expr_cluster
  DefaultAssay(scData) <- "SCT"
  Idents(scData) <- "expr_cluster"
  return(list(
    scData = scData,
    sweep = do.call(rbind, sweepRes),
    lambda = lambda
  ))
}

##' @title Positive markers for a grouping
spatial_group_markers <- function(scData, groupBy, param) {
  Idents(scData) <- groupBy
  markers <- FindAllMarkers(
    object = scData,
    assay = "SCT",
    test.use = param$DE.method,
    only.pos = TRUE,
    min.pct = param$min.pct,
    return.thresh = param$pvalue_allMarkers,
    logfc.threshold = param$logfc.threshold,
    verbose = FALSE
  )
  if (nrow(markers) == 0) {
    return(data.frame())
  }
  cm <- markers[, c("gene", "cluster", "pct.1", "pct.2", "avg_log2FC", "p_val_adj")]
  cm$diff_pct <- abs(cm$pct.1 - cm$pct.2)
  cm <- cm[cm$p_val_adj < param$pvalue_allMarkers, ]
  cm <- cm[order(cm$cluster, -cm$avg_log2FC), ]
  rownames(cm) <- NULL
  return(cm)
}

##' @title Read a region marker table
##' @param markerFile TSV with columns region and gene; "" uses the built-in
##'   mouse brain panel.
spatial_region_markers <- function(markerFile, species) {
  if (!ezIsSpecified(markerFile) || markerFile == "") {
    if (species != "Mouse") {
      return(NULL)
    }
    markerFile <- system.file(
      "extdata/regionMarkers_mouseBrain.tsv",
      package = "ezRun",
      mustWork = TRUE
    )
  }
  tab <- utils::read.delim(markerFile, stringsAsFactors = FALSE)
  stopifnot(all(c("region", "gene") %in% colnames(tab)))
  return(split(tab$gene, tab$region))
}

##' @title Domain labels of the spots
##' @return the annotated domain labels if available, else the BANKSY domains.
spatial_domain_labels <- function(scData) {
  if ("domain_label" %in% colnames(scData[[]])) {
    return(scData$domain_label)
  }
  return(scData$banksy_domain)
}

##' @title Annotate BANKSY domains with anatomical regions
##' @description Scores every spot for each region marker set (UCell), then
##'   averages the z-scored region scores per domain and labels each domain
##'   with its best region (see \code{assign_domain_labels}).
spatial_annotate_domains <- function(scData, param) {
  species <- getSpecies(param$refBuild)
  regionSets <- spatial_region_markers(param$regionMarkerFile, species)
  if (is.null(regionSets)) {
    stop("no region marker file given and no built-in panel for ", species)
  }
  regionSets <- lapply(regionSets, intersect, rownames(scData))
  regionSets <- regionSets[lengths(regionSets) >= 2]
  ## UCell column names must be syntactic; keep the readable region names
  regionLabels <- stats::setNames(names(regionSets), make.names(names(regionSets)))
  names(regionSets) <- names(regionLabels)
  ucell <- UCell::ScoreSignatures_UCell(
    GetAssayData(scData, assay = "SCT", layer = "data"),
    features = regionSets,
    name = "_region"
  )
  ucell <- as.data.frame(ucell)
  rownames(ucell) <- colnames(scData)
  scData <- AddMetaData(scData, ucell)
  scoreCols <- paste0(names(regionSets), "_region")
  scores <- as.matrix(scData[[scoreCols]])
  colnames(scores) <- regionLabels[names(regionSets)]
  z <- scale(scores)
  domainMeans <- apply(z, 2, function(x) tapply(x, scData$banksy_domain, mean))
  if (is.null(dim(domainMeans))) {
    domainMeans <- matrix(domainMeans, nrow = 1, dimnames = list(levels(scData$banksy_domain), colnames(z)))
  }
  labels <- assign_domain_labels(domainMeans, minMargin = 0.25)
  ## unique, readable labels: "D3 Cortex L5"
  labels$domain_label <- paste(labels$domain, labels$label)
  domainLabel <- factor(
    labels$domain_label[match(scData$banksy_domain, labels$domain)],
    levels = labels$domain_label[order(as.integer(sub("^D", "", labels$domain)))]
  )
  scData <- AddMetaData(
    scData,
    data.frame(
      domain_label = domainLabel,
      spot_region = colnames(z)[max.col(z, ties.method = "first")],
      row.names = colnames(scData)
    )
  )
  return(list(
    scData = scData,
    domainScores = domainMeans,
    domainLabels = labels,
    regionSets = regionSets
  ))
}

##' @title Spatially variable genes and spatial expression patterns
spatial_svg <- function(scData, param) {
  adj <- spatial_get_adjacency(scData, rings = 1)
  dat <- GetAssayData(scData, assay = "SCT", layer = "data")
  detected <- Matrix::rowSums(dat > 0) >= max(10, 0.01 * ncol(dat))
  dat <- dat[detected, , drop = FALSE]
  moran <- spatial_morans_i(dat, adj)
  moran$isVariableFeature <- moran$gene %in% VariableFeatures(scData)
  topN <- as.integer(param$nnsvgTopN)
  nnsvgRes <- NULL
  if (ezIsSpecified(param$computeNnSVG) && param$computeNnSVG && topN > 0) {
    genes <- head(moran$gene[!is.na(moran$moransI)], topN)
    y <- as.matrix(dat[genes, , drop = FALSE])
    coords <- cbind(scData$pxl_col_fullres, scData$pxl_row_fullres)
    nn <- nnSVG::nnSVG(
      y,
      spatial_coords = coords,
      n_threads = as.integer(param$cores)
    )
    nnsvgRes <- data.frame(
      gene = rownames(nn),
      nnSVG_rank = nn[, "rank"],
      nnSVG_LR_stat = nn[, "LR_stat"],
      nnSVG_prop_sv = nn[, "prop_sv"],
      nnSVG_padj = nn[, "padj"]
    )
    moran <- merge(moran, nnsvgRes, by = "gene", all.x = TRUE, sort = FALSE)
    moran <- moran[order(moran$moransI, decreasing = TRUE), ]
  }
  ## Group the top SVGs into spatial patterns: cluster their smoothed,
  ## scaled profiles and summarise every pattern by its mean profile.
  nPatternGenes <- min(500, sum(moran$fdr < 0.05, na.rm = TRUE))
  patterns <- NULL
  if (nPatternGenes >= 20) {
    genes <- moran$gene[seq_len(nPatternGenes)]
    sm <- spatial_smooth(dat[genes, , drop = FALSE], adj)
    sm <- t(scale(t(sm)))
    sm[!is.finite(sm)] <- 0
    k <- min(as.integer(param$svgPatterns), nPatternGenes %/% 5)
    hc <- stats::hclust(stats::as.dist(1 - stats::cor(t(sm))), method = "ward.D2")
    cl <- stats::cutree(hc, k = k)
    ## order patterns by size
    ord <- order(table(cl), decreasing = TRUE)
    cl <- match(cl, as.integer(names(table(cl)))[ord])
    names(cl) <- genes
    patternScores <- sapply(seq_len(k), function(p) {
      colMeans(sm[names(cl)[cl == p], , drop = FALSE])
    })
    colnames(patternScores) <- paste0("svg_pattern_", seq_len(k))
    scData <- AddMetaData(scData, as.data.frame(patternScores))
    moran$svgPattern <- NA_integer_
    moran$svgPattern[match(names(cl), moran$gene)] <- cl
    patterns <- list(k = k, genes = split(names(cl), cl))
  }
  rownames(moran) <- NULL
  return(list(scData = scData, table = moran, patterns = patterns))
}

##' @title Cell-type deconvolution with RCTD
spatial_rctd <- function(scData, featInfo, param) {
  refPath <- resolve_rctd_reference_path(param)
  reference <- load_rctd_reference(refPath)
  counts <- GetAssayData(scData, assay = "Spatial", layer = "counts")
  geneIds <- stats::setNames(featInfo$gene_id, rownames(featInfo))
  counts <- match_rctd_genes(counts, geneIds, rownames(reference@counts))
  coords <- data.frame(
    x = scData$pxl_col_fullres,
    y = scData$pxl_row_fullres,
    row.names = colnames(scData)
  )
  res <- run_rctd_spots(
    counts,
    coords,
    reference,
    cores = as.integer(param$cores),
    umiMin = as.numeric(param$rctdUMImin),
    doubletMode = "full"
  )
  weights <- res$weights
  ## spots that RCTD dropped (too few UMIs) get NA
  full <- matrix(
    NA_real_,
    nrow = ncol(scData),
    ncol = ncol(weights),
    dimnames = list(colnames(scData), colnames(weights))
  )
  full[rownames(weights), ] <- weights
  cellTypes <- make.unique(gsub("[ _]+", "-", colnames(full)))
  colnames(full) <- cellTypes
  wMat <- t(full)
  wMat[is.na(wMat)] <- 0
  scData[["RCTD"]] <- SeuratObject::CreateAssay5Object(data = methods::as(wMat, "dgCMatrix"))
  dominant <- colnames(full)[max.col(replace(full, is.na(full), -1), ties.method = "first")]
  dominant[rowSums(is.na(full)) > 0] <- NA
  scData$rctd_dominant <- factor(dominant)
  scData$rctd_dominant_weight <- apply(full, 1, max)
  composition <- stats::aggregate(
    full,
    by = list(domain = as.character(spatial_domain_labels(scData))),
    FUN = mean,
    na.rm = TRUE
  )
  return(list(
    scData = scData,
    weights = full,
    composition = composition,
    reference = basename(refPath)
  ))
}

##' @title Domain neighbourhoods and cell-type co-localisation
spatial_neighbourhood <- function(scData, weights = NULL) {
  adj <- spatial_get_adjacency(scData, rings = 1)
  domains <- spatial_domain_labels(scData)
  enrichment <- spatial_adjacency_enrichment(domains, adj, nPerm = 500)
  coloc <- NULL
  if (!is.null(weights)) {
    ok <- stats::complete.cases(weights)
    w <- weights[ok, , drop = FALSE]
    ## ignore cell types that are essentially absent from this section
    w <- w[, colMeans(w) > 0.01, drop = FALSE]
    coloc <- spatial_colocalisation(w, adj[ok, ok])
  }
  return(list(domainAdjacency = enrichment, colocalisation = coloc))
}

##' @title Pathway (PROGENy) and TF (DoRothEA) activities per spot
spatial_activity <- function(scData, param) {
  species <- getSpecies(param$refBuild)
  if (!species %in% c("Human", "Mouse")) {
    stop("pathway/TF activity is only available for human and mouse")
  }
  mat <- as.matrix(GetAssayData(scData, assay = "SCT", layer = "data"))
  if (species == "Mouse") {
    progenyNet <- progeny::getModel(organism = "Mouse", top = 500)
    utils::data("dorothea_mm", package = "dorothea", envir = environment())
    dorotheaNet <- get("dorothea_mm")
  } else {
    progenyNet <- progeny::getModel(organism = "Human", top = 500)
    utils::data("dorothea_hs", package = "dorothea", envir = environment())
    dorotheaNet <- get("dorothea_hs")
  }
  progenyNet <- data.frame(
    target = rep(rownames(progenyNet), ncol(progenyNet)),
    source = rep(colnames(progenyNet), each = nrow(progenyNet)),
    mor = as.vector(as.matrix(progenyNet))
  )
  progenyNet <- progenyNet[progenyNet$mor != 0, ]
  dorotheaNet <- dorotheaNet[dorotheaNet$confidence %in% c("A", "B", "C"), ]
  dorotheaNet <- data.frame(
    source = dorotheaNet$tf,
    target = dorotheaNet$target,
    mor = dorotheaNet$mor
  )
  to_assay <- function(acts) {
    m <- tidyr::pivot_wider(
      acts[, c("source", "condition", "score")],
      names_from = "condition",
      values_from = "score"
    )
    rn <- m$source
    m <- as.matrix(m[, colnames(scData)])
    rownames(m) <- gsub("_", "-", rn)
    m[is.na(m)] <- 0
    SeuratObject::CreateAssay5Object(data = methods::as(m, "dgCMatrix"))
  }
  progenyAct <- decoupleR::run_ulm(mat, progenyNet, minsize = 5)
  dorotheaAct <- decoupleR::run_ulm(mat, dorotheaNet, minsize = 10)
  scData[["progeny"]] <- to_assay(progenyAct)
  scData[["dorothea"]] <- to_assay(dorotheaAct)
  domains <- spatial_domain_labels(scData)
  domain_means <- function(assayName) {
    m <- as.matrix(GetAssayData(scData, assay = assayName, layer = "data"))
    t(apply(m, 1, function(x) tapply(x, domains, mean)))
  }
  return(list(
    scData = scData,
    progenyPerDomain = domain_means("progeny"),
    dorotheaPerDomain = domain_means("dorothea")
  ))
}

##' @title Distance-constrained ligand-receptor analysis between domains
##' @description CellChat in spatial mode: communication probabilities are
##'   only computed between groups within the interaction range, which makes
##'   the result more plausible than non-spatial CellChat. Spots are
##'   multi-cell mixtures, so this is illustrative at domain level.
spatial_communication <- function(scData, param, spotInfo) {
  library(CellChat)
  domains <- droplevels(spatial_domain_labels(scData))
  ## CellChat needs group names without special characters
  groupNames <- make.names(levels(domains))
  meta <- data.frame(
    labels = factor(groupNames[as.integer(domains)], levels = groupNames),
    row.names = colnames(scData)
  )
  dataInput <- GetAssayData(scData, assay = "SCT", layer = "data")
  coords <- data.frame(
    imagerow = scData$pxl_row_fullres,
    imagecol = scData$pxl_col_fullres,
    row.names = colnames(scData)
  )
  sf <- spotInfo$scaleFactors
  scaleFactors <- list(
    spot.diameter = 55,
    spot = sf$spot_diameter_fullres,
    fiducial = sf$fiducial_diameter_fullres,
    hires = sf$tissue_hires_scalef,
    lowres = sf$tissue_lowres_scalef
  )
  cellchat <- createCellChat(
    object = dataInput,
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = coords,
    scale.factors = scaleFactors
  )
  species <- getSpecies(param$refBuild)
  db <- if (species == "Mouse") CellChatDB.mouse else CellChatDB.human
  cellchat@DB <- subsetDB(db, search = "Secreted Signaling")
  cellchat <- subsetData(cellchat)
  future::plan("sequential")
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = 0.1,
    distance.use = TRUE,
    interaction.length = as.numeric(param$cellchatInteractionLength),
    scale.distance = 0.01
  )
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  comm <- subsetCommunication(cellchat)
  pathways <- subsetCommunication(cellchat, slot.name = "netP")
  ## back to readable group labels
  lookup <- stats::setNames(levels(domains), groupNames)
  for (col in c("source", "target")) {
    comm[[col]] <- lookup[as.character(comm[[col]])]
    pathways[[col]] <- lookup[as.character(pathways[[col]])]
  }
  return(list(
    cellchat = cellchat,
    interactions = comm,
    pathways = pathways,
    groupLookup = lookup
  ))
}

##' @title Azimuth label transfer (kept for non-brain tissues)
spatial_azimuth <- function(scData, param) {
  library(Azimuth)
  reference <- gsub(" \\(human\\)| \\(mouse\\)", "", param$Azimuth)
  environment(MyDietSeurat) <- asNamespace("Seurat")
  assignInNamespace("DietSeurat", MyDietSeurat, ns = "Seurat")
  tmp <- scData
  tmp[["RNA"]] <- JoinLayers(tmp[["Spatial"]])
  tmp <- Azimuth::RunAzimuth(tmp, reference, assay = "RNA")
  aziNames <- setdiff(colnames(tmp@meta.data), colnames(scData@meta.data))
  aziNames <- grep("^predicted\\.", aziNames, value = TRUE)
  aziNames <- grep("score$", aziNames, value = TRUE, invert = TRUE)
  aziResults <- tmp@meta.data[, aziNames, drop = FALSE]
  colnames(aziResults) <- sub("^predicted\\.", "Azimuth.", colnames(aziResults))
  return(aziResults)
}

##' @title Write a Loupe Browser file
##' @description loupeR converts every factor column and assay of the object,
##'   which takes very long on the full object; a slim copy with the counts,
##'   the clusterings and the UMAPs is enough for Loupe.
spatial_write_cloupe <- function(scData, sampleName) {
  keepCols <- intersect(
    c("expr_cluster", "banksy_domain", "domain_label", "rctd_dominant"),
    colnames(scData[[]])
  )
  slim <- CreateSeuratObject(
    counts = GetAssayData(scData, assay = "Spatial", layer = "counts"),
    assay = "Spatial",
    meta.data = scData[[]][, keepCols, drop = FALSE]
  )
  for (r in intersect(c("umap", "umap.banksy"), SeuratObject::Reductions(scData))) {
    slim[[r]] <- scData[[r]]
  }
  loupeR::create_loupe_from_seurat(slim, output_name = sampleName)
}
