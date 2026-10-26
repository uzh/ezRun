###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## Cell-cell communication and niche models shared by the spatial apps
## (SpatialSeurat for Visium SD, VisiumHDSeurat). The functions take light
## inputs - a log-normalised matrix, coordinates in microns, labels and a
## spot x cell-type weight matrix - and never touch the Seurat object, so a
## Visium HD job does not copy its 100k-bin object.
##
## Cost rules, measured on the p37810 mouse brain sections:
## - CellChat v2 works at group level and finds the distance between groups
##   with nearest-neighbour search, so it scales to Visium HD; its cost is the
##   nboot permutations, which run in parallel.
## - SpatialCellChat and CARD keep an N x N spot distance matrix: fine for
##   Visium SD (~2,500 spots), on Visium HD only after aggregation into
##   pseudo-spots (aggregate_pseudospots).
## - mistyR fits one random forest per target; the paraview uses all spot
##   pairs unless nn is set.
## CellChat and SpatialCellChat export more than 100 functions with the same
## names, so both are always called with `::` and never attached.

##' @title Spot coordinates in microns
##' @param x,y full-resolution pixel coordinates.
##' @param micronsPerPixel scale (Visium SD: 55 / spot_diameter_fullres;
##'   Visium HD: microns_per_pixel of scalefactors_json.json).
##' @param spots spot names.
##' @return data.frame with x and y in microns, row names = spots.
spatial_coords_um <- function(x, y, micronsPerPixel, spots) {
  data.frame(
    x = as.numeric(x) * micronsPerPixel,
    y = as.numeric(y) * micronsPerPixel,
    row.names = spots
  )
}

##' @title Aggregate spots or bins into square pseudo-spots
##' @description One sparse indicator product: counts are summed, weights and
##'   coordinates averaged over the members. Rows of \code{weights} that are
##'   NA (spots RCTD did not fit) do not count towards the average.
##' @param counts gene x spot matrix (counts, or NULL).
##' @param coordsUm data.frame with x and y in microns, row names = spots.
##' @param sizeUm side length of a pseudo-spot in microns.
##' @param weights spot x cell-type matrix, or NULL.
##' @param minMembers pseudo-spots with fewer members are dropped.
##' @return list with counts, coords (microns), weights, nMembers and the
##'   pseudo-spot of every input spot (\code{membership}, NA when dropped).
aggregate_pseudospots <- function(
  counts = NULL,
  coordsUm,
  sizeUm,
  weights = NULL,
  minMembers = 1
) {
  gx <- floor(coordsUm$x / sizeUm)
  gy <- floor(coordsUm$y / sizeUm)
  unit <- factor(paste0("ps_", gx, "_", gy))
  ind <- Matrix::sparseMatrix(
    i = seq_along(unit),
    j = as.integer(unit),
    x = 1,
    dims = c(length(unit), nlevels(unit)),
    dimnames = list(rownames(coordsUm), levels(unit))
  )
  nMembers <- Matrix::colSums(ind)
  keep <- nMembers >= minMembers
  ind <- ind[, keep, drop = FALSE]
  nMembers <- nMembers[keep]
  coords <- as.matrix(Matrix::t(ind) %*% as.matrix(coordsUm[, c("x", "y")])) / nMembers
  res <- list(
    coords = as.data.frame(coords),
    nMembers = nMembers,
    membership = stats::setNames(
      ifelse(keep[as.integer(unit)], as.character(unit), NA_character_),
      rownames(coordsUm)
    )
  )
  if (!is.null(counts)) {
    res$counts <- counts[, rownames(coordsUm), drop = FALSE] %*% ind
  }
  if (!is.null(weights)) {
    w <- as.matrix(weights[rownames(coordsUm), , drop = FALSE])
    fitted <- stats::complete.cases(w)
    w[!fitted, ] <- 0
    nFitted <- as.vector(Matrix::t(ind) %*% as.numeric(fitted))
    wSum <- as.matrix(Matrix::t(ind) %*% w)
    wAvg <- wSum / nFitted
    wAvg[nFitted == 0, ] <- NA
    res$weights <- wAvg
  }
  return(res)
}

##' @title Stratified subsample of spots per label
##' @param labels factor or character, named by spot.
##' @param maxPerGroup maximal spots per label.
##' @return the names (or indices when unnamed) of the kept spots, in input order.
stratified_subsample <- function(labels, maxPerGroup, seed = 42) {
  idx <- seq_along(labels)
  withr::with_seed(seed, {
    keep <- unlist(lapply(split(idx, as.character(labels)), function(i) {
      if (length(i) > maxPerGroup) i[sort(sample.int(length(i), maxPerGroup))] else i
    }), use.names = FALSE)
  })
  keep <- sort(keep)
  if (!is.null(names(labels))) {
    return(names(labels)[keep])
  }
  return(keep)
}

##' @title Evaluate an expression with a temporary future plan
##' @description Forked workers when \code{workers > 1}. The previous plan is
##'   restored afterwards, so later steps keep the app's plan.
with_future_plan <- function(workers, expr) {
  workers <- as.integer(workers)
  oldPlan <- if (workers > 1) {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::sequential)
  }
  on.exit(future::plan(oldPlan), add = TRUE)
  force(expr)
}

##' @title Call a function at the top level of a fresh R process
##' @description SpatialCellChat registers a global progress handler before
##'   every parallel loop (progressr::handlers(global = TRUE)). R refuses
##'   that while tryCatch/withCallingHandlers are on the stack - always the
##'   case inside an app step, and inside callr. A separate Rscript calls the
##'   function at its top level; arguments and result go through qs2 files.
##' @param fun function to call (its namespace is loaded by the child).
##' @param args list of arguments.
##' @return the value of \code{fun}.
run_in_rscript <- function(fun, args, workDir = tempfile("rscript_")) {
  dir.create(workDir, recursive = TRUE)
  on.exit(unlink(workDir, recursive = TRUE), add = TRUE)
  inFile <- file.path(workDir, "input.qs2")
  outFile <- file.path(workDir, "result.qs2")
  scriptFile <- file.path(workDir, "run.R")
  qs2::qs_save(list(fun = fun, args = args), inFile)
  writeLines(
    c(
      paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
      paste0("input <- qs2::qs_read(", deparse(inFile), ")"),
      "result <- do.call(input$fun, input$args)",
      paste0("qs2::qs_save(result, ", deparse(outFile), ")")
    ),
    scriptFile
  )
  ezSystem(paste(file.path(R.home("bin"), "Rscript"), scriptFile))
  qs2::qs_read(outFile)
}

##' @title Matrix of the CellChat ligand-receptor database for a species
cellchat_db <- function(species, search, pkg = "CellChat") {
  dbName <- if (species == "Mouse") "CellChatDB.mouse" else "CellChatDB.human"
  e <- new.env()
  utils::data(list = dbName, package = pkg, envir = e)
  db <- get(dbName, envir = e)
  subsetFun <- getExportedValue(pkg, "subsetDB")
  if (pkg == "SpatialCellChat") {
    return(subsetFun(db, search = search, non_protein = FALSE))
  }
  return(subsetFun(db, search = search))
}

##' @title Distance scale for CellChat v2
##' @description CellChat refuses scaled group distances below 1 and suggests
##'   [1, 2]; the right scale depends on the platform and on how the groups
##'   are mixed. Computes the group distances once (nearest-neighbour search,
##'   cheap) and scales the smallest to 1.1.
cellchat_scale_distance <- function(coords, labels, interactionRange, tol, contactRange) {
  meta <- data.frame(group = labels, samples = factor("sample1"))
  d <- CellChat::computeRegionDistance(
    coordinates = coords,
    meta = meta,
    interaction.range = interactionRange,
    ratio = 1,
    tol = tol,
    contact.dependent = TRUE,
    contact.range = contactRange
  )$d.spatial
  diag(d) <- NaN
  dMin <- suppressWarnings(min(d, na.rm = TRUE))
  if (!is.finite(dMin) || dMin <= 0) {
    return(0.01)
  }
  return(1.1 / dMin)
}

##' @title Distance-constrained CellChat v2 between groups of spots
##' @param data log-normalised gene x spot matrix (symbols).
##' @param labels factor of spot groups, named by spot.
##' @param coordsUm data.frame with x and y in microns.
##' @param spotSizeUm spot or bin diameter (tolerance = half of it).
##' @param contactRangeUm range of contact-dependent signalling.
##' @param interactionRangeUm maximal distance of secreted signalling.
##' @param species "Mouse" or "Human".
##' @param dbSearch CellChatDB categories.
##' @param nboot permutations for the p-values (run in parallel).
##' @param minCells groups with fewer spots are left out.
##' @return list with the slim CellChat object, the L-R and pathway tables
##'   (readable group names) and the group lookup.
run_cellchat_v2 <- function(
  data,
  labels,
  coordsUm,
  spotSizeUm,
  contactRangeUm,
  interactionRangeUm = 250,
  species = "Mouse",
  dbSearch = "Secreted Signaling",
  nboot = 20,
  workers = 1,
  minCells = 10
) {
  labels <- as.factor(labels)
  tab <- table(labels)
  small <- names(tab)[tab < minCells]
  if (length(small) > 0) {
    ezLog("CellChat: leaving out groups with < ", minCells, " spots: ", paste(small, collapse = ", "))
  }
  spots <- names(labels)[!labels %in% small]
  labels <- droplevels(labels[spots])
  if (nlevels(labels) < 2) {
    stop("CellChat needs at least two groups with ", minCells, " or more spots")
  }
  ## CellChat needs group names without special characters
  groupNames <- make.names(levels(labels), unique = TRUE)
  lookup <- stats::setNames(levels(labels), groupNames)
  meta <- data.frame(
    labels = factor(groupNames[as.integer(labels)], levels = groupNames),
    samples = factor("sample1"),
    row.names = spots
  )
  coords <- as.matrix(coordsUm[spots, c("x", "y")])
  tol <- spotSizeUm / 2
  cellchat <- CellChat::createCellChat(
    object = data[, spots, drop = FALSE],
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = coords,
    spatial.factors = data.frame(ratio = 1, tol = tol)
  )
  cellchat@DB <- cellchat_db(species, dbSearch)
  cellchat <- CellChat::subsetData(cellchat)
  scaleDistance <- cellchat_scale_distance(
    coords, meta$labels, interactionRangeUm, tol, contactRangeUm
  )
  ezLog(sprintf("CellChat: %d spots, %d groups, scale.distance %.4g", length(spots), nlevels(labels), scaleDistance))
  cellchat <- with_future_plan(workers, {
    cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
    cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
    CellChat::computeCommunProb(
      cellchat,
      type = "truncatedMean",
      trim = 0.1,
      distance.use = TRUE,
      interaction.range = interactionRangeUm,
      scale.distance = scaleDistance,
      contact.dependent = TRUE,
      contact.range = contactRangeUm,
      nboot = nboot
    )
  })
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = minCells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  comm <- CellChat::subsetCommunication(cellchat)
  pathways <- CellChat::subsetCommunication(cellchat, slot.name = "netP")
  ## back to readable group labels
  for (col in c("source", "target")) {
    comm[[col]] <- unname(lookup[as.character(comm[[col]])])
    pathways[[col]] <- unname(lookup[as.character(pathways[[col]])])
  }
  return(list(
    cellchat = cellchat,
    interactions = comm,
    pathways = pathways,
    groupLookup = lookup,
    params = list(
      nSpots = length(spots), nGroups = nlevels(labels), nboot = nboot,
      interactionRange = interactionRangeUm, contactRange = contactRangeUm,
      scaleDistance = scaleDistance, db = paste(dbSearch, collapse = ", "),
      version = as.character(utils::packageVersion("CellChat"))
    )
  ))
}

##' @title Cell types for a weights-based analysis
##' @description Keeps the cell types with a mean weight of at least
##'   \code{minMean} that dominate at least \code{minSpots} spots.
select_weight_types <- function(weights, minMean = 0.01, minSpots = 10) {
  w <- weights[stats::complete.cases(weights), , drop = FALSE]
  dominant <- colnames(w)[max.col(w, ties.method = "first")]
  nDominant <- table(factor(dominant, levels = colnames(w)))
  colnames(w)[colMeans(w) >= minMean & as.vector(nDominant) >= minSpots]
}

##' @title SpatialCellChat for spot data with cell-type weights
##' @description Infers communication between neighbouring spots, then
##'   averages it to cell types weighted by the deconvolution
##'   (computeAvgCommunProb_Visium). Every L-R pair costs an N x N outer
##'   product of ligand and receptor expression, so the input is capped at
##'   \code{maxSpots} (Visium HD goes through aggregate_pseudospots first)
##'   and the defaults keep the data sparse and the pairs few: no ALRA
##'   imputation (dense data made Rep1 run > 25 min for 2,059 pairs),
##'   secreted signalling, ligand and receptor both variable.
##' @param data log-normalised gene x spot matrix (symbols).
##' @param weights spot x cell-type matrix (rows sum to 1, NA = not fitted).
##' @param minWeight smaller weights are set to 0, as in the SpatialCellChat
##'   tutorial.
##' @return list with the L-R and pathway tables, per-spot sender and
##'   receiver scores of the top pathways and the parameters.
run_spatialcellchat <- function(
  data,
  weights,
  coordsUm,
  spotSizeUm,
  species = "Mouse",
  dbSearch = "Secreted Signaling",
  interactionRangeUm = 250,
  contactRangeUm = 10,
  scaleDistance = 0.2,
  nboot = 20,
  workers = 1,
  minWeight = 0.1,
  minSpots = 10,
  maxSpots = 15000,
  nTopPathways = 4,
  impute = FALSE,
  variableBoth = TRUE
) {
  spots <- rownames(weights)[stats::complete.cases(weights)]
  spots <- intersect(spots, colnames(data))
  if (length(spots) > maxSpots) {
    stop(sprintf(
      "SpatialCellChat keeps an N x N spot matrix: %d spots exceed the limit of %d (aggregate first)",
      length(spots), maxSpots
    ))
  }
  types <- select_weight_types(weights[spots, , drop = FALSE], minSpots = minSpots)
  if (length(types) < 2) {
    stop("fewer than two cell types dominate at least ", minSpots, " spots")
  }
  w <- as.matrix(weights[spots, types, drop = FALSE])
  w[w < minWeight] <- 0
  spots <- spots[rowSums(w) > 0]
  w <- w[spots, , drop = FALSE]
  typeNames <- make.names(types, unique = TRUE)
  lookup <- stats::setNames(types, typeNames)
  colnames(w) <- typeNames
  labels <- factor(typeNames[max.col(w, ties.method = "first")], levels = typeNames)
  meta <- data.frame(labels = labels, samples = factor("sample1"), row.names = spots)
  coords <- as.matrix(coordsUm[spots, c("x", "y")])
  chat <- SpatialCellChat::createSpatialCellChat(
    object = data[, spots, drop = FALSE],
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = coords,
    spatial.factors = list(ratio = 1, tol = spotSizeUm / 2)
  )
  chat@DB <- cellchat_db(species, dbSearch, pkg = "SpatialCellChat")
  chat <- with_future_plan(workers, {
    chat <- SpatialCellChat::subsetData(chat)
    if (impute) {
      chat <- SpatialCellChat::preProcessing(chat)
    }
    chat <- SpatialCellChat::identifyOverExpressedGenes(chat, selection.method = "wilcox", do.grid = FALSE)
    chat <- SpatialCellChat::identifyOverExpressedInteractions(chat, variable.both = variableBoth)
    chat <- SpatialCellChat::computeCommunProb(
      chat,
      distance.use = TRUE,
      scale.distance = scaleDistance,
      contact.dependent = TRUE,
      interaction.range = interactionRangeUm,
      contact.range = contactRangeUm
    )
    chat <- SpatialCellChat::filterProbability(chat, nboot = nboot)
    chat <- SpatialCellChat::filterCommunication(chat, min.cells = NULL, min.links = 10, min.cells.sr = 10)
    SpatialCellChat::computeAvgCommunProb_Visium(
      chat,
      cell.type.decomposition = w,
      avg.type = "avg",
      nboot = nboot,
      do.permutation = TRUE
    )
  })
  chat <- SpatialCellChat::filterCommunication(chat, min.cells = 10, min.links = NULL, min.cells.sr = NULL)
  comm <- SpatialCellChat::subsetCommunication(chat)
  chat <- SpatialCellChat::computeCommunProbPathway(chat)
  chat <- SpatialCellChat::aggregateNet(chat)
  pathways <- SpatialCellChat::subsetCommunication(chat, slot.name = "netP")
  for (col in c("source", "target")) {
    comm[[col]] <- unname(lookup[as.character(comm[[col]])])
    pathways[[col]] <- unname(lookup[as.character(pathways[[col]])])
  }
  spotScores <- spatialcellchat_spot_scores(chat, pathways, nTopPathways)
  return(list(
    interactions = comm,
    pathways = pathways,
    spotScores = spotScores,
    params = list(
      nSpots = length(spots), nTypes = length(types), nboot = nboot,
      interactionRange = interactionRangeUm, contactRange = contactRangeUm,
      scaleDistance = scaleDistance, minWeight = minWeight,
      db = paste(dbSearch, collapse = ", "),
      version = as.character(utils::packageVersion("SpatialCellChat"))
    )
  ))
}

##' @title Per-spot sender and receiver scores of the top pathways
##' @description Sums the spot-to-spot communication probabilities of all
##'   L-R pairs of a pathway: outgoing = row sums, incoming = column sums.
##' @return long data.frame (spot, pathway, outgoing, incoming) or NULL.
spatialcellchat_spot_scores <- function(chat, pathways, nTop = 4) {
  probList <- chat@net$tmp$prob.cell
  if (is.null(probList) || nrow(pathways) == 0) {
    return(NULL)
  }
  strength <- tapply(pathways$prob, pathways$pathway_name, sum)
  top <- names(sort(strength, decreasing = TRUE))[seq_len(min(nTop, length(strength)))]
  lrPathway <- chat@DB$interaction[names(probList), "pathway_name"]
  spots <- rownames(chat@images$coordinates)
  do.call(rbind, lapply(top, function(p) {
    mats <- probList[which(lrPathway == p)]
    if (length(mats) == 0) {
      return(NULL)
    }
    total <- Reduce(`+`, mats)
    data.frame(
      spot = spots,
      pathway = p,
      outgoing = Matrix::rowSums(total),
      incoming = Matrix::colSums(total)
    )
  }))
}

##' @title MISTy spatial dependencies between cell types
##' @description Explains every cell type's weight from the other types in the
##'   same spot (intraview), the direct neighbours (juxtaview) and the wider
##'   neighbourhood (paraview, Gaussian kernel). The gain in R2 over the
##'   intraview shows how much the surroundings add.
##' @param weights spot x cell-type matrix (NA rows = not fitted).
##' @param coordsUm data.frame with x and y in microns.
##' @param spacingUm centre-to-centre distance of neighbouring spots.
##' @param paraRangeUm length scale of the paraview kernel.
##' @param nn paraview neighbours (NULL = all spot pairs; set it for large
##'   inputs).
##' @return list with improvements, contributions and importances (readable
##'   cell-type names) and the parameters.
run_misty_weights <- function(
  weights,
  coordsUm,
  spacingUm,
  paraRangeUm = 300,
  nn = NULL,
  workers = 1,
  cvFolds = 5,
  numTrees = 100,
  minMean = 0.01,
  seed = 42
) {
  w <- weights[stats::complete.cases(weights), , drop = FALSE]
  types <- colnames(w)[colMeans(w) >= minMean &
    apply(w, 2, function(x) length(unique(x))) > 2 * cvFolds]
  if (length(types) < 3) {
    stop("MISTy needs at least three cell types with a mean weight of ", minMean)
  }
  w <- w[, types, drop = FALSE]
  typeNames <- make.names(types, unique = TRUE)
  lookup <- stats::setNames(types, typeNames)
  colnames(w) <- typeNames
  pos <- as.data.frame(coordsUm[rownames(w), c("x", "y")])
  resDir <- tempfile("misty_")
  on.exit(unlink(resDir, recursive = TRUE), add = TRUE)
  with_future_plan(workers, {
    views <- mistyR::create_initial_view(tibble::as_tibble(w))
    views <- mistyR::add_juxtaview(views, positions = pos, neighbor.thr = 1.5 * spacingUm, verbose = FALSE)
    views <- mistyR::add_paraview(views, positions = pos, l = paraRangeUm, zoi = spacingUm, nn = nn, verbose = FALSE)
    mistyR::run_misty(
      views,
      results.folder = resDir,
      seed = seed,
      cv.folds = cvFolds,
      num.trees = numTrees
    )
  })
  res <- mistyR::collect_results(resDir)
  rename <- function(x) {
    hit <- x %in% names(lookup)
    x[hit] <- lookup[x[hit]]
    x
  }
  improvements <- as.data.frame(res$improvements)
  improvements$target <- rename(improvements$target)
  contributions <- as.data.frame(res$contributions)
  contributions$target <- rename(contributions$target)
  importances <- as.data.frame(res$importances)
  importances$Predictor <- rename(importances$Predictor)
  importances$Target <- rename(importances$Target)
  importances$sample <- NULL
  return(list(
    improvements = improvements,
    contributions = contributions,
    importances = importances,
    params = list(
      nSpots = nrow(w), nTypes = length(types), spacing = spacingUm,
      paraRange = paraRangeUm, nn = nn, cvFolds = cvFolds, numTrees = numTrees,
      version = as.character(utils::packageVersion("mistyR"))
    )
  ))
}

##' @title Consensus ligand-receptor ranking with liana
##' @description Fast methods only (no CellPhoneDB permutations). Not
##'   spatial by itself; \code{cellAdj} (source, target, adjacency) restricts
##'   the pairs to groups that touch.
##' @param counts,data raw and log-normalised gene x spot matrices.
##' @param labels factor of spot groups, named by spot.
##' @return the aggregated liana table (aggregate_rank = consensus rank).
run_liana_fast <- function(
  counts,
  data,
  labels,
  species = "Mouse",
  cellAdj = NULL,
  methods = c("natmi", "connectome", "logfc", "sca"),
  minCells = 10
) {
  labels <- as.factor(labels)
  tab <- table(labels)
  spots <- names(labels)[labels %in% names(tab)[tab >= minCells]]
  labels <- droplevels(labels[spots])
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts[, spots, drop = FALSE], logcounts = data[, spots, drop = FALSE]),
    colData = S4Vectors::DataFrame(group = labels, row.names = spots)
  )
  resource <- if (species == "Mouse") "MouseConsensus" else "Consensus"
  ## OmnipathR (used by liana) logs into ./omnipathr-log, the report folder
  if (!dir.exists("omnipathr-log")) {
    on.exit(unlink("omnipathr-log", recursive = TRUE), add = TRUE)
  }
  res <- liana::liana_wrap(
    sce,
    method = methods,
    resource = resource,
    idents_col = "group",
    min_cells = minCells,
    cell.adj = cellAdj,
    verbose = FALSE
  )
  agg <- liana::liana_aggregate(res, verbose = FALSE)
  return(as.data.frame(agg))
}
