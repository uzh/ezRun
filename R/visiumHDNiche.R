###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## Niche and communication steps of the VisiumHDSeurat app. A Visium HD
## sample has ~100k (16 um) to ~450k (8 um) bins, so nothing here runs on
## every bin at the N x N level:
## - CellChat v2 works on group means; it gets at most cellchatMaxPerGroup
##   bins per RCTD cell type / BANKSY niche.
## - MISTy and SpatialCellChat run on square pseudo-spots
##   (aggregate_pseudospots), which keeps HD's resolution of the cell-type
##   weights but bounds the number of units.
## The functions take the inputs from visiumhd_niche_inputs() and leave the
## Seurat object alone.

##' @title Bin (or cell) size in microns from the binSize parameter
visiumhd_bin_um <- function(binSize) {
  um <- suppressWarnings(as.numeric(sub(".*square_0*([0-9]+)um.*", "\\1", binSize)))
  if (is.na(um)) {
    um <- 10 # segmented cells: a typical cell diameter
  }
  um
}

##' @title Microns per full-resolution pixel of a Visium HD sample
##' @description From scalefactors_json.json when Space Ranger wrote
##'   \code{microns_per_pixel}; otherwise from the bin grid (median
##'   nearest-neighbour distance = one bin).
visiumhd_microns_per_pixel <- function(spaceRangerDir, binSize, coordsPx) {
  candidates <- c(
    file.path(spaceRangerDir, binSize, "spatial", "scalefactors_json.json"),
    file.path(spaceRangerDir, "spatial", "scalefactors_json.json")
  )
  for (f in candidates[file.exists(candidates)]) {
    sf <- jsonlite::fromJSON(f)
    if (!is.null(sf$microns_per_pixel)) {
      return(as.numeric(sf$microns_per_pixel))
    }
  }
  idx <- sample.int(nrow(coordsPx), min(5000, nrow(coordsPx)))
  nn <- RANN::nn2(coordsPx[, c("x", "y")], coordsPx[idx, c("x", "y")], k = 2)$nn.dists[, 2]
  visiumhd_bin_um(binSize) / stats::median(nn)
}

##' @title Extract the inputs of the niche and communication steps once
##' @param weights bins x cell types (RCTD, rows sum to 1), or NULL.
##' @param resultsDf RCTD doublet-mode table (spot_class, first_type), or
##'   NULL.
##' @return list with data, counts, coordsUm, weights, firstType, niche,
##'   binUm and species.
visiumhd_niche_inputs <- function(scData, assay, weights, resultsDf, param, input) {
  coords <- SeuratObject::GetTissueCoordinates(scData)
  if ("cell" %in% colnames(coords)) {
    rownames(coords) <- coords$cell
  }
  colnames(coords)[1:2] <- c("x", "y")
  coords <- coords[colnames(scData), c("x", "y")]
  umPerPx <- visiumhd_microns_per_pixel(
    input$getFullPaths("SpaceRangerDir"),
    param$binSize,
    coords
  )
  firstType <- NULL
  if (!is.null(resultsDf)) {
    keep <- rownames(resultsDf)[resultsDf$spot_class != "reject"]
    firstType <- stats::setNames(as.character(resultsDf[keep, "first_type"]), keep)
  }
  list(
    data = SeuratObject::GetAssayData(scData, assay = assay, layer = "data"),
    counts = SeuratObject::GetAssayData(scData, assay = assay, layer = "counts"),
    coordsUm = spatial_coords_um(coords$x, coords$y, umPerPx, colnames(scData)),
    micronsPerPixel = umPerPx,
    weights = weights,
    firstType = firstType,
    niche = stats::setNames(as.character(scData$banksy_cluster), colnames(scData)),
    binUm = visiumhd_bin_um(param$binSize),
    species = getSpecies(param$refBuild)
  )
}

##' @title Run an optional niche/communication step
##' @param hd inputs from visiumhd_niche_inputs (NULL when they failed).
##' @param expr expression returning the step result (evaluated only when
##'   the step runs).
##' @param needsRctd the step needs the RCTD weights.
visiumhd_step <- function(name, flag, hd, expr, status, needsRctd = TRUE) {
  if (!(ezIsSpecified(flag) && isTRUE(as.logical(flag)))) {
    return(skip_step(name, "disabled by parameter", status))
  }
  if (is.null(hd)) {
    return(skip_step(name, "the spot coordinates could not be prepared", status))
  }
  if (needsRctd && is.null(hd$weights)) {
    return(skip_step(name, "needs the RCTD cell-type deconvolution", status))
  }
  run_optional_step(name, expr, status)
}

##' @title CellChat v2 between RCTD cell types and between BANKSY niches
visiumhd_communication <- function(hd, param) {
  maxPerGroup <- as.integer(param$cellchatMaxPerGroup)
  run <- function(labels) {
    labels <- labels[!is.na(labels)]
    keep <- stratified_subsample(labels, maxPerGroup)
    run_cellchat_v2(
      hd$data,
      factor(labels[keep]),
      hd$coordsUm,
      spotSizeUm = hd$binUm,
      contactRangeUm = hd$binUm,
      interactionRangeUm = as.numeric(param$cellchatInteractionLength),
      species = hd$species,
      nboot = as.integer(param$cellchatNboot)
    )
  }
  res <- list(niche = run(hd$niche))
  if (!is.null(hd$firstType)) {
    res$celltype <- run(hd$firstType)
  }
  res
}

##' @title RCTD weights on square pseudo-spots
visiumhd_pseudospots <- function(hd, sizeUm, withCounts = FALSE) {
  fitted <- rownames(hd$weights)
  aggregate_pseudospots(
    counts = if (withCounts) hd$counts[, fitted, drop = FALSE] else NULL,
    coordsUm = hd$coordsUm[fitted, , drop = FALSE],
    sizeUm = sizeUm,
    weights = hd$weights,
    minMembers = 2
  )
}

##' @title MISTy on pseudo-spots of RCTD weights
visiumhd_misty <- function(hd, param) {
  size <- as.numeric(param$nicheGridSize)
  ps <- visiumhd_pseudospots(hd, size)
  ezLog(sprintf("MISTy: %d bins -> %d pseudo-spots of %g um", nrow(hd$weights), nrow(ps$coords), size))
  ## nn neighbours per paraview; empty or 0 = all pseudo-spot pairs
  nn <- suppressWarnings(as.integer(param$mistyNeighbours))
  if (length(nn) != 1 || is.na(nn) || nn <= 0) {
    nn <- NULL
  }
  res <- run_misty_weights(
    ps$weights,
    ps$coords,
    spacingUm = size,
    paraRangeUm = as.numeric(param$mistyParaviewRange),
    nn = nn,
    workers = as.integer(param$cores)
  )
  res$params$gridSize <- size
  res$params$nBins <- nrow(hd$weights)
  res
}

##' @title SpatialCellChat on pseudo-spots
##' @description Counts are summed per pseudo-spot and log-normalised; the
##'   weights are averaged. The default 100 um grid gives Visium-SD-like
##'   units, for which SpatialCellChat was designed.
visiumhd_spatialcellchat <- function(hd, param) {
  size <- as.numeric(param$spatialCellChatGridSize)
  ps <- visiumhd_pseudospots(hd, size, withCounts = TRUE)
  ezLog(sprintf("SpatialCellChat: %d bins -> %d pseudo-spots of %g um", nrow(hd$weights), nrow(ps$coords), size))
  data <- Seurat::NormalizeData(ps$counts, verbose = FALSE)
  ## in its own Rscript: see run_in_rscript
  res <- run_in_rscript(run_spatialcellchat, list(
    data = data,
    weights = ps$weights,
    coordsUm = ps$coords,
    spotSizeUm = size,
    species = hd$species,
    interactionRangeUm = as.numeric(param$cellchatInteractionLength),
    nboot = as.integer(param$cellchatNboot),
    workers = as.integer(param$cores)
  ))
  res$params$gridSize <- size
  ## spot scores refer to pseudo-spots: keep their coordinates for the maps
  res$pseudospotCoords <- ps$coords
  res
}

##' @title liana between RCTD cell types on the CellChat subsample
visiumhd_liana <- function(hd, param) {
  labels <- hd$firstType
  keep <- stratified_subsample(labels, as.integer(param$cellchatMaxPerGroup))
  res <- run_liana_fast(
    hd$counts[, keep, drop = FALSE],
    hd$data[, keep, drop = FALSE],
    factor(labels[keep]),
    species = hd$species
  )
  res <- res[order(res$aggregate_rank), , drop = FALSE]
  utils::head(res, 1000)
}
