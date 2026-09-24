###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

addCellCycleToSCE <- function(sce, refBuild, BPPARAM) {
  counts <- counts(sce)
  rownames(counts) <- rowData(sce)$ID
  cellPhase <- getCellCycle(counts, refBuild, BPPARAM)
  if (!is.null(cellPhase)) {
    colData(sce)$CellCycle <- cellPhase$Phase
    colData(sce)$CellCycleG1 <- cellPhase$G1
    colData(sce)$CellCycleS <- cellPhase$S
    colData(sce)$CellCycleG2M <- cellPhase$G2M
  }
  return(sce)
}

addCellCycleToSeurat <- function(scData, refBuild, BPPARAM, assay = "RNA",
                                 method = c("cyclone", "seurat")) {
  if (match.arg(method) == "seurat") {
    return(addSeuratCellCycle(scData, refBuild, assay))
  }
  counts <- GetAssayData(scData, layer = "counts", assay = assay)
  metaFeatures <- scData[[assay]]@meta.data
  if ("gene_id" %in% names(metaFeatures)) {
    rownames(counts) <- metaFeatures$gene_id
  } else {
    rownames(counts) <- metaFeatures$ensemblID
  }
  cellPhase <- getCellCycle(counts, refBuild, BPPARAM)
  if (!is.null(cellPhase)) {
    cellcycleInfo = data.frame(
      CellCycle = cellPhase$Phase,
      CellCycleG1 = cellPhase$G1,
      CellCycleS = cellPhase$S,
      CellCycleG2M = cellPhase$G2M,
      CC.Difference = cellPhase$S - cellPhase$G2M,
      row.names = colnames(scData)
    )
    scData <- AddMetaData(scData, metadata = cellcycleInfo)
  }
  return(scData)
}

## Seurat's marker-module scoring (cc.genes.updated.2019). 19 s on 556k
## VisiumHD bins, where cyclone took 160 min at 4 cores (p41757 Batch2).
## Same columns as the cyclone path except CellCycleG1, which it has no score for.
addSeuratCellCycle <- function(scData, refBuild, assay) {
  genes <- Seurat::cc.genes.updated.2019
  species <- getSpecies(refBuild)
  if (species == "Mouse") {
    genes <- lapply(genes, function(g) paste0(substr(g, 1, 1), tolower(substring(g, 2))))
  } else if (species != "Human") {
    return(scData)
  }
  cc <- Seurat::CellCycleScoring(scData, s.features = genes$s.genes,
                                 g2m.features = genes$g2m.genes, assay = assay)
  AddMetaData(scData, data.frame(
    CellCycle = cc$Phase,
    CellCycleS = cc$S.Score,
    CellCycleG2M = cc$G2M.Score,
    CC.Difference = cc$S.Score - cc$G2M.Score,
    row.names = colnames(scData)
  ))
}

getCellCycle <- function(counts, refBuild, BPPARAM) {
  require(scran)

  species <- sub("\\/.*", "", refBuild)
  trainDataFile <- switch(
    species,
    Homo_sapiens = system.file(
      "exdata",
      "human_cycle_markers.rds",
      package = "scran"
    ),
    Mus_musculus = system.file(
      "exdata",
      "mouse_cycle_markers.rds",
      package = "scran"
    )
  )
  if (is.null(trainDataFile) || !file.exists(trainDataFile)) {
    return(NULL)
  } else {
    trainData <- readRDS(trainDataFile)
    cellCycleData <- cyclone(counts, trainData, BPPARAM = BPPARAM)
    cellPhase <- tibble(Name = colnames(counts), Phase = cellCycleData$phases)
    cellPhase <- bind_cols(cellPhase, cellCycleData$scores)
    return(cellPhase)
  }
}

getPerplexity <- function(n) {
  ifelse(n > 200, 30, 10)
}

SingleCorPlot <- function(
  data,
  col.by = NULL,
  cols = NULL,
  pt.size = NULL,
  smooth = FALSE,
  rows.highlight = NULL,
  legend.title = NULL,
  na.value = 'grey50',
  span = NULL
) {
  pt.size <- pt.size <- pt.size %||% AutoPointSize(data = data)
  orig.names <- colnames(x = data)
  names.plot <- colnames(x = data) <- gsub(
    pattern = '-',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  names.plot <- colnames(x = data) <- gsub(
    pattern = ':',
    replacement = '.',
    x = colnames(x = data),
    fixed = TRUE
  )
  if (ncol(x = data) < 2) {
    msg <- "Too few variables passed"
    if (ncol(x = data) == 1) {
      msg <- paste0(msg, ', only have ', colnames(x = data)[1])
    }
    stop(msg, call. = FALSE)
  }
  plot.cor <- round(x = cor(x = data[, 1], y = data[, 2]), digits = 2)
  if (!is.null(x = rows.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = rows.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size,
      cols.highlight = 'red',
      col.base = 'black',
      pt.size = pt.size
    )
    cols <- highlight.info$color
    col.by <- factor(
      x = highlight.info$highlight,
      levels = rev(x = highlight.info$plot.order)
    )
    plot.order <- order(col.by)
    data <- data[plot.order, ]
    col.by <- col.by[plot.order]
  }
  if (!is.null(x = col.by)) {
    data$colors <- col.by
  }
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = names.plot[1], y = names.plot[2])
  ) +
    labs(
      x = orig.names[1],
      y = orig.names[2],
      title = plot.cor,
      color = legend.title
    )
  if (smooth) {
    # density <- kde2d(x = data[, names.plot[1]], y = data[, names.plot[2]], h = Bandwidth(data = data[, names.plot]), n = 200)
    # density <- data.frame(
    #   expand.grid(
    #     x = density$x,
    #     y = density$y
    #   ),
    #   density = as.vector(x = density$z)
    # )
    plot <- plot +
      stat_density2d(
        mapping = aes(fill = ..density..^0.25),
        geom = 'tile',
        contour = FALSE,
        n = 200,
        h = Bandwidth(data = data[, names.plot])
      ) +
      # geom_tile(
      #   mapping = aes_string(
      #     x = 'x',
      #     y = 'y',
      #     fill = 'density'
      #   ),
      #   data = density
      # ) +
      scale_fill_continuous(low = 'white', high = 'dodgerblue4') +
      guides(fill = FALSE)
  }
  if (!is.null(x = col.by)) {
    plot <- plot +
      geom_point(
        mapping = aes_string(color = 'colors'),
        position = 'jitter',
        size = pt.size
      )
  } else {
    plot <- plot + geom_point(position = 'jitter', size = pt.size)
  }
  if (!is.null(x = cols)) {
    cols.scale <- if (
      length(x = cols) == 1 && cols %in% rownames(x = brewer.pal.info)
    ) {
      scale_color_brewer(palette = cols)
    } else {
      scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + cols.scale
    if (!is.null(x = rows.highlight)) {
      plot <- plot + guides(color = FALSE)
    }
  }
  plot <- plot + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
  if (!is.null(x = span)) {
    plot <- plot +
      geom_smooth(
        mapping = aes_string(x = names.plot[1], y = names.plot[2]),
        method = 'loess',
        span = span
      )
  }
  return(plot)
}


VariableFeaturePlot_mod <- function(
  sce,
  cols = c('black', 'red'),
  pt.size = 1,
  log = NULL
) {
  if (length(x = cols) != 2) {
    stop("'cols' must be of length 2")
  }
  hvf.info <- as.data.frame(rowData(sce))
  colnames(hvf.info) <- gsub(("sct."), "", colnames(hvf.info))
  var.status <- c('no', 'yes')[unlist(x = hvf.info[, "variable"]) + 1]
  hvf.info <- hvf.info[, c("gmean", "residual_variance")]
  axis.labels <- c('Average Expression', 'Residual Variance')

  plot <- SingleCorPlot(
    data = hvf.info,
    col.by = var.status,
    pt.size = pt.size
  )
  plot <- plot +
    labs(title = NULL, x = axis.labels[1], y = axis.labels[2]) +
    scale_color_manual(
      labels = paste(
        c('Non-variable', 'Variable'),
        'count:',
        table(var.status)
      ),
      values = cols
    )

  plot <- plot + scale_x_log10()

  return(plot)
}

RidgePlot.sce <- function(sce, feature, yaxis) {
  data = data.frame(
    feature = logcounts(sce)[feature, ],
    yaxis = as.character(colData(sce)[, yaxis]),
    row.names = colnames(sce)
  )
  ggplot(data, aes(x = feature, y = yaxis, fill = yaxis)) +
    labs(x = "Expression level", y = yaxis, title = feature, fill = NULL) +
    geom_density_ridges() +
    theme_ridges() +
    theme(legend.position = "none")
  return(plot)
}

cellsProportion <- function(object, groupVar1, groupVar2) {
  if (is(object, "SingleCellExperiment")) {
    cellCounts <- table(
      colData(object)[, groupVar1],
      colData(object)[, groupVar2]
    )
  } else {
    #it is a Seurat object then
    cellCounts <- table(
      object@meta.data[, groupVar1],
      object@meta.data[, groupVar2]
    )
  }

  cellPerc <- sweep(cellCounts, 2, colSums(cellCounts), "/")
  colnames(cellPerc) <- paste0(colnames(cellPerc), "_fraction")
  table <- cbind(cellCounts, cellPerc)
  table <- round(table, digits = 4)
  total <- apply(table, 2, sum)
  table <- cbind(rownames(cellCounts), table)
  rownames(table) <- NULL
  colnames(table)[1] <- groupVar1
  table <- rbind(table, c("Total", total))
  return(table)
}

##' @title Resolve the species from a refBuild string
##' @description Maps a `refBuild` to "Human", "Mouse" or "other". This gates the
##'   CellMarker2/AUCell annotation section, and "other" means NO annotation at
##'   all, silently - so a parsing miss costs a run its annotation with no error.
##'
##'   Matching is on whole path SEGMENTS, not a prefix. The former
##'   `startsWith(refBuild, "Homo_sapiens")` missed three spellings that occur in
##'   live gStore results: an absolute reference path (13 mouse runs, p40923 and
##'   p40924, Jan-Mar 2026), and a bare assembly name (3 runs). Segment matching
##'   also keeps `Mus_minutoides` - which has a reference installed here - from
##'   being read as mouse, which a `grepl("Mus_musculus", ...)` fix would not.
##'
##'   A `Chimera_*` reference is deliberately "other": the matrix mixes two
##'   genomes and neither species' marker sets are valid on it.
##' @param refBuild character(1), e.g. "Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-03"
##' @return "Human", "Mouse" or "other"
getSpecies <- function(refBuild) {
  if (length(refBuild) != 1L || is.na(refBuild) || !nzchar(refBuild)) {
    return("other")
  }
  segments <- setdiff(strsplit(refBuild, "/", fixed = TRUE)[[1]], "")
  if (!length(segments)) {
    return("other")
  }
  ## a mixed-genome reference is not either of its parents
  if (any(startsWith(segments, "Chimera"))) {
    return("other")
  }
  if ("Homo_sapiens" %in% segments) {
    return("Human")
  }
  if ("Mus_musculus" %in% segments) {
    return("Mouse")
  }
  ## refBuild given as a bare assembly, with no genus directory
  if (any(grepl("^GRCh[0-9]", segments))) {
    return("Human")
  }
  if (any(grepl("^GRCm[0-9]", segments))) {
    return("Mouse")
  }
  "other"
}

# Decide whether the Pan-Human Azimuth (CloudAzimuth) step should run, and say
# why not when it should not. This is the ONLY species gate: neither ScSeuratApp
# nor CloudAzimuth itself checks, so without this a mouse job silently ships its
# matrix to the human model. The reference is keyed on HUMAN gene symbols, so on
# a mouse object almost nothing matches (Gpr4 vs GPR4) and what comes back is a
# confident-looking annotation computed from a near-empty matrix.
# Returns list(run, reason); the reason is logged so a skip is visible.
panHumanAzimuthPlan <- function(param) {
  enabled <- ezIsSpecified(param$AzimuthPanHuman) &&
    (isTRUE(param$AzimuthPanHuman) ||
      identical(tolower(as.character(param$AzimuthPanHuman)[1]), "true"))
  if (!enabled) {
    return(list(run = FALSE, reason = "disabled by parameter"))
  }
  if (!ezIsSpecified(param$refBuild)) {
    return(list(
      run = FALSE,
      reason = "no refBuild given, cannot confirm the dataset is human"
    ))
  }
  species <- getSpecies(param$refBuild)
  if (!identical(species, "Human")) {
    return(list(
      run = FALSE,
      reason = paste0(
        "Pan-Human Azimuth is human-only and refBuild resolves to species '",
        species,
        "'"
      )
    ))
  }
  return(list(run = TRUE, reason = "human dataset"))
}

# panhumanpy's label refiner does not return a label when its 8 hierarchy heads
# contradict each other: it returns a SENTINEL (annotate_tools.py:1503 returns the
# string "False"; annotate.py:545 masks None/False/"False"/"false"). CellRanger
# fills that sentinel before writing its CSV (fill_unannotated, annotate.py:548);
# AzimuthAPI::CloudAzimuth does not, so the raw sentinel lands in the delivered
# object. On p42258/o42614 that was 3,504 of 17,777 cells (19.7%) reading "False"
# in azimuth_broad/medium/fine -- the third most common value in all three, which
# exploreSC then lists as a browsable cell type.
#
# We forward-fill from the parent tier (annotate.py:551's fill_unannotated=FALSE
# branch) rather than writing "Not Confidently Annotated": every tier keeps a real
# cell type, just a coarser one, and no cell drops out of a plot legend. A sentinel
# broad falls back to the first component of full_hierarchical_labels, which was
# verified equal to CellRanger's broad_cell_type on 15,163/15,163 cells.
#
# NOT a sentinel: "Unassigned" is a trained reject class (one of 13 level-0 output
# neurons), carrying a median confidence of 0.90 -- it must survive untouched.
fillAzimuthSentinels <- function(meta) {
  isSentinel <- function(x) {
    is.na(x) | as.character(x) %in% c("False", "false", "None", "")
  }
  # order matters: medium inherits the filled broad, fine the filled medium
  hierPart <- function(i) {
    if (!"full_hierarchical_labels" %in% colnames(meta)) {
      return(NULL)
    }
    vapply(
      strsplit(as.character(meta$full_hierarchical_labels), "|", fixed = TRUE),
      function(p) if (length(p) >= i) p[[i]] else NA_character_,
      character(1)
    )
  }
  fillFrom <- function(col, parent) {
    if (!col %in% colnames(meta) || is.null(parent)) {
      return(invisible(NULL))
    }
    bad <- isSentinel(meta[[col]])
    if (any(bad)) {
      meta[[col]] <<- as.character(meta[[col]])
      meta[[col]][bad] <<- as.character(parent)[bad]
    }
  }
  fillFrom("azimuth_broad", hierPart(1))
  fillFrom(
    "azimuth_medium",
    if ("azimuth_broad" %in% colnames(meta)) meta$azimuth_broad else hierPart(1)
  )
  fillFrom(
    "azimuth_fine",
    if ("azimuth_medium" %in% colnames(meta)) meta$azimuth_medium else hierPart(1)
  )
  meta
}

# CellRanger >= 10.1.0 runs the SAME Pan-Human Azimuth model locally and by
# default, writing outs/cell_types/Azimuth/cell_types.csv. Measured on
# p42258/o42614 over an identical 15,163-cell matrix, CellRanger's vendored ONNX
# fork and upstream panhumanpy 1.0.0 agree BIT-FOR-BIT on confidence (Pearson
# r = 1.0000, max |diff| = 0.0000) and on full_hierarchical_labels (15,163/15,163).
# So when the file is there, calling CloudAzimuth recomputes an identical answer
# AND ships the expression matrix to azimuthapi.satijalab.org to do it.
#
# Returns a data.frame keyed by barcode using the column names the report expects,
# or NULL. NULL must mean "fall back to the API" -- never a partial object, which
# is the state the guard at app-ScSeurat.R:807 exists to prevent.
readCellRangerPanHuman <- function(countMatrixPath) {
  if (is.null(countMatrixPath) || length(countMatrixPath) != 1L ||
    is.na(countMatrixPath) || !nzchar(countMatrixPath)) {
    return(NULL)
  }
  # CountMatrix is <outs>/filtered_feature_bc_matrix (dir) or ...matrix.h5 (file);
  # dirname() lands on <outs> either way.
  csv <- file.path(
    dirname(countMatrixPath), "cell_types", "Azimuth", "cell_types.csv"
  )
  if (!file.exists(csv)) {
    return(NULL)
  }
  ann <- try(
    read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE),
    silent = TRUE
  )
  if (inherits(ann, "try-error") || !is.data.frame(ann) || nrow(ann) == 0L) {
    return(NULL)
  }
  needed <- c(
    "barcode", "broad_cell_type", "coarse_cell_type", "fine_cell_type",
    "full_hierarchical_labels", "final_level_softmax_prob"
  )
  if (!all(needed %in% colnames(ann))) {
    return(NULL)
  }
  # final_level_labels is the deepest component of the hierarchy string; verified
  # equal to panhumanpy's own final_level_labels on 15,163/15,163 cells.
  finalLevel <- vapply(
    strsplit(as.character(ann$full_hierarchical_labels), "|", fixed = TRUE),
    function(p) if (length(p)) p[[length(p)]] else NA_character_,
    character(1)
  )
  out <- data.frame(
    full_hierarchical_labels = ann$full_hierarchical_labels,
    final_level_labels = finalLevel,
    final_level_confidence = as.numeric(ann$final_level_softmax_prob),
    azimuth_broad = ann$broad_cell_type,
    azimuth_medium = ann$coarse_cell_type,
    azimuth_fine = ann$fine_cell_type,
    azimuth_label = finalLevel,
    stringsAsFactors = FALSE
  )
  # full_consistent_hierarchy is deliberately NOT reconstructed: deriving it from
  # the "Not Confidently Annotated" label matches the real flag on only
  # 15,055/15,163 cells, and multiOmicsUtils.R only uses it to EXCLUDE a column
  # from group.by lists, so absent is both honest and harmless.
  rownames(out) <- ann$barcode
  out
}

geneMeansCluster <- function(object) {
  if (is(object, "SingleCellExperiment")) {
    tr_cnts <- expm1(logcounts(object))
    group = object$ident
  } else {
    tr_cnts <- expm1(GetAssayData(object, layer = "data", assay = "SCT"))
    group = Idents(object)
  }
  geneMeans <- rowsum(DelayedArray::t(tr_cnts), group = group)
  geneMeans <- sweep(
    geneMeans,
    1,
    STATS = table(group)[rownames(geneMeans)],
    FUN = "/"
  )
  geneMeans <- log1p(t(geneMeans))
  colnames(geneMeans) <- paste("cluster", colnames(geneMeans), sep = "_")
  return(geneMeans)
}

cellsLabelsWithAUC <- function(
  counts,
  species,
  tissue,
  minGsSize = 3,
  BPPARAM = NULL
) {
  if (species == "other") {
    return(NULL)
  }
  # Return NULL if tissue is not specified
  if (!ezIsSpecified(tissue)) {
    futile.logger::flog.info(
      "No tissue specified for CellMarker2, skipping AUC annotation"
    )
    return(NULL)
  }
  geneSets <- createCellMarker2_GeneSets(species, tissue, minGsSize)
  if (is.null(geneSets) || length(geneSets) == 0) {
    return(NULL)
  }
  cells_rankings <- AUCell_buildRankings(
    counts,
    plotStats = FALSE,
    BPPARAM = BPPARAM,
    splitByBlocks = TRUE
  )
  cells_AUC <- AUCell_calcAUC(
    geneSets,
    cells_rankings,
    verbose = FALSE,
    nCores = ifelse(is.null(BPPARAM), 1, BPPARAM$workers)
  )
  return(cells_AUC)
}


## old function uses CellMarker v1, gene sets
createGeneSets <- function(species, tissue) {
  tissue <- unlist(strsplit(tissue, ","))
  cell_markers <- read.table(
    "/srv/GT/databases/scGeneSets/all_cell_markers.txt",
    sep = "\t",
    header = TRUE
  )
  cell_markers <- cell_markers[
    cell_markers$speciesType == species &
      cell_markers$tissueType %in% tissue,
  ]
  if (nrow(cell_markers) == 0) {
    #stop(sprintf("No cell markers found for %s: %s", species, paste(tissue, collapse=", ")))
    warning(sprintf(
      "No cell markers found for %s: %s",
      species,
      paste(tissue, collapse = ", ")
    ))
    return(NULL)
  }
  geneSetList <- strsplit(cell_markers$geneSymbol, ",")
  geneSetList <- lapply(geneSetList, function(gs) {
    gs <- gs[!is.na(gs)]
    gs <- gsub("^ ", "", gsub(" $", "", gs))
    gs <- gsub("[", "", gs, fixed = TRUE)
    gs <- gsub("]", "", gs, fixed = TRUE)
    gs <- gsub("11-Sep", "SEPTIN11", gs)
    gs <- setdiff(gs, c("NA", ""))
  })
  ## merge the genesets from the same cell type
  geneSetArray = tapply(
    geneSetList,
    cell_markers$cellName,
    function(x) {
      unique(unlist(x))
    },
    simplify = FALSE
  )
  ## conver the array  returned by tapply to a list
  geneSetList = lapply(geneSetArray, function(gs) {
    gs
  })
  return(geneSetList)
}


createCellMarker2_GeneSets <- function(species, tissue, minGsSize = 3) {
  tissue <- unlist(strsplit(tissue, ","))
  cell_markers <- ezRead.table(
    "/srv/GT/databases/scGeneSets/CellMarker_2.0-2023-09-27/Cell_marker_All.txt",
    row.names = NULL
  )
  cell_markers <- cell_markers[
    cell_markers$species == species &
      cell_markers$tissue_class %in% tissue,
  ]
  if (nrow(cell_markers) == 0) {
    #stop(sprintf("No cell markers found for %s: %s", species, paste(tissue, collapse=", ")))
    warning(sprintf(
      "No cell markers found for %s: %s",
      species,
      paste(tissue, collapse = ", ")
    ))
    return(NULL)
  }
  geneSetList <- split(cell_markers$Symbol, cell_markers$cell_name)
  geneSetList <- geneSetList[sapply(geneSetList, length) >= minGsSize]
  return(geneSetList)
}


cellsLabelsWithSingleR <- function(
  logCounts,
  current_clusters,
  refDataName,
  BPPARAM = SerialParam()
) {
  if (!ezIsSpecified(refDataName) || refDataName == "none") {
    return(NULL)
  }
  library(SingleR)
  singlerResultsList <- list()
  for (nm in refDataName) {
    ref <- eval(parse(text = paste0('celldex::', nm, "()")))
    singlerResultsList[[nm]] <- list()
    singlerResultsList[[nm]][["single.fine"]] <- SingleR(
      test = logCounts,
      ref = ref,
      labels = ref$label.fine,
      BPPARAM = BPPARAM
    )
    singlerResultsList[[nm]][["cluster.fine"]] <- SingleR(
      test = logCounts,
      ref = ref,
      labels = ref$label.fine,
      clusters = current_clusters,
      BPPARAM = BPPARAM
    )
  }
  return(singlerResultsList)
}

filterCellsAndGenes <- function(object, param) {
  UseMethod("filterCellsAndGenes", object)
}

filterCellsAndGenes.Seurat <- function(scData, param) {
  library(scater)
  library(Seurat)

  # Cells filtering
  scData <- PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
  scData <- PercentageFeatureSet(
    scData,
    "(?i)^RPS|^RPL",
    col.name = "percent_riboprot"
  )
  if (grepl("Spatial", param$appName)) {
    assay <- "Spatial"
    att_nCounts <- "nCount_Spatial"
    att_nGenes <- "nFeature_Spatial"
  } else {
    att_nCounts <- "nCount_RNA"
    att_nGenes <- "nFeature_RNA"
    assay <- "RNA"
  }

  if (is.na(param$nreads)) {
    if (assay == 'Spatial') {
      qc.lib <- scData@meta.data[, 'nCount_Spatial_SpotSweeper_outliers']
    } else {
      qc.lib <- isOutlier(
        scData@meta.data[, att_nCounts],
        log = TRUE,
        nmads = param$nmad,
        type = "lower"
      )
    }
  } else {
    qc.lib <- scData@meta.data[, att_nCounts] < param$nreads
  }

  if (is.na(param$ngenes)) {
    if (assay == 'Spatial') {
      qc.nexprs <- scData@meta.data[, 'nFeature_Spatial_SpotSweeper_outliers']
    } else {
      qc.nexprs <- isOutlier(
        scData@meta.data[, att_nGenes],
        nmads = param$nmad,
        log = TRUE,
        type = "lower"
      )
    }
  } else {
    qc.nexprs <- scData@meta.data[, att_nGenes] < param$ngenes
  }

  if (is.na(param$perc_mito)) {
    if (assay == 'Spatial') {
      qc.mito <- scData@meta.data[, "percent_mito_SpotSweeper_outliers"]
    } else {
      qc.mito <- isOutlier(
        scData@meta.data[, "percent_mito"],
        nmads = param$nmad,
        log = TRUE,
        type = "lower"
      )
    }
  } else {
    qc.mito <- scData@meta.data[, "percent_mito"] > param$perc_mito
  }

  if (is.na(param$perc_ribo)) {
    if (assay == 'Spatial') {
      qc.ribo <- rep(FALSE, nrow(scData@meta.data))
    } else {
      qc.ribo <- isOutlier(
        scData@meta.data[, "percent_riboprot"],
        nmads = param$nmad,
        type = "higher"
      )
    }
  } else {
    qc.ribo <- scData@meta.data[, "percent_riboprot"] > param$perc_ribo
  }

  discard <- qc.lib | qc.nexprs | qc.mito | qc.ribo
  scData$discard <- discard
  scData$qc.lib <- qc.lib
  scData$qc.nexprs <- qc.nexprs
  scData$qc.mito <- qc.mito
  scData$qc.ribo <- qc.ribo
  scData.unfiltered <- scData
  if (any(discard)) {
    scData <- scData[, -which(discard)]
  }

  # Genes filtering
  ## remove low expressed genes
  num.cells <- param$cellsFraction * ncol(scData) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  cellsPerGene <- Matrix::rowSums(
    GetAssayData(scData, layer = "counts") >= param$nUMIs
  )
  is.expressed <- cellsPerGene >= num.cells
  cellsPerGeneFraction <- data.frame(
    frac = cellsPerGene / ncol(scData),
    row.names = rownames(cellsPerGene)
  )
  scData <- scData[is.expressed, ]
  return(list(
    scData.unfiltered = scData.unfiltered,
    scData = scData,
    cellsPerGeneFraction = cellsPerGeneFraction
  ))
}

##' @title FindClusters without Seurat's per-singleton loop
##' @description Drop-in for `Seurat::FindClusters()` (one resolution). Seurat's
##'   `GroupSingletons` loops over singletons x clusters with a name-indexed
##'   sparse subset per pair; on a degenerate SNN graph (p26168 TMA4 8 um:
##'   14,007 singletons, 43 clusters) that ran 32 h on one core. Here Seurat
##'   labels them "singleton" and [reassignSingletons()] applies the same rule
##'   in one sparse product.
##' @param object Seurat object with `graph.name` computed.
##' @param resolution numeric(1).
##' @param graph.name,cluster.name as in `FindClusters()`.
##' @param ... passed to `FindClusters()`.
##' @return the object, with the cluster column, `Idents` and (when Seurat set
##'   it) `seurat_clusters` free of the "singleton" label.
findClustersFast <- function(object, resolution, graph.name = NULL,
                             cluster.name = NULL, ...) {
  graph.name <- graph.name %||% paste0(DefaultAssay(object), "_snn")
  cluster.name <- cluster.name %||% paste0(graph.name, "_res.", resolution)
  ## Seurat 5.5.1's parallel branch (nbrOfWorkers() > 1) drops
  ## group.singletons and runs the slow loop anyway; one resolution gains
  ## nothing from it.
  oplan <- future::plan("sequential")
  on.exit(future::plan(oplan), add = TRUE)
  object <- FindClusters(object, resolution = resolution,
                         graph.name = graph.name, cluster.name = cluster.name,
                         group.singletons = FALSE, ...)
  ids <- setNames(as.character(object[[cluster.name, drop = TRUE]]),
                  colnames(object))
  single <- !is.na(ids) & ids == "singleton"
  if (!any(single)) {
    return(object)
  }
  futile.logger::flog.info(
    "%s: reassigning %d singletons (a large count means a near-random SNN graph)",
    cluster.name, sum(single)
  )
  ids[!is.na(ids)] <- reassignSingletons(ids[!is.na(ids)], object[[graph.name]])
  levs <- as.character(sort(as.integer(unique(na.omit(ids)))))
  newIds <- setNames(factor(ids, levels = levs), colnames(object))
  hadSeuratClusters <- identical(as.character(object$seurat_clusters),
                                 as.character(object[[cluster.name, drop = TRUE]]))
  object[[cluster.name]] <- newIds
  Idents(object) <- newIds
  if (hadSeuratClusters) {
    object$seurat_clusters <- newIds
  }
  object
}

##' @title Assign singleton cells to their best-connected cluster
##' @description The rule of Seurat's `GroupSingletons`: each singleton joins
##'   the cluster with the highest mean SNN weight to it, and a tie is broken
##'   exactly as Seurat does (`set.seed(1); sample(tied, 1)`, clusters in order
##'   of first appearance). That tie-break matters: on p26168 TMA4 all 14,007
##'   singletons had zero weight to every cluster (their edges only reach other
##'   singletons), so each was a 43-way tie and all went to one cluster. Seurat
##'   grows clusters as it goes, which can only change a near-tied call.
##' @param ids named character; singletons carry the label "singleton".
##' @param snn cells x cells SNN graph whose dimnames include `names(ids)`.
##' @return `ids` with every "singleton" replaced by a cluster label.
reassignSingletons <- function(ids, snn) {
  single <- ids == "singleton"
  if (!any(single) || all(single)) {
    return(ids)
  }
  cl <- factor(ids[!single], levels = unique(ids[!single]))
  member <- Matrix::sparseMatrix(
    i = match(names(cl), rownames(snn)), j = as.integer(cl), x = 1,
    dims = c(nrow(snn), nlevels(cl))
  )
  conn <- as.matrix(snn[names(ids)[single], , drop = FALSE] %*% member)
  conn <- sweep(conn, 2, tabulate(as.integer(cl), nlevels(cl)), "/")
  best <- levels(cl)[max.col(conn, ties.method = "first")]
  isMax <- conn == apply(conn, 1, max)
  for (r in which(rowSums(isMax) > 1)) {
    best[r] <- withr::with_seed(1, sample(levels(cl)[isMax[r, ]], 1))
  }
  ids[single] <- best
  ids
}
