###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch




addCellCycleToSCE <- function(sce, refBuild, BPPARAM){
  counts <- counts(sce)
  rownames(counts) <- rowData(sce)$ID
  cellPhase <- getCellCycle(counts, refBuild, BPPARAM)
  if (!is.null(cellPhase)){
    colData(sce)$CellCycle <- cellPhase$Phase
    colData(sce)$CellCycleG1 <- cellPhase$G1
    colData(sce)$CellCycleS <- cellPhase$S
    colData(sce)$CellCycleG2M <- cellPhase$G2M
  }
  return(sce)
}

addCellCycleToSeurat <- function(scData, refBuild, BPPARAM, assay = "RNA"){
  counts <- GetAssayData(scData, layer="counts", assay = assay)
  metaFeatures <- scData[[assay]]@meta.data
  if ("gene_id" %in% names(metaFeatures)) {
    rownames(counts) <- metaFeatures$gene_id
  } else {
    rownames(counts) <- metaFeatures$ensemblID
  }
  cellPhase <- getCellCycle(counts, refBuild, BPPARAM)
  if (!is.null(cellPhase)){
    cellcycleInfo = data.frame(CellCycle = cellPhase$Phase, CellCycleG1 = cellPhase$G1, CellCycleS = cellPhase$S, CellCycleG2M = cellPhase$G2M, row.names = colnames(scData))
    scData <- AddMetaData(scData, metadata = cellcycleInfo)
  }
  return(scData)
}


getCellCycle <- function(counts, refBuild, BPPARAM){
  require(scran)

  species <- sub("\\/.*", "", refBuild)
  trainDataFile <- switch(species,
                      Homo_sapiens=system.file("exdata","human_cycle_markers.rds", 
                                               package = "scran"),
                      Mus_musculus=system.file("exdata","mouse_cycle_markers.rds", 
                                               package = "scran"))
  if (is.null(trainDataFile) || !file.exists(trainDataFile)){
    return(NULL)
  } else {
    trainData <- readRDS(trainDataFile)
    cellCycleData <- cyclone(counts, trainData, BPPARAM = BPPARAM)
    cellPhase <- tibble(Name = colnames(counts),
                        Phase = cellCycleData$phases)
    cellPhase <- bind_cols(cellPhase, cellCycleData$scores)
    return(cellPhase)
  }
}

getPerplexity <- function(n){
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
    plot <- plot + stat_density2d(
      mapping = aes(fill = ..density.. ^ 0.25),
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
    plot <- plot + geom_point(
      mapping = aes_string(color = 'colors'),
      position = 'jitter',
      size = pt.size
    )
  } else {
    plot <- plot + geom_point(position = 'jitter', size = pt.size)
  }
  if (!is.null(x = cols)) {
    cols.scale <- if (length(x = cols) == 1 && cols %in% rownames(x = brewer.pal.info)) {
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
    plot <- plot + geom_smooth(
      mapping = aes_string(x = names.plot[1], y = names.plot[2]),
      method = 'loess',
      span = span
    )
  }
  return(plot)
}


VariableFeaturePlot_mod <- function(sce, cols = c('black', 'red'), pt.size = 1, log = NULL) {
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
      labels = paste(c('Non-variable', 'Variable'), 'count:', table(var.status)),
      values = cols
    )
  
  plot <- plot + scale_x_log10()
  
  return(plot)
}

RidgePlot.sce <- function(sce, feature, yaxis) {
   data= data.frame(feature = logcounts(sce)[feature,], yaxis = as.character(colData(sce)[,yaxis]), row.names = colnames(sce))
   ggplot(data, aes(x = feature, y = yaxis, fill = yaxis)) +
    labs(x = "Expression level", y = yaxis, title = feature, fill = NULL) +
    geom_density_ridges() +
    theme_ridges() + 
    theme(legend.position = "none")
  return(plot)
}

cellsProportion <- function(object, groupVar1, groupVar2) {
  if(is(object, "SingleCellExperiment"))
     cellCounts <- table(colData(object)[,groupVar1], colData(object)[,groupVar2])
  else  #it is a Seurat object then
    cellCounts <- table(object@meta.data[,groupVar1], object@meta.data[,groupVar2])
  
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

getSpecies <- function(refBuild) {
  if(startsWith(refBuild, "Homo_sapiens")){
    species <- "Human"
  }else if(startsWith(refBuild, "Mus_musculus")){
    species <- "Mouse"
  } else {
  species <- "other"
  }
  return(species)
}

geneMeansCluster <- function(object) {
  if(is(object, "SingleCellExperiment")) {
    tr_cnts <- expm1(logcounts(object))
    group=object$ident
 } else {
    tr_cnts <- expm1(GetAssayData(object, layer="data", assay = "SCT"))
    group=Idents(object)
 }
  geneMeans <- rowsum(DelayedArray::t(tr_cnts), group=group)
  geneMeans <- sweep(geneMeans, 1, STATS=table(group)[rownames(geneMeans)], FUN="/")
  geneMeans <- log1p(t(geneMeans))
  colnames(geneMeans) <- paste("cluster", colnames(geneMeans), sep="_")
return(geneMeans)
}

cellsLabelsWithAUC <- function(counts, species, tissue, minGsSize = 3, BPPARAM=NULL) {
  if (species == "other"){
    return(NULL)
  }
  geneSets <- createCellMarker2_GeneSets(species, tissue, minGsSize)
  if(is.null(geneSets) || length(geneSets) == 0){ 
    return(NULL)
  }
  cells_rankings <- AUCell_buildRankings(counts, plotStats=FALSE, BPPARAM=BPPARAM, splitByBlocks=TRUE)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, verbose = FALSE, 
                   nCores = ifelse(is.null(BPPARAM), 1, BPPARAM$workers))
  return(cells_AUC)
}


## old function uses CellMarker v1, gene sets
createGeneSets <- function(species, tissue) {
  tissue <- unlist(strsplit(tissue, ","))
  cell_markers <- read.table("/srv/GT/databases/scGeneSets/all_cell_markers.txt", sep = "\t", header = TRUE)
  cell_markers <- cell_markers[cell_markers$speciesType == species &
                                 cell_markers$tissueType %in% tissue, ]
  if (nrow(cell_markers) == 0) {
    #stop(sprintf("No cell markers found for %s: %s", species, paste(tissue, collapse=", ")))
    warning(sprintf("No cell markers found for %s: %s", species, paste(tissue, collapse=", ")))
    return(NULL)
  }
  geneSetList <- strsplit(cell_markers$geneSymbol, ",")
  geneSetList <- lapply(geneSetList, function(gs){
    gs <- gs[!is.na(gs)]
    gs <- gsub("^ ", "", gsub(" $", "", gs))
    gs <- gsub("[", "", gs, fixed = TRUE)
    gs <- gsub("]", "", gs, fixed = TRUE)
    gs <- gsub("11-Sep", "SEPTIN11", gs)
    gs <- setdiff(gs, c("NA", ""))
  })
  ## merge the genesets from the same cell type
  geneSetArray = tapply(geneSetList, cell_markers$cellName, 
                        function(x){unique(unlist(x))}, simplify = FALSE)
  ## conver the array  returned by tapply to a list
  geneSetList = lapply(geneSetArray, function(gs){gs})
  return(geneSetList)
}


createCellMarker2_GeneSets <- function(species, tissue, minGsSize=3) {
  tissue <- unlist(strsplit(tissue, ","))
  cell_markers <- ezRead.table("/srv/GT/databases/scGeneSets/CellMarker_2.0-2023-09-27/Cell_marker_All.txt", row.names = NULL)
  cell_markers <- cell_markers[cell_markers$species == species & 
                                 cell_markers$tissue_class %in% tissue, ]
  if (nrow(cell_markers) == 0) {
    #stop(sprintf("No cell markers found for %s: %s", species, paste(tissue, collapse=", ")))
    warning(sprintf("No cell markers found for %s: %s", species, paste(tissue, collapse=", ")))
    return(NULL)
  }
  geneSetList <- split(cell_markers$Symbol, cell_markers$cell_name)
  geneSetList <- geneSetList[sapply(geneSetList, length) >= minGsSize]
  return(geneSetList)
}


cellsLabelsWithSingleR <- function(logCounts, current_clusters, refDataName, BPPARAM = SerialParam()) {
  if (!ezIsSpecified(refDataName) || refDataName == "none"){
    return(NULL)
  }
  library(SingleR)
  singlerResultsList <- list()
  for (nm in refDataName){
    ref <- eval(parse(text=paste0('celldex::', nm, "()")))
    singlerResultsList[[nm]] <- list()
    singlerResultsList[[nm]][["single.fine"]] <- SingleR(test = logCounts, ref = ref, labels = ref$label.fine, BPPARAM = BPPARAM)
    singlerResultsList[[nm]][["cluster.fine"]] <- SingleR(test = logCounts, ref = ref, labels = ref$label.fine, clusters = current_clusters, BPPARAM = BPPARAM)
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
  scData <- PercentageFeatureSet(scData, "(?i)^RPS|^RPL", col.name = "percent_riboprot")
  if(grepl("Spatial", param$appName)) {
    assay <- "Spatial"
    att_nCounts <- "nCount_Spatial"
    att_nGenes <- "nFeature_Spatial"
  } else {
    att_nCounts <- "nCount_RNA"
    att_nGenes <- "nFeature_RNA"
    assay <- "RNA"
  }
  
  if (is.na(param$nreads)) {
    qc.lib <- isOutlier(scData@meta.data[,att_nCounts], log = TRUE, nmads = param$nmad, type = "lower")
  } else {
    qc.lib <- scData@meta.data[,att_nCounts] < param$nreads
  }
  if (is.na(param$ngenes)) {
    qc.nexprs <- isOutlier(scData@meta.data[,att_nGenes], nmads = param$nmad, log = TRUE, type = "lower")
  } else {
    qc.nexprs <- scData@meta.data[,att_nGenes] < param$ngenes
  }
  if (is.na(param$perc_mito)) {
    qc.mito <- isOutlier(scData@meta.data[,"percent_mito"], nmads = param$nmad, type = "higher")
  } else {
    qc.mito <- scData@meta.data[,"percent_mito"] > param$perc_mito
  }
  if (is.na(param$perc_ribo )) {
    qc.ribo <- isOutlier(scData@meta.data[,"percent_riboprot"], nmads = param$nmad, type = "higher")
  } else {
    qc.ribo <- scData@meta.data[,"percent_riboprot"] > param$perc_ribo
  }
  
  discard <- qc.lib | qc.nexprs | qc.mito | qc.ribo
  scData$discard <- discard
  scData$qc.lib <- qc.lib
  scData$qc.nexprs <- qc.nexprs
  scData$qc.mito <- qc.mito
  scData$qc.ribo <- qc.ribo
  scData.unfiltered <- scData
  if(any(discard))
    scData <- scData[, -which(discard)]
  
  # Genes filtering
  ## remove low expressed genes
  num.cells <- param$cellsFraction * ncol(scData) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  cellsPerGene <- Matrix::rowSums(GetAssayData(scData, layer="counts") >= param$nUMIs)
  is.expressed <- cellsPerGene >= num.cells
  cellsPerGeneFraction <- data.frame(frac = cellsPerGene/ncol(scData), row.names = rownames(cellsPerGene))
  scData <- scData[is.expressed,]
  return(list(scData.unfiltered = scData.unfiltered, scData = scData, cellsPerGeneFraction = cellsPerGeneFraction))
}

filterCellsAndGenes.SingleCellExperiment <- function(sce, param) {
  library(scater)
  library(Matrix)
  
  # Cells filtering
  mito.genes <- grep("^MT.", rownames(sce), ignore.case = TRUE)
  ribo.genes <- grep("^RPS|^RPL", rownames(sce), ignore.case = TRUE)
  
  sce <- addPerCellQC(sce, subsets = list(Mito = mito.genes, Ribo = ribo.genes))
  
  if (is.na(param$nreads)) {
    qc.lib <- isOutlier(sce$sum, log = TRUE, nmads = param$nmad, type = "lower")
  } else {
    qc.lib <- sce$sum < as.double(param$nreads)
  }
  if (is.na(param$ngenes)) {
    qc.nexprs <- isOutlier(sce$detected, nmads = param$nmad, log = TRUE, type = "lower")
  } else {
    qc.nexprs <- sce$detected < as.double(param$ngenes)
  }
  if (is.na(param$perc_mito)) {
    qc.mito <- isOutlier(sce$subsets_Mito_percent, nmads = param$nmad, type = "higher")
  } else {
    qc.mito <- sce$subsets_Mito_percent > as.double(param$perc_mito)
  }
  
  if (is.na(param$perc_ribo)) {
    qc.ribo <- isOutlier(sce$subsets_Ribo_percent, nmads = param$nmad, type = "higher")
  } else {
    qc.ribo <- sce$subsets_Ribo_percent > as.double(param$perc_ribo)
  }
  
  discard <- qc.lib | qc.nexprs | qc.mito | qc.ribo
  sce$discard <- discard
  sce$qc.lib <- qc.lib
  sce$qc.nexprs <- qc.nexprs
  sce$qc.mito <- qc.mito
  sce$qc.ribo <- qc.ribo
  sce.unfiltered <- sce
  sce <- sce[, !discard]
  
  # Genes filtering
  num.cells <- param$cellsFraction * ncol(sce) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  is.expressed <- Matrix::rowSums(counts(sce) >= param$nUMIs) >= num.cells
  sce <- sce[is.expressed, ]
  rowData(sce.unfiltered)$is.expressed <- is.expressed
  
  return(list(sce.unfiltered = sce.unfiltered, sce = sce))
}
