###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

getCellCycle <- function(counts, refBuild){
  require(scran)
  # The training data is only available for Hsap and Mmus Ensembl
  if(startsWith(refBuild, "Homo_sapiens")){
    species <- "human"
    hasTrainData <- TRUE
  }else if(startsWith(refBuild, "Mus_musculus")){
    species <- "mouse"
    hasTrainData <- TRUE
  }else{
    hasTrainData <- FALSE
  }
  if(isTRUE(hasTrainData)){
    trainData <- readRDS(system.file("exdata",
                                     paste0(species, "_cycle_markers.rds"), 
                                     package = "scran", mustWork=TRUE))
    cellCycleData <- cyclone(counts, trainData)
    cellPhase <- tibble(Name = colnames(counts),
                        Phase = cellCycleData$phases)
    cellPhase <- bind_cols(cellPhase, cellCycleData$scores)
  }else{
    cellPhase <- tibble(Name = colnames(counts),
                        Phase = NA)
  }
  return(cellPhase)
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

cellsProportion <- function(sce, groupVar1, groupVar2){
  library(tidyverse)
  toTable <- tibble(names(summary(colData(sce)[,groupVar1])))  %>%
    `colnames<-`(groupVar1) 
  cellCountsByVar <- tibble(colData(sce)[,groupVar2],
                            as.character(colData(sce)[,groupVar1])) %>%
    `colnames<-`(c(groupVar2, groupVar1)) %>%
    group_by_at(c(groupVar2, groupVar1)) %>% summarise(n()) %>%
    spread(groupVar2, `n()`, fill=0)
  cellPercByVar2 <- select(cellCountsByVar, -groupVar1) %>%
    as.matrix()
  rownames(cellPercByVar2) <- cellCountsByVar[[groupVar1]]
  cellPercByVar2 <- sweep(cellPercByVar2, 2, colSums(cellPercByVar2), "/")
  colnames(cellPercByVar2) <- paste0(colnames(cellPercByVar2), "_fraction")
  toTable <- left_join(toTable, cellCountsByVar, by=groupVar1) %>%
    left_join(as_tibble(cellPercByVar2, rownames=groupVar1), by=groupVar1)
  total <- bind_cols("Total", summarise_at(toTable, setdiff(colnames(toTable), groupVar1),sum))
  colnames(total)[1] <- groupVar1
  toTable <- bind_rows(toTable,total)
  return(toTable)
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

geneMeansCluster <- function(sce) {
tr_cnts <- expm1(logcounts(sce))
means <- rowsum(t(as.matrix(tr_cnts)), group=colData(sce)[,"ident"])
means <- sweep(means, 1, STATS=table(colData(sce)[,"ident"])[rownames(means)], FUN="/")
means <- log1p(t(means))
colnames(means) <- paste("cluster", colnames(means), sep="_")
return(means)
}

cellsLabelsWithAUC <- function(scData, species, tissue, minGsSize = 3) {
  library(AUCell)
  if (species == "other")
    return(NULL)
  geneSets <- createGeneSets(species, tissue)
  expressionMatrix <- GetAssayData(scData, slot = "counts")
  cells_rankings <- AUCell_buildRankings(expressionMatrix, plotStats=FALSE)
  cells_AUC <- tryCatch({AUCell_calcAUC(geneSets[sapply(geneSets, length) >= minGsSize], cells_rankings, verbose = FALSE)},error = function(e) NULL)
  return(cells_AUC)
}


createGeneSets <- function(species, tissue) {
  tissue <- unlist(strsplit(tissue, ","))
  cell_markers <- read.table("/srv/GT/databases/scGeneSets/all_cell_markers.txt", sep = "\t", header = TRUE)
  cell_markers <- cell_markers[cell_markers$speciesType == species & 
                                 cell_markers$tissueType %in% tissue, ]
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

cellsLabelsWithSingleR <- function(counts, current_clusters, species) {
  library(SingleR)
  if(species == "Human"){
    reference <- HumanPrimaryCellAtlasData()
    singler.results.single <- SingleR(test = counts, ref = reference, 
                                      labels = reference$label.main, method="single", de.method = "wilcox")
    singler.results.cluster <- SingleR(test = counts, ref = reference, 
                                       labels = reference$label.main, method="cluster", clusters=current_clusters, de.method = "wilcox")
  }else {
    reference <- celldex::MouseRNAseqData()
    singler.results.single <- SingleR(test = counts, ref = reference, labels = reference$label.main)
    singler.results.cluster <- SingleR(test = counts, ref = reference, labels = reference$label.main, method="cluster", clusters=current_clusters)
  }
  
  return(list(singler.results.single=singler.results.single, singler.results.cluster=singler.results.cluster))
}