###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

getCellCycle <- function(counts, refBuild){
  require(scran)
  require(tibble)
  require(dplyr)
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

RidgePlot.sce <- function(sce, feature) {
  data= data.frame(logcounts(sce)[feature,], cluster = colData(sce)[,"ident"], row.names = colnames(sce))
  if(grepl("-", feature))
    feature <- gsub("-", "", feature)
  colnames(data)[1] = feature
  
  plot <- ggplot(data, aes_string(x = feature, y = "cluster", fill = "cluster")) +
    labs(x = "Expression level", y = "Identity", title = feature, fill = NULL) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))+
    geom_density_ridges() +
    theme_ridges() +
    theme(legend.position = "none")
  return(plot)
}
