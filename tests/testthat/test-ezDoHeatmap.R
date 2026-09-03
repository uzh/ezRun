## Regression tests for the marker heatmap that lost its gene names.
##
## Seurat::DoHeatmap draws one raster column per cell, twice: the expression
## body (raster = TRUE) and the group bar (always), plus ceiling(nCells * 0.0025)
## blank columns per group. Cairo refuses a raster wider than 32767 px and then
## silently stops drawing for the rest of the page:
##   body too wide                  -> all-white PNG            (p39836/o43013, 2026-09-02)
##   raster = FALSE, bar too wide   -> tiles, then nothing: no gene names, no group
##                                     bar, no cluster labels, no legend (p39836, 2026-09-03)
## ezDoHeatmap caps the number of cells so the widest raster stays below the limit.

suppressPackageStartupMessages(library(Seurat))

## Count distinct bytes in the decompressed IDAT stream. A single-colour PNG
## decompresses to one repeated byte, so 1 == blank.
pngDistinctBytes <- function(path) {
  raw <- readBin(path, "raw", file.size(path))
  i <- 9L
  idat <- raw(0)
  while (i < length(raw)) {
    len <- sum(as.integer(raw[i:(i + 3L)]) * c(2^24, 2^16, 2^8, 1))
    if (rawToChar(raw[(i + 4L):(i + 7L)]) == "IDAT") idat <- c(idat, raw[(i + 8L):(i + 7L + len)])
    i <- i + 12L + len
  }
  length(unique(memDecompress(idat, "gzip")))
}

## Ink in the strip left of the panel, where ggplot draws the y-axis (gene) labels.
## 12-character gene names at 192 dpi push the panel edge past 200 px, so the
## first 120 px hold nothing but label glyphs.
leftMarginInk <- function(path, width = 120L) {
  img <- png::readPNG(path)
  sum(img[, seq_len(width), 1:3] < 0.9)
}

renderToPng <- function(plot) {
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 2304, height = 2880, res = 192)  # knitr's device settings
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  path
}

makeObject <- function(nCells, nGenes = 3L, nGroups = 15L) {
  set.seed(38)
  counts <- matrix(rpois(nGenes * nCells, 5), nrow = nGenes,
                   dimnames = list(sprintf("MARKERGENE%02d", seq_len(nGenes)),
                                   paste0("Cell", seq_len(nCells))))
  obj <- Seurat::CreateSeuratObject(as(counts, "dgCMatrix"))
  obj <- Seurat::SetAssayData(obj, layer = "scale.data",
                              new.data = scale(t(as.matrix(counts))) |> t())
  Seurat::Idents(obj) <- factor(rep_len(seq_len(nGroups), nCells))
  obj
}

## The group bar is DoHeatmap's annotation_raster: one column per drawn cell plus padding.
groupBarWidth <- function(plot) {
  isBar <- vapply(plot$layers, function(l) inherits(l$geom, "GeomRasterAnn"), logical(1))
  ncol(plot$layers[[which(isBar)[1]]]$geom_params$raster)
}

## 32562 cells + 15 groups * ceiling(32562 * 0.0025) = 33792 columns > 32767.
test_that("the blankness check can actually detect a blank PNG", {
  obj <- makeObject(32562L)
  expect_equal(pngDistinctBytes(renderToPng(
    DoHeatmap(obj, features = rownames(obj)))), 1L)
})

test_that("the gene-name check can detect the labels going missing", {
  skip_if_not_installed("png")
  obj <- makeObject(32562L)
  ## raster = FALSE draws the tiles, then the oversized group bar kills the device
  expect_equal(leftMarginInk(renderToPng(
    DoHeatmap(obj, features = rownames(obj), raster = FALSE))), 0L)
  ## the same check does see the labels when the bar fits
  expect_gt(leftMarginInk(renderToPng(
    DoHeatmap(makeObject(1000L), features = rownames(obj)))), 0L)
})

test_that("ezDoHeatmap keeps gene names, group bar and legend above the limit", {
  skip_if_not_installed("png")
  obj <- makeObject(32562L)
  p <- ezDoHeatmap(obj, features = rownames(obj))
  expect_lt(groupBarWidth(p), 32767L)
  path <- renderToPng(p)
  expect_gt(pngDistinctBytes(path), 1L)
  expect_gt(leftMarginInk(path), 0L)
})

test_that("ezDoHeatmap draws every cell and keeps the raster path below the limit", {
  obj <- makeObject(5000L)
  p <- ezDoHeatmap(obj, features = rownames(obj))
  expect_equal(groupBarWidth(p), 5000L + 15L * ceiling(5000 * 0.0025))
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomRaster")
})

test_that("an explicit cells argument is honoured", {
  obj <- makeObject(40000L)
  cells <- colnames(obj)[1:2000]
  expect_equal(groupBarWidth(ezDoHeatmap(obj, features = rownames(obj), cells = cells)),
               2000L + 15L * ceiling(2000 * 0.0025))
})
