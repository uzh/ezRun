## Regression test for the silently blank marker heatmap (p39836/o43013, 2026-09-02).
## Seurat's DoHeatmap(raster = TRUE) draws one raster column per cell; cairo drops
## a raster wider than 32767 px with no error and no warning, so the report
## embedded an all-white PNG. ezDoHeatmap must switch to the vector renderer.

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

renderToPng <- function(plot) {
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 2304, height = 2880)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  path
}

makeObject <- function(nCells, nGenes = 3L, nGroups = 15L) {
  set.seed(38)
  counts <- matrix(rpois(nGenes * nCells, 5), nrow = nGenes,
                   dimnames = list(paste0("Gene", seq_len(nGenes)),
                                   paste0("Cell", seq_len(nCells))))
  obj <- Seurat::CreateSeuratObject(as(counts, "dgCMatrix"))
  obj <- Seurat::SetAssayData(obj, layer = "scale.data",
                              new.data = scale(t(as.matrix(counts))) |> t())
  Seurat::Idents(obj) <- factor(rep_len(seq_len(nGroups), nCells))
  obj
}

test_that("the blankness check can actually detect a blank PNG", {
  ## Positive control: bare DoHeatmap above the cairo raster limit IS blank.
  ## 32562 cells + 15 groups * ceiling(32562*0.0025) = 33792 > 32767.
  obj <- makeObject(32562L)
  expect_equal(pngDistinctBytes(renderToPng(
    Seurat::DoHeatmap(obj, features = rownames(obj)))), 1L)
})

test_that("ezDoHeatmap renders above the cairo raster limit", {
  obj <- makeObject(32562L)
  expect_gt(pngDistinctBytes(renderToPng(
    ezDoHeatmap(obj, features = rownames(obj)))), 1L)
})

test_that("ezDoHeatmap keeps the fast raster path below the limit", {
  obj <- makeObject(5000L)
  geomOf <- function(p) class(p$layers[[1]]$geom)[1]
  expect_equal(geomOf(ezDoHeatmap(obj, features = rownames(obj))), "GeomRaster")
  expect_equal(geomOf(ezDoHeatmap(makeObject(32562L), features = rownames(obj))), "GeomTile")
})

test_that("an explicit raster argument is respected", {
  obj <- makeObject(1000L)
  expect_equal(class(ezDoHeatmap(obj, features = rownames(obj),
                                 raster = FALSE)$layers[[1]]$geom)[1], "GeomTile")
})
