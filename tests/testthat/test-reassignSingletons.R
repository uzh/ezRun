## reassignSingletons / findClustersFast replace Seurat's GroupSingletons, whose
## singletons x clusters loop ran 32 h on p26168 TMA4 (14,007 singletons).

skip_if_not_installed("Seurat")
suppressPackageStartupMessages(library(Seurat))

## Clusters of unequal size with dense within-cluster weights, plus s
## singletons each tied mostly to one planted cluster. Every singleton also has
## weak ties to all of the largest cluster, so "highest SUM" picks the large
## cluster while Seurat's "highest MEAN" picks the planted one.
.plantedSnn <- function(sizes = c(20, 40, 60, 80, 160), s = 30, seed = 1) {
  set.seed(seed)
  k <- length(sizes)
  m <- sum(sizes) / k
  n <- sum(sizes) + s
  cells <- paste0("c", seq_len(n))
  clus <- rep(seq_len(k) - 1L, times = sizes)
  truth <- sample(seq_len(k) - 1L, s, replace = TRUE)
  w <- matrix(0, n, n, dimnames = list(cells, cells))
  for (j in seq_len(k)) {
    idx <- which(clus == j - 1L)
    w[idx, idx] <- runif(length(idx)^2, 0.3, 1)
  }
  for (i in seq_len(s)) {
    row <- k * m + i
    w[row, ] <- runif(n, 0, 0.02)
    w[row, which(clus == k - 1L)] <- 0.1
    home <- which(clus == truth[i])
    w[row, home] <- runif(length(home), 0.3, 0.6)
  }
  w <- pmax(w, t(w))
  diag(w) <- 1
  list(snn = Matrix::Matrix(w, sparse = TRUE), clus = clus, truth = truth,
       cells = cells, k = k, m = m, s = s)
}

test_that("reassignSingletons matches Seurat's GroupSingletons", {
  p <- .plantedSnn()
  ## Seurat's input: each singleton is its own community
  seuratIds <- setNames(c(as.character(p$clus), paste0("s", seq_len(p$s))), p$cells)
  ref <- Seurat:::GroupSingletons(seuratIds, p$snn, group.singletons = TRUE,
                                  verbose = FALSE)
  ours <- setNames(c(as.character(p$clus), rep("singleton", p$s)), p$cells)
  ours <- reassignSingletons(ours, p$snn)
  expect_identical(unname(as.character(ours)), unname(as.character(ref)))
  expect_identical(unname(ours[-seq_len(p$k * p$m)]), as.character(p$truth))
})

test_that("zero-connectivity singletons follow Seurat's tie-break", {
  ## TMA4 shape: singletons linked only to each other, never to a cluster.
  ## 9 clusters: set.seed(1); sample(9, 1) is 9, so "first cluster" would fail.
  p <- .plantedSnn(sizes = rep(12, 9), s = 0)
  n0 <- length(p$cells)
  nS <- 25
  cells <- c(p$cells, paste0("z", seq_len(nS)))
  w <- matrix(0, length(cells), length(cells), dimnames = list(cells, cells))
  w[seq_len(n0), seq_len(n0)] <- as.matrix(p$snn)
  zi <- n0 + seq_len(nS)
  w[zi, zi] <- 0.5
  diag(w) <- 1
  snn <- Matrix::Matrix(w, sparse = TRUE)
  seuratIds <- setNames(c(as.character(p$clus), paste0("s", seq_len(nS))), cells)
  ref <- Seurat:::GroupSingletons(seuratIds, snn, group.singletons = TRUE,
                                  verbose = FALSE)
  ours <- reassignSingletons(
    setNames(c(as.character(p$clus), rep("singleton", nS)), cells), snn
  )
  expect_identical(unname(as.character(ours)), unname(as.character(ref)))
})

test_that("reassignSingletons is fast on a TMA4-shaped graph", {
  set.seed(2)
  n <- 40000; nSingle <- 5000; k <- 40
  nnz <- n * 30
  snn <- Matrix::sparseMatrix(i = sample(n, nnz, TRUE), j = sample(n, nnz, TRUE),
                              x = runif(nnz), dims = c(n, n))
  snn <- snn + Matrix::t(snn)
  cells <- paste0("c", seq_len(n))
  dimnames(snn) <- list(cells, cells)
  ids <- setNames(as.character(sample(0:(k - 1), n, TRUE)), cells)
  ids[sample(n, nSingle)] <- "singleton"
  elapsed <- system.time(res <- reassignSingletons(ids, snn))[["elapsed"]]
  expect_false(any(res == "singleton"))
  expect_lt(elapsed, 5)
})

test_that("findClustersFast leaves no singleton label anywhere", {
  set.seed(3)
  m <- matrix(rpois(100 * 400, 2), 100, 400,
              dimnames = list(paste0("g", 1:100), paste0("c", 1:400)))
  o <- CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE))
  o <- NormalizeData(o, verbose = FALSE)
  VariableFeatures(o) <- rownames(o)
  o <- ScaleData(o, verbose = FALSE)
  o <- RunPCA(o, npcs = 20, verbose = FALSE)
  o <- FindNeighbors(o, dims = 1:20, verbose = FALSE)
  ## pure noise at a high resolution: Louvain leaves singletons (positive control)
  raw <- FindClusters(o, resolution = 8, group.singletons = FALSE, verbose = FALSE)
  expect_true("singleton" %in% raw$seurat_clusters)
  ## the app runs under plan("multicore"), where Seurat 5.5.1 ignores
  ## group.singletons = FALSE: our reassignment must still be the one that runs
  oplan <- future::plan("multicore", workers = 2)
  withr::defer(future::plan(oplan))
  logged <- character()
  futile.logger::flog.appender(function(line) logged <<- c(logged, line))
  withr::defer(futile.logger::flog.appender(futile.logger::appender.console()))
  fast <- findClustersFast(o, resolution = 8, verbose = FALSE)
  expect_true(any(grepl("reassigning", logged)))
  expect_identical(future::nbrOfWorkers(), 2L)
  expect_false("singleton" %in% fast$seurat_clusters)
  expect_false("singleton" %in% fast$RNA_snn_res.8)
  expect_false("singleton" %in% Idents(fast))
  expect_identical(levels(fast$seurat_clusters),
                   as.character(sort(as.integer(levels(fast$seurat_clusters)))))
})
