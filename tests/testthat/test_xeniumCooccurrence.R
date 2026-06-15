test_that("computeCelltypeCooccurrence recovers clustered-vs-segregated structure", {
  skip_if_not_installed("dbscan")
  skip_if_not_installed("Matrix")
  # Two tight clusters far apart: A at origin, B ~100 um away.
  coords <- rbind(
    c(0, 0), c(1, 0), c(0, 1),       # A
    c(100, 100), c(101, 100), c(100, 101) # B
  )
  ct <- c("A", "A", "A", "B", "B", "B")
  res <- computeCelltypeCooccurrence(coords, ct, radius = 5, n_perm = 200,
                                     min_cells = 1, seed = 1)

  expect_equal(res$levels, c("A", "B"))
  expect_equal(dim(res$observed), c(2L, 2L))
  # No A-B edges across the 100 um gap
  expect_equal(res$observed["A", "B"], 0)
  # Same-type proximity enriched, cross-type depleted (both methods)
  expect_gt(res$log2enrich["A", "A"], 0)
  expect_lt(res$log2enrich["A", "B"], 0)
  expect_gt(res$zscore["A", "A"], 0)
  expect_lte(res$zscore["A", "B"], 0)
  # Matrices symmetric
  expect_equal(res$observed, t(res$observed))
  expect_equal(res$log2enrich, t(res$log2enrich))
})

test_that("Giotto log2 and squidpy z agree on the sign of enrichment", {
  skip_if_not_installed("dbscan")
  set.seed(3)
  # Checkerboard-ish: interleave so A and B are forced neighbours.
  g <- expand.grid(x = 1:8, y = 1:8)
  ct <- ifelse((g$x + g$y) %% 2 == 0, "A", "B")
  res <- computeCelltypeCooccurrence(as.matrix(g), ct, radius = 1.5,
                                     n_perm = 300, min_cells = 1, seed = 2)
  # On a checkerboard, A-B are the 4-connected neighbours (enriched),
  # A-A / B-B are not. Both statistics must agree on the sign for every pair.
  expect_true(all(sign(res$log2enrich) == sign(res$zscore) |
                    res$zscore == 0))
  expect_gt(res$log2enrich["A", "B"], 0)
  expect_lt(res$log2enrich["A", "A"], 0)
})

test_that("rare cell types are dropped via min_cells", {
  skip_if_not_installed("dbscan")
  coords <- rbind(c(0, 0), c(1, 0), c(0, 1), c(2, 2), c(50, 50))
  ct <- c("A", "A", "A", "B", "C") # B has 1, C has 1
  res <- computeCelltypeCooccurrence(coords, ct, radius = 5, n_perm = 50,
                                     min_cells = 2, seed = 1)
  expect_false("B" %in% res$levels)
  expect_false("C" %in% res$levels)
  expect_setequal(res$dropped_types, c("B", "C"))
})

test_that("subsampling respects max_cells", {
  skip_if_not_installed("dbscan")
  set.seed(1)
  coords <- matrix(runif(2000), ncol = 2)
  ct <- rep(c("A", "B"), 500)
  res <- computeCelltypeCooccurrence(coords, ct, radius = 0.1, n_perm = 20,
                                     min_cells = 1, max_cells = 200, seed = 1)
  expect_lte(res$n_cells, 200)
  expect_equal(res$n_total, 1000)
})

test_that("no edges and <2 types are handled without error", {
  skip_if_not_installed("dbscan")
  # radius too small -> no neighbours
  coords <- rbind(c(0, 0), c(10, 10), c(20, 20), c(30, 30))
  ct <- c("A", "A", "B", "B")
  res0 <- computeCelltypeCooccurrence(coords, ct, radius = 1e-6, n_perm = 10,
                                      min_cells = 1, seed = 1)
  expect_equal(res0$n_edges, 0)
  expect_null(res0$observed)
  # single type
  res1 <- computeCelltypeCooccurrence(coords, rep("A", 4), radius = 100,
                                      n_perm = 10, min_cells = 1)
  expect_null(res1$log2enrich)
  expect_equal(res1$levels, "A")
})

test_that("results are reproducible for a fixed seed", {
  skip_if_not_installed("dbscan")
  coords <- rbind(c(0, 0), c(1, 0), c(0, 1), c(100, 100), c(101, 100), c(100, 101))
  ct <- c("A", "A", "A", "B", "B", "B")
  r1 <- computeCelltypeCooccurrence(coords, ct, radius = 5, n_perm = 100, min_cells = 1, seed = 7)
  r2 <- computeCelltypeCooccurrence(coords, ct, radius = 5, n_perm = 100, min_cells = 1, seed = 7)
  expect_equal(r1$zscore, r2$zscore)
  expect_equal(r1$pval, r2$pval)
})
