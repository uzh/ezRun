## Cheap helpers of R/spatialCommunication.R that the spatial apps rely on to
## keep the cell-cell communication and niche steps affordable.

test_that("spatial_coords_um scales pixels to microns", {
  co <- spatial_coords_um(c(0, 100), c(10, 20), 0.5, c("a", "b"))
  expect_equal(co$x, c(0, 50))
  expect_equal(co$y, c(5, 10))
  expect_identical(rownames(co), c("a", "b"))
})

test_that("aggregate_pseudospots sums counts and averages fitted weights", {
  ## four 16 um bins in one 48 um square, one bin in the next square
  co <- data.frame(x = c(1, 17, 33, 1, 60), y = c(1, 1, 1, 17, 1), row.names = paste0("b", 1:5))
  counts <- Matrix::Matrix(rbind(g1 = c(1, 2, 3, 4, 5), g2 = c(0, 1, 0, 1, 0)), sparse = TRUE)
  colnames(counts) <- rownames(co)
  w <- rbind(c(1, 0), c(0, 1), c(NA, NA), c(0.5, 0.5), c(0.2, 0.8))
  dimnames(w) <- list(rownames(co), c("A", "B"))
  agg <- aggregate_pseudospots(counts, co, sizeUm = 48, weights = w)
  expect_equal(ncol(agg$counts), 2)
  first <- agg$membership[["b1"]]
  expect_equal(unname(agg$counts[, first]), c(10, 2))
  expect_equal(unname(agg$nMembers[first]), 4)
  ## b3 has no weights and must not dilute the average
  expect_equal(unname(agg$weights[first, ]), c(0.5, 0.5))
  expect_equal(unname(agg$weights[agg$membership[["b5"]], ]), c(0.2, 0.8))
  expect_equal(unname(agg$coords[first, "x"]), mean(c(1, 17, 33, 1)))
})

test_that("aggregate_pseudospots drops sparse pseudo-spots", {
  co <- data.frame(x = c(1, 2, 100), y = c(1, 2, 1), row.names = c("a", "b", "c"))
  agg <- aggregate_pseudospots(coordsUm = co, sizeUm = 50, minMembers = 2)
  expect_equal(length(agg$nMembers), 1)
  expect_true(is.na(agg$membership[["c"]]))
})

test_that("stratified_subsample caps every group and keeps small ones", {
  labels <- setNames(rep(c("big", "small"), c(100, 5)), paste0("s", 1:105))
  keep <- stratified_subsample(labels, maxPerGroup = 20)
  expect_equal(as.vector(table(labels[keep])), c(20, 5))
  expect_identical(keep, stratified_subsample(labels, maxPerGroup = 20))
  expect_false(is.unsorted(match(keep, names(labels))))
})

test_that("with_future_plan restores the previous plan", {
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  inside <- with_future_plan(2, future::nbrOfWorkers())
  expect_equal(inside, 2)
  expect_equal(future::nbrOfWorkers(), 1)
  expect_error(with_future_plan(2, stop("boom")), "boom")
  expect_equal(future::nbrOfWorkers(), 1)
})

test_that("select_weight_types keeps abundant, dominating types", {
  set.seed(1)
  w <- cbind(A = runif(50, 0.5, 1), B = runif(50, 0, 0.5), C = 0.001)
  w[1:20, "B"] <- 0.99
  w <- w / rowSums(w)
  expect_identical(select_weight_types(w, minSpots = 10), c("A", "B"))
})

test_that("dominant_type skips spots without weights", {
  w <- rbind(s1 = c(A = 0.7, B = 0.3), s2 = c(NA, NA), s3 = c(0.1, 0.9))
  d <- dominant_type(w)
  expect_identical(names(d), c("s1", "s3"))
  expect_identical(as.character(d), c("A", "B"))
  expect_identical(levels(d), c("A", "B"))
})
