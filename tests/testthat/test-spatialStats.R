## Synthetic Visium-like grid: even rows use even columns, odd rows odd columns
make_hex_grid <- function(nRow = 20, nColHalf = 20) {
  grid <- expand.grid(row = seq_len(nRow) - 1, k = seq_len(nColHalf) - 1)
  grid$col <- 2 * grid$k + grid$row %% 2
  grid$name <- paste0("s", seq_len(nrow(grid)))
  grid
}

test_that("hex graph gives 6 neighbours inside and 18 with two rings", {
  g <- make_hex_grid()
  adj1 <- spatial_hex_graph(g$row, g$col, g$name, rings = 1)
  adj2 <- spatial_hex_graph(g$row, g$col, g$name, rings = 2)
  inner <- g$row %in% 3:16 & g$k %in% 3:16
  expect_true(all(Matrix::rowSums(adj1)[inner] == 6))
  expect_true(all(Matrix::rowSums(adj2)[inner] == 18))
  expect_true(Matrix::isSymmetric(adj1))
  expect_equal(sum(Matrix::diag(adj1)), 0)
})

test_that("Moran's I separates a spatial gradient from noise", {
  g <- make_hex_grid()
  adj <- spatial_hex_graph(g$row, g$col, g$name)
  set.seed(1)
  mat <- rbind(
    gradient = g$col + rnorm(nrow(g), sd = 1),
    noise = rnorm(nrow(g)),
    constant = rep(1, nrow(g))
  )
  colnames(mat) <- g$name
  res <- spatial_morans_i(mat, adj)
  expect_gt(res["gradient", "moransI"], 0.8)
  expect_lt(abs(res["noise", "moransI"]), 0.1)
  expect_lt(res["gradient", "fdr"], 1e-10)
  expect_true(is.na(res["constant", "moransI"]))
  expect_equal(rownames(res)[1], "gradient")
})

test_that("adjacency enrichment finds planted touching domains", {
  g <- make_hex_grid()
  adj <- spatial_hex_graph(g$row, g$col, g$name)
  ## three vertical stripes: A | B | C; A and C never touch
  lab <- cut(g$col, breaks = c(-1, 12, 26, 100), labels = c("A", "B", "C"))
  res <- spatial_adjacency_enrichment(lab, adj, nPerm = 100)
  expect_gt(res$zScore["A", "A"], 5)
  expect_equal(res$observed["A", "C"], 0)
  expect_lt(res$zScore["A", "C"], 0)
  expect_gt(res$zScore["A", "B"], res$zScore["A", "C"])
})

test_that("smoothness, ARI and colocalisation behave", {
  g <- make_hex_grid()
  adj <- spatial_hex_graph(g$row, g$col, g$name)
  stripes <- ifelse(g$col < 20, "left", "right")
  set.seed(2)
  random <- sample(stripes)
  expect_gt(
    mean(spatial_label_smoothness(stripes, adj)),
    mean(spatial_label_smoothness(random, adj)) + 0.3
  )
  expect_equal(adjusted_rand_index(stripes, stripes), 1)
  expect_lt(abs(adjusted_rand_index(stripes, random)), 0.1)
  scores <- cbind(
    left = as.numeric(g$col < 20),
    right = as.numeric(g$col >= 20),
    alsoLeft = as.numeric(g$col < 20) * runif(nrow(g))
  )
  co <- spatial_colocalisation(scores, adj)
  expect_lt(co$lagged["left", "right"], -0.5)
  expect_gt(co$lagged["left", "alsoLeft"], 0.3)
})

test_that("domain labels use margin and fall back to mixed", {
  m <- rbind(
    d1 = c(Cortex = 1.5, Striatum = 0.2, Thalamus = -0.5),
    d2 = c(Cortex = 0.6, Striatum = 0.5, Thalamus = 0.1),
    d3 = c(Cortex = -0.2, Striatum = -0.3, Thalamus = -0.1)
  )
  res <- assign_domain_labels(m, minMargin = 0.25)
  expect_equal(res$label[1], "Cortex")
  expect_match(res$label[2], "^mixed \\(Cortex/Striatum\\)")
  expect_match(res$label[3], "^mixed")
})
