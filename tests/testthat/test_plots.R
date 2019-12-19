context("Tests the functions from plots.R")

file = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
ds = EzDataset$new(file=file, dataRoot=system.file(".", package="ezRun", mustWork = TRUE))
cond = ezConditionsFromDataset(ds$meta)
types = data.frame(matrix(rep(1:10, each=10), 10))

test_that("Tests getSampleColors(), getSamplePch(), getSampleLty() and getBlueRedScale()", {
  colors = getSampleColors(cond)
  pch = getSamplePch(cond)
  lty = getSampleLty(cond)
  scale = getBlueRedScale()
  expect_is(colors, "character")
  expect_true(all(grepl("#", colors)))
  expect_is(pch, "integer")
  expect_is(lty, "integer")
  expect_is(scale, "character")
  expect_true(all(grepl("#", scale)))
})

test_that("Tests ezColorLegend() and ezLegend()", {
  ret = ezColorLegend()
  expect_is(ret, "numeric")
  ret2 = ezLegend(1:3)
  expect_is(ret2, "list")
})

test_that("Tests ezVolcano()", {
  ret = ezVolcano(log2Ratio=1:100, pValue=rep(10^(-4:5), each=10),
                  isPresent=1:50, types=types, colors=rainbow(ncol(types)))#, pch=16, legendPos="bottomleft")
  expect_is(ret, "list")
  expect_is(ret$x, "numeric")
})

test_that("Tests ezSmoothScatter()", {
  ret = ezSmoothScatter(y=data.frame(a=1:10, b=21:30, c=41:50))
  expect_null(ret)
})

test_that("Tests ezScatter() and ezXYScatter()", {
  x = runif(n=1000)
  y = runif(n=1000)
  isPresent = x > 0.2 & y > 0.2
  ret = ezScatter(y=data.frame(a=1:10, b=21:30, c=41:50))
  ret2 = ezXYScatter(x, y, isPresent=isPresent)
  expect_null(ret)
  expect_null(ret2)
})

test_that("Tests ezAllPairScatter()", {
  ret = ezAllPairScatter(x=matrix(1:10,5))
  expect_null(ret)
})

test_that("Tests ezCorrelationPlot()", {
  ret = ezCorrelationPlot(z=matrix(1:100,10))
  expect_is(ret, "numeric")
})

test_that("Tests intHist()", {
  ret = intHist(1:10)
  expect_is(ret, "histogram")
})

test_that("Tests ezHeatmap()", {
  x = matrix(1:100, 10)
  colnames(x) = letters[1:10]
  ret = ezHeatmap(x)
  expect_is(ret, "list")
})

test_that("Tests colorClusterLabel()", {
  hc = hclust(dist(USArrests), "ave")
  dend = as.dendrogram(hc)
  labels = colorClusterLabels(dend, rainbow(6))
  expect_is(labels, "dendrogram")
})
