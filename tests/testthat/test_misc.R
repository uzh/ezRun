context("Tests miscellaneous functions: gage.r; igv.r; presentAnalysis.R; system.R")

param = ezParam(userParam = list('refBuild' = 'Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Version-2013-03-07'))
m1 = matrix(1:20,5)
rownames(m1) = letters[1:5]
colnames(m1) = letters[6:9]


test_that("Tests some functions from igv.r",{
  file = ezIgvTemplateFile()
  expect_is(file, "character")
  expect_true(ezIsAbsolutePath(file))
  genome = getIgvGenome(param)
  expect_is(genome, "character")
  expect_false(grepl("/", genome))
})

test_that("Tests ezPresentFlags() from presentAnalysis.R", {
  flags = ezPresentFlags(m1, param=list(useSigThresh=TRUE, sigThresh=10))
  expect_is(flags, "matrix")
  expect_is(flags[1, 1], "logical")
})

test_that("Tests ezSystem() and ezThreads() from system.R", {
  cmd = "echo Hello world!"
  res = ezSystem(cmd)
  expect_equal(res, 0)
  threads = ezThreads()
  expect_is(threads, "integer")
})
