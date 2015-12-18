context("Test annotation and gtf: annotation.r; gff.r; go-analysis.R")

param = ezParam()
param$ezRef@refFeatureFile = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
param$ezRef@refAnnotationFile = ""
fp = "/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/WholeGenomeFasta/genome.fa"
param$ezRef@refFastaFile = fp
seqAnno = writeAnnotationFromGtf(param)

test_that("Tests functions in annotation.r", {
  featureAnno = ezFeatureAnnotation(param, rownames(seqAnno), "gene")
  expect_is(seqAnno, "data.frame")
  expect_identical(names(seqAnno), names(featureAnno))
  geneMapping = getGeneMapping(param, seqAnno)
  expect_is(geneMapping, "array")
  hasMapping = hasGeneMapping(param, seqAnno)
  expect_true(hasMapping)
})

test_that("", {
  
})

test_that("Tests annotation functions related to GO (from go-analysis.R)", {
  do = doGo(param, seqAnno)
  expect_is(do, "logical")
  has = hasGoAnnotation(seqAnno)
  expect_is(has, "logical")
  expect_null(getGOparents("noGoTerm"))
  parents = getGOparents("GO:0034767")
  expect_is(parents, "character")
  expect_true(all(grepl("GO", names(parents))))
  added = addGoParents(c("GO:0034767", "GO:0034768"), "BP")
  expect_is(added, "list")
  expect_true(all(grepl("GO", unlist(added))))
  goList = goStringsToList(c("GO:0034762; GO:0034763", "GO:0034764; GO:0034765"))
  expect_is(goList, "list")
  expect_equal(length(unlist(goList)), 4)
})
