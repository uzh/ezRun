context("Test annotation and gtf: annotation.r; gff.r; go-analysis.R; ngsReferenceFiles.r")

library(rtracklayer)
param <- list()
param[['refBuild']] <- 'Saccharomyces_cerevisiae/Ensembl/R64/Annotation/Release_98-2019-12-03'
param[['refFeatureFile']] <- 'genes.gtf'
param <- ezParam(param)
gtfGR <- import(param$ezRef@refFeatureFile)
gtfDF <- ezReadGff(param$ezRef@refFeatureFile)
gtf <- ezLoadFeatures(param, param$ezRef@refFeatureFile)
txAnno <- ezFeatureAnnotation(param, dataFeatureType="transcript")

test_that("Tests functions in annotation.r", {
  geneAnno <- ezFeatureAnnotation(param, dataFeatureType="gene")
  expect_is(txAnno, "data.frame")
  expect_setequal(rownames(geneAnno), gtfGR$gene_id)
  expect_setequal(rownames(txAnno), na.omit(gtfGR$transcript_id))
  
  geneMapping <- getGeneMapping(param, txAnno)
  expect_setequal(geneMapping, gtfGR$gene_id)
  expect_setequal(names(geneMapping), na.omit(gtfGR$transcript_id))
  
  hasMapping <- hasGeneMapping(param, txAnno)
  expect_true(hasMapping)
})

test_that("Tests basic functions in gff.r", {
  expect_is(gtf, "data.frame")
  
  transcripts <- ezGffAttributeField(gtfDF$attributes, field="transcript_id", 
                                     attrsep="; ", valuesep=" ")
  expect_is(transcripts, "character")
  expect_identical(length(gtfGR), length(transcripts))
  
  newAttribute <- ezBuildAttributeField(mcols(gtfGR)[ , c("gene_id", "gene_name", "transcript_id")],
                                       "gtf")
  expect_is(newAttribute, "character")
  expect_identical(length(gtfGR), length(newAttribute))
})

test_that("Tests more functions from gff.r", {
  exonNumber = getExonNumber(gtf)
  expect_is(exonNumber, "integer")
  expect_identical(length(exonNumber), nrow(gtf))
  
  transcriptDB = ezTranscriptDbFromRef(param$ezRef)
  expect_is(transcriptDB, "TxDb")
  
  tr2Gene = getTranscript2Gene(gtf)
  expect_is(tr2Gene, "character")
  
  counts = countIsoformsPerGene(tr2Gene)
  expect_is(counts, "table")
  
  transcriptSequences = getTranscriptSequences(param)
  expect_is(transcriptSequences, "DNAStringSet")
  
  gcAndWidth = getTranscriptGcAndWidth(param)
  expect_is(gcAndWidth$gc, "numeric")
  expect_is(gcAndWidth$featWidth, "integer")
  
  ensemblTypes = getEnsemblTypes(gtf)
  expect_is(ensemblTypes, "character")
})

test_that("Tests annotation functions related to GO (from go-analysis.R)", {
  do <- doGo(param, txAnno)
  expect_is(do, "logical")
  
  has <- hasGoAnnotation(txAnno)
  expect_is(has, "logical")
  expect_null(getGOparents("noGoTerm"))
  
  parents <- getGOparents("GO:0034767")
  expect_is(parents, "character")
  expect_true(all(grepl("GO", names(parents))))
  
  added <- addGoParents(c("GO:0034767", "GO:0034768"), "BP")
  expect_is(added, "list")
  expect_true(all(grepl("GO", unlist(added))))
  
  goList <- goStringsToList(c("GO:0034762; GO:0034763", "GO:0034764; GO:0034765"))
  expect_is(goList, "list")
  expect_equal(length(unlist(goList)), 4)
})

test_that("Tests cleanGenomeFiles() from ngsReferenceFiles.r", {
  cleanedGenome <- cleanGenomeFiles(param$ezRef@refFastaFile,
                                    param$ezRef@refFeatureFile)
  expect_is(cleanedGenome$genomeSeq, "DNAStringSet")
  expect_is(cleanedGenome$gtf, "GRanges")
})
