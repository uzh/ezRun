test_that("getTargetTypeCounts uses the downstream windows for the width of mRNA Downstream 2kb", {
  ## one protein-coding transcript at 501-1500 on a 10 kb chromosome: the 2 kb
  ## promoter is trimmed at the chromosome start to 500 bases, the 2 kb
  ## downstream window (1501-3500) is not trimmed
  gff <- data.frame(
    seqid = "chr1",
    source = "protein_coding",
    type = "exon",
    start = 501L,
    end = 1500L,
    score = ".",
    strand = "+",
    phase = ".",
    attributes = 'gene_id "g1"; transcript_id "t1"; ID "t1-e1";',
    transcript_id = "t1",
    gene_id = "g1",
    stringsAsFactors = FALSE
  )
  reads <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(600L, 2000L), width = 50L),
    strand = "+",
    seqlengths = c(chr1 = 10000L)
  )
  param <- list(strandMode = "sense")
  result <- getTargetTypeCounts(param, gff, reads)
  expect_equal(result["mRNA Promoter 2kb", "width"], 500)
  expect_equal(result["mRNA Downstream 2kb", "width"], 2000)
  expect_equal(result["mRNA Downstream 2kb", "count"], 1)
})
