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

test_that("getTargetTypeCounts takes promoter and downstream on the annotated strand in every strandMode", {
  ## one protein-coding transcript at 4001-5000 on a 10 kb chromosome; two
  ## reads fall in its true 2 kb promoter and one in its true 2 kb downstream
  ## window, so swapped windows would swap the counts
  make_gff <- function(geneStrand) {
    data.frame(
      seqid = "chr1",
      source = "protein_coding",
      type = "exon",
      start = 4001L,
      end = 5000L,
      score = ".",
      strand = geneStrand,
      phase = ".",
      attributes = 'gene_id "g1"; transcript_id "t1"; ID "t1-e1";',
      transcript_id = "t1",
      gene_id = "g1",
      stringsAsFactors = FALSE
    )
  }
  readStarts <- list(
    "+" = c(2500L, 3000L, 6000L),
    "-" = c(6000L, 6500L, 3000L)
  )
  readStrands <- list(
    sense = c("+" = "+", "-" = "-"),
    antisense = c("+" = "-", "-" = "+"),
    both = c("+" = "+", "-" = "-")
  )
  for (strandMode in names(readStrands)) {
    for (geneStrand in c("+", "-")) {
      reads <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(start = readStarts[[geneStrand]], width = 50L),
        strand = readStrands[[strandMode]][[geneStrand]],
        seqlengths = c(chr1 = 10000L)
      )
      result <- getTargetTypeCounts(
        list(strandMode = strandMode),
        make_gff(geneStrand),
        reads
      )
      label <- paste0("strandMode=", strandMode, ", gene strand ", geneStrand)
      expect_equal(result["mRNA Promoter 2kb", "count"], 2, label = label)
      expect_equal(result["mRNA Downstream 2kb", "count"], 1, label = label)
      expect_equal(result["mRNA Promoter 2kb", "width"], 2000, label = label)
      expect_equal(result["mRNA Downstream 2kb", "width"], 2000, label = label)
    }
  }
})
