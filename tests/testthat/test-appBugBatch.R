context("App bug batch helpers")

test_that("rmatsDartsBayesFdr is the running mean of 1 - posterior per event type", {
  post <- c(0.99, 0.5, 0.9, 0.95, 0.2, 0.8, NA)
  eventType <- c("SE", "SE", "SE", "RI", "RI", "SE", "SE")
  pass <- c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  fdr <- rmatsDartsBayesFdr(post, eventType, pass)
  ## SE passing: 0.99, 0.9, 0.5 -> 0.01, (0.01 + 0.1) / 2, (0.01 + 0.1 + 0.5) / 3
  expect_equal(fdr[c(1, 3, 2)], c(0.01, 0.055, 0.61 / 3))
  ## RI passing: 0.95, 0.2 -> 0.05, (0.05 + 0.8) / 2
  expect_equal(fdr[c(4, 5)], c(0.05, 0.425))
  ## below coverage or no posterior -> NA
  expect_true(is.na(fdr[6]))
  expect_true(is.na(fdr[7]))
  ## monotone in the posterior within an event type
  se <- which(eventType == "SE" & pass & !is.na(post))
  ord <- se[order(post[se], decreasing = TRUE)]
  expect_true(all(diff(fdr[ord]) >= 0))
  expect_true(all(fdr[!is.na(fdr)] >= 0 & fdr[!is.na(fdr)] <= 1))
})

test_that("isMainChromosome works for UCSC and Ensembl names", {
  keep <- c("chr1", "chr22", "chrX", "chrY", "1", "22", "X", "Y", "2L", "chr2R")
  drop <- c(
    "chrM", "chrMT", "M", "MT", "chrEBV", "chrUn_KI270302v1", "chr1_KI270706v1_random",
    "chr6_GL000250v2_alt", "KI270728.1", "GL000009.2", "chrUn", "Un", "chr4_random"
  )
  expect_true(all(isMainChromosome(keep)))
  expect_false(any(isMainChromosome(drop)))
})

test_that("bamDuplicatesRemoved reads @PG CL fields", {
  pg <- function(cl) paste0("@PG\tID:x\tPN:x\tVN:1\tCL:", cl)
  hdr <- c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:1000")
  expect_false(bamDuplicatesRemoved(headerLines = hdr))
  expect_false(bamDuplicatesRemoved(headerLines = c(hdr, pg("samtools markdup in.bam out.bam"))))
  expect_false(bamDuplicatesRemoved(headerLines = c(hdr, pg("samtools view -b -F 4 -o out.bam in.bam"))))
  expect_false(bamDuplicatesRemoved(headerLines = c(
    hdr, pg("MarkDuplicates --INPUT in.bam --REMOVE_DUPLICATES false")
  )))
  expect_true(bamDuplicatesRemoved(headerLines = c(hdr, pg("samtools view -b -@ 4 -F 1024 -o out.bam in.bam"))))
  expect_true(bamDuplicatesRemoved(headerLines = c(hdr, pg("samtools view -b -F 0x400 in.bam"))))
  expect_true(bamDuplicatesRemoved(headerLines = c(hdr, pg("samtools view -b -F 1028 in.bam"))))
  expect_true(bamDuplicatesRemoved(headerLines = c(hdr, pg("samtools view -b --exclude-flags DUP in.bam"))))
  expect_true(bamDuplicatesRemoved(headerLines = c(hdr, pg("samtools markdup -r -@ 4 in.bam out.bam"))))
  expect_true(bamDuplicatesRemoved(headerLines = c(hdr, pg("samtools rmdup in.bam out.bam"))))
  expect_true(bamDuplicatesRemoved(headerLines = c(
    hdr, pg("MarkDuplicates INPUT=in.bam OUTPUT=out.bam REMOVE_DUPLICATES=true")
  )))
  expect_true(bamDuplicatesRemoved(headerLines = c(
    hdr, pg("MarkDuplicates --INPUT in.bam --OUTPUT out.bam --REMOVE_DUPLICATES true")
  )))
})

test_that("import_macsPeak reads narrowPeak and broadPeak files", {
  narrowFile <- tempfile(fileext = ".narrowPeak")
  writeLines(c(
    "chr1\t99\t200\tpeak_1\t50\t.\t4.5\t10.2\t8.1\t40",
    "chr2\t499\t700\tpeak_2\t80\t.\t6.0\t12.0\t9.5\t100"
  ), narrowFile)
  narrow <- import_macsPeak(narrowFile)
  expect_equal(length(narrow), 2)
  expect_equal(GenomicRanges::start(narrow), c(100L, 500L))
  expect_equal(narrow$fold_enrichment, c(4.5, 6.0))
  expect_equal(narrow$peakScore, c(40L, 100L))

  broadFile <- tempfile(fileext = ".broadPeak")
  writeLines(c(
    "track name=broad",
    "chr1\t999\t5000\tbroad_1\t30\t.\t2.5\t5.1\t3.2",
    "chr3\t0\t800\tbroad_2\t20\t.\t1.8\t4.0\t2.2"
  ), broadFile)
  broad <- import_macsPeak(broadFile)
  expect_equal(length(broad), 2)
  expect_equal(broad$fold_enrichment, c(2.5, 1.8))
  expect_equal(broad$negLog10_qvalue, c(3.2, 2.2))
  expect_null(broad$peakScore)

  badFile <- tempfile(fileext = ".bed")
  writeLines("chr1\t1\t10\tx", badFile)
  expect_error(import_macsPeak(badFile), "expected a MACS")
})

test_that("DiffPeak HOMER annotation BED is 0-based like writeHomerBed", {
  peaks <- data.frame(
    Chr = c("chr1", "chr2"), Start = c(101, 2001), End = c(400, 2500),
    Strand = c("+", "-"), row.names = c("peak_1", "peak_2")
  )
  bed <- homerAnnotationBed(peaks)
  expect_equal(bed$Start, c(100L, 2000L))
  expect_equal(bed$End, c(400L, 2500L))
  expect_equal(bed$Names, c("peak_1", "peak_2"))
  homerFile <- tempfile(fileext = ".bed")
  writeHomerBed(data.frame(peaks, peakId = rownames(peaks)), homerFile)
  expect_equal(read.delim(homerFile, header = FALSE)$V2, bed$Start)
})

test_that("GatkJoinGenoTypes helpers: contig rename map and dbNSFP build check", {
  hdr <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=248956422>",
    "##contig=<ID=chrX,length=156040895>",
    "##contig=<ID=chrM,length=16569>",
    "##contig=<ID=KI270728.1,length=1872759>",
    "##INFO=<ID=chrom_note,Number=1,Type=String,Description=\"keep chr here\">"
  )
  map <- vcfChrRenameMap(hdr)
  expect_equal(map$old, c("chr1", "chrX", "chrM"))
  expect_equal(map$new, c("1", "X", "M"))
  expect_equal(nrow(vcfChrRenameMap(c("##contig=<ID=1,length=10>"))), 0)

  hg38 <- "Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_37-2021-05-04"
  hg19 <- "Homo_sapiens/UCSC/hg19/Annotation/Version-2014"
  expect_true(dbNsfpBuildCheck(hg38, "dbNSFP4.3_GRCh38.txt.gz")$use)
  expect_false(dbNsfpBuildCheck(hg19, "dbNSFP4.3_GRCh38.txt.gz")$use)
  expect_false(dbNsfpBuildCheck(hg38, "")$use)
  expect_true(dbNsfpBuildCheck(hg38, "dbNSFP4.3a")$use)
})

test_that("Mutect2 normal sample is read from the @RG SM field", {
  hdr <- c(
    "@HD\tVN:1.6",
    "@RG\tID:lane1\tSM:Normal_A\tLB:lib1\tPL:ILLUMINA",
    "@RG\tID:lane2\tSM:Normal_A\tLB:lib1\tPL:ILLUMINA"
  )
  expect_equal(bamReadGroupSamples(hdr), "Normal_A")
  expect_equal(
    bamReadGroupSamples(c(hdr, "@RG\tID:lane3\tPL:ILLUMINA\tSM:Other")),
    c("Normal_A", "Other")
  )
  expect_equal(length(bamReadGroupSamples("@HD\tVN:1.6")), 0)
})

test_that("Meme numeric parameters use exact names and defaults", {
  expect_equal(memeNumericParam(list(), "distanceToTSS", 2000), 2000)
  expect_equal(memeNumericParam(list(distanceToTSS = "500"), "distanceToTSS", 2000), 500)
  ## no partial matching of 'distance' etc.
  expect_equal(memeNumericParam(list(distanceToTSSx = "5"), "distanceToTSS", 2000), 2000)
  expect_error(memeNumericParam(list(minFold = "abc"), "minFold", 4), "must be a number")
})

test_that("VisiumQC slide type is HD, SD or mixed", {
  expect_equal(visiumQcType(c(TRUE, TRUE)), "HD")
  expect_equal(visiumQcType(c(FALSE, FALSE)), "SD")
  expect_equal(visiumQcType(c(TRUE, FALSE)), "mixed")
})

test_that("TargetedSeqQC BED filter keeps only main-chromosome regions", {
  bedFile <- tempfile(fileext = ".bed")
  writeLines(c(
    "track name=design",
    "chr1\t100\t200\tt1",
    "chrUn_KI270302v1\t1\t50\tt2",
    "chr2\t300\t400\tt3"
  ), bedFile)
  outFile <- tempfile(fileext = ".bed")
  .filterBedToContigs(bedFile, c("chr1", "chr2"), outFile)
  expect_equal(readLines(outFile), c("chr1\t100\t200\tt1", "chr2\t300\t400\tt3"))
  expect_error(.filterBedToContigs(bedFile, c("1", "2"), outFile), "main chromosomes")
})
