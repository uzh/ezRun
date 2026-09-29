## methods_facts() of the DNA mapping apps: BWA, Bowtie2 and DnaBamStats.
## Gated facts follow the typed parameters; anchors tie each literal value a fact states
## to the code of the function the fact names (deparse drops comments and fact text).

codeOfFn <- function(...) paste(unlist(lapply(list(...), deparse)), collapse = "\n")
has <- function(facts, pattern, fixed = TRUE) any(grepl(pattern, facts, fixed = fixed))
paramDefaults <- function() readLines(system.file("extdata/EZ_PARAM_DEFAULTS.txt", package = "ezRun", mustWork = TRUE))

## the BWA_dup run (p2220 BWA_2026-09-25--10-50-59), typed as ezParam would
bwaRun <- list(paired = TRUE, algorithm = "mem", cmdOptions = "", markDuplicates = TRUE, dupDistance = 2500,
               trimAdapter = TRUE, average_qual = 0, poly_x_min_len = 10, length_required = 18,
               nReads = -1, subsampleReads = 1)
## the Bowtie2_fixed run (p43268 o43311_Bowtie2_2026-09-24--17-03-05)
bowtieRun <- list(paired = TRUE, cmdOptions = "--no-unal -u 50000000", markDuplicates = TRUE, dupDistance = 2500,
                  generateBigWig = TRUE, secondRef = "", trimAdapter = TRUE, average_qual = 20,
                  poly_x_min_len = 10, length_required = 30, nReads = 5e7, subsampleReads = 1)
## the DnaBamStats_trout run (p33783 o43373_DNABamStats_2026-09-28--10-21-16)
dnaRun <- list(paired = TRUE, pixelDist = 2500, runQualimap = TRUE, runPicard = TRUE,
               keepProperPairsOnly = TRUE, refFeatureFile = "")

test_that("BWA facts: unconditional, gated and carrying the run's values", {
  app <- EzAppBWA$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 3)
  expect_true(has(f0, "added no other BWA option"))
  expect_true(has(f0, "no mapping-quality, proper-pair or duplicate filter"))
  expect_false(has(f0, "MarkDuplicates") || has(f0, "bwa mem was given") || has(f0, "bwa sampe") ||
                 has(f0, "subsampled"))

  f <- app$methods_facts(bwaRun)
  expect_lte(length(f), 12)
  expect_true(has(f, "bwa mem was given the read group"))
  expect_true(has(f, "PL ILLUMINA"))
  expect_true(has(f, "REMOVE_DUPLICATES=false"))
  expect_true(has(f, "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"))
  expect_true(has(f, "allIllumina-forTrimmomatic"))            # methodsFastpFacts, trimAdapter true
  expect_false(has(f, "bwa sampe") || has(f, "subsampled"))
  expect_true(has(app$methods_facts(modifyList(bwaRun, list(dupDistance = 100))), "OPTICAL_DUPLICATE_PIXEL_DISTANCE=100"))

  ## aln: sampe for paired, samse for single-end; no mem read-group fact
  aln <- app$methods_facts(modifyList(bwaRun, list(algorithm = "aln")))
  expect_true(has(aln, "bwa sampe"))
  expect_false(has(aln, "bwa mem was given"))
  expect_true(has(app$methods_facts(modifyList(bwaRun, list(algorithm = "aln", paired = FALSE))), "bwa samse"))
  ## markDuplicates off: no Picard, says not marked
  off <- app$methods_facts(modifyList(bwaRun, list(markDuplicates = FALSE)))
  expect_false(has(off, "MarkDuplicates"))
  expect_true(has(off, "Duplicates were not marked (markDuplicates false)"))
})

test_that("Bowtie2 facts: unconditional, gated and carrying the run's values", {
  app <- EzAppBowtie2$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 3)
  expect_true(has(f0, "set no alignment mode or scoring option itself"))
  expect_true(has(f0, "LB RGLB_<sample>, PL illumina, PU RGPU_<sample>"))
  expect_false(has(f0, "bigWig") || has(f0, "MarkDuplicates") || has(f0, "secondRef") || has(f0, "subsampled"))

  f <- app$methods_facts(bowtieRun)
  expect_lte(length(f), 12)
  bw <- grep("bigWig", f, value = TRUE)
  expect_length(bw, 1)
  expect_match(bw, "readGAlignmentPairs", fixed = TRUE)
  expect_match(bw, "no CPM or other normalisation", fixed = TRUE)
  expect_true(has(f, "subsampled to 50,000,000 reads"))
  expect_true(has(f, "seed 123"))
  expect_true(has(f, "Read1 and Read2 separately"))
  expect_true(has(f, "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"))
  expect_false(has(f, "secondRef"))

  expect_false(has(app$methods_facts(modifyList(bowtieRun, list(generateBigWig = FALSE))), "bigWig"))
  expect_match(grep("bigWig", app$methods_facts(modifyList(bowtieRun, list(paired = FALSE))), value = TRUE),
               "(readGAlignments)", fixed = TRUE)
  expect_false(has(app$methods_facts(modifyList(bowtieRun, list(nReads = -1))), "subsampled"))
  sub <- app$methods_facts(modifyList(bowtieRun, list(nReads = -1, subsampleReads = 4, paired = FALSE)))
  expect_true(has(sub, "Read Count / 4 reads"))
  expect_false(has(sub, "Read1 and Read2"))
  sr <- app$methods_facts(modifyList(bowtieRun, list(secondRef = "/srv/GT/x/extra.fa")))
  expect_true(has(sr, "bowtie2-build --seed 42"))
  expect_true(has(app$methods_facts(modifyList(bowtieRun, list(markDuplicates = FALSE))), "Duplicates were not marked"))
})

test_that("DnaBamStats facts: unconditional, gated and carrying the run's values", {
  app <- EzAppDnaBamStats$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 2)
  expect_true(has(f0, "counting the mapped alignment records per read name"))
  expect_true(has(f0, "mapping rate below 70%, average coverage below 10x, duplication rate above 50%"))
  expect_false(has(f0, "Qualimap") || has(f0, "Picard") || has(f0, "ATACseqQC") || has(f0, "proper pair"))

  f <- app$methods_facts(dnaRun)
  expect_gte(length(f), 4)
  expect_lte(length(f), 12)
  expect_true(has(f, "only paints chromosome limits"))
  expect_true(has(f, "not Qualimap's 'number of mapped reads'"))
  expect_true(has(f, "OPTICAL_DUPLICATE_PIXEL_DISTANCE equal to pixelDist (2500)"))
  expect_true(has(f, "READ_PAIR_OPTICAL_DUPLICATES / READ_PAIRS_EXAMINED"))
  expect_true(has(f, "ATACseqQC fragSizeDist"))
  expect_true(has(f, "only first mates in a proper pair"))

  g <- app$methods_facts(modifyList(dnaRun, list(runQualimap = FALSE, runPicard = FALSE, paired = FALSE)))
  expect_false(has(g, "Qualimap") || has(g, "Picard") || has(g, "MarkDuplicates") || has(g, "ATACseqQC") ||
                 has(g, "proper pair"))
  np <- app$methods_facts(modifyList(dnaRun, list(keepProperPairsOnly = FALSE)))
  expect_true(has(np, "only first mates were counted"))
  expect_false(has(np, "in a proper pair"))
  expect_true(has(app$methods_facts(modifyList(dnaRun, list(pixelDist = 100))), "equal to pixelDist (100)"))
})

test_that("mapping facts are anchored in the code they describe", {
  ## BWA
  bwa <- codeOfFn(ezMethodBWA)
  expect_match(bwa, "PL:ILLUMINA", fixed = TRUE)
  expect_match(bwa, "\"sampe\"", fixed = TRUE)
  expect_match(bwa, "\"samse\"", fixed = TRUE)
  expect_match(bwa, "operation = \"mark\"", fixed = TRUE)
  expect_match(bwa, "\"-R\",\\s*readGroupOpt,\\s*\"-t\"", perl = TRUE)
  expect_no_match(bwa, "-q |filteroutBam|removeDuplicates", perl = TRUE)
  dup <- codeOfFn(dupBam)
  expect_match(dup, "OPTICAL_DUPLICATE_PIXEL_DISTANCE=", fixed = TRUE)
  expect_match(dup, "if_else\\(operation ==\\s*\"mark\",\\s*\"false\",\\s*\"true\"\\)", perl = TRUE)
  expect_match(dup, "Rsamtools::indexBam(outBam)", fixed = TRUE)
  expect_match(dup, "sub(\".bam$\", \"_metrics.txt\", outBam)", fixed = TRUE)
  expect_true(any(grepl("^dupDistance\tnumeric\t2500\t", paramDefaults())))
  ## Bowtie2
  bt <- codeOfFn(ezMethodBowtie2)
  expect_match(bt, "--rg PL:illumina", fixed = TRUE)
  expect_match(bt, "--rg LB:RGLB_", fixed = TRUE)
  expect_match(bt, "--rg PU:RGPU_", fixed = TRUE)
  expect_match(bt, "method = \"Bioconductor\"", fixed = TRUE)
  expect_match(bt, "Coverage_", fixed = TRUE)
  expect_match(bt, "reads_per_chromosome * readSize/chrLengths", fixed = TRUE)
  bw <- codeOfFn(bam2bw)
  expect_match(bw, "readGAlignmentPairs(file)", fixed = TRUE)
  expect_match(bw, "export.bw(cov, destination)", fixed = TRUE)
  expect_match(codeOfFn(getBowtie2Reference), "--seed 42", fixed = TRUE)
  ## subsampling (ezMethodFastpTrim -> ezMethodSubsampleFastq -> subsampleFastqFile)
  expect_match(codeOfFn(subsampleFastqFile), "seed = 123L", fixed = TRUE)
  expect_match(codeOfFn(subsampleFastqFile), "FastqSampler(inFile", fixed = TRUE)
  trim <- codeOfFn(ezMethodFastpTrim)
  expect_match(trim, "min(param$nReads, totalReads)", fixed = TRUE)
  expect_match(trim, "1/param$subsampleReads * totalReads", fixed = TRUE)
  ## DnaBamStats
  expect_match(codeOfFn(getBamMultiMatchingFromNH), "tag = \"NH\"", fixed = TRUE)
  qn <- codeOfFn(getBamMultiMatchingFromQnames)
  expect_match(qn, "isProperPair = ezBamFlagKeepOnly(param$keepProperPairsOnly)", fixed = TRUE)
  expect_match(qn, "isFirstMateRead = TRUE", fixed = TRUE)
  expect_no_match(qn, "isSecondaryAlignment|isSupplementary", perl = TRUE)
  expect_match(codeOfFn(getBamMultiMatching), "nReads - sum(result)", fixed = TRUE)
  expect_true(any(grepl("^keepProperPairsOnly\tlogical\tTRUE\t", paramDefaults())))
  qm <- codeOfFn(get_dna_qualimap_stats)
  expect_match(qm, "\" -c -nt \"", fixed = TRUE)
  expect_no_match(qm, "-gff|refFeatureFile", perl = TRUE)
  expect_match(qm, "\"^.*number of reads = \"", fixed = TRUE)
  expect_match(qm, "qualimapReads/2", fixed = TRUE)
  expect_match(qm, "100 * qualimapReads/nReads", fixed = TRUE)
  expect_match(codeOfFn(get_dna_qualimap_multi_sample_summary), "qualimap multi-bamqc", fixed = TRUE)
  expect_match(codeOfFn(dna_bamstats_can_reuse_dup_metrics), "existingPixelDist == as.numeric(pixelDist)", fixed = TRUE)
  pic <- codeOfFn(get_dna_picard_dup_stats)
  expect_match(pic, "REMOVE_DUPLICATES=false", fixed = TRUE)
  expect_match(pic, "dna_bamstats_can_reuse_dup_metrics(dupMetricsFile, param$pixelDist)", fixed = TRUE)
  expect_match(codeOfFn(parse_dna_picard_dup_metrics_file), "100 * optDuplicates/readPairs", fixed = TRUE)
  expect_match(codeOfFn(get_dna_paired_end_plots), "fragSizeDist(bamFile, sampleName)", fixed = TRUE)
  lc <- codeOfFn(ez_dna_lib_complexity)
  expect_match(lc, "set.seed(42)", fixed = TRUE)
  expect_match(lc, "times = 100", fixed = TRUE)
  expect_match(lc, "seq(5, 20, by = 5)", fixed = TRUE)
  expect_match(codeOfFn(get_dna_bamstats_qc_thresholds), "Threshold = c(70, 10, 50, 20)", fixed = TRUE)
})
