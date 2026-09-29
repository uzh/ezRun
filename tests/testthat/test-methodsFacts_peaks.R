## methods_facts() of the peak apps: NfCoreAtacSeq, MACS3, PeakCombiner.
## Gated facts follow the typed parameters (and the input dataset); anchors tie each
## literal value a fact states to the code of the function the fact names.

codeOfFn <- function(...) paste(unlist(lapply(list(...), function(f) deparse(f, width.cutoff = 500L))), collapse = "\n")
has <- function(facts, pattern, fixed = TRUE) any(grepl(pattern, facts, fixed = fixed))
withInput <- function(p, df) { attr(p, "input") <- df; p }

test_that("NfCoreAtacSeq facts: unconditional, gated and carrying the run's values", {
  app <- EzAppNfCoreAtacSeq$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 4)
  expect_true(has(f0, "-profile apptainer"))
  expect_true(has(f0, "80% of the summed sequence lengths"))
  expect_true(has(f0, "software_versions.yml"))
  expect_true(has(f0, "ezRun default 2.1.2"))
  expect_false(has(f0, "qcMode") || has(f0, "narrow_peak") || has(f0, "keepBams"))

  run <- list(pipelineVersion = "2.1.2", peakStyle = "broad", varStabilizationMethod = "vst",
              qcMode = TRUE, qcReadsPerSample = 1e7, keepBams = FALSE, grouping = "Condition",
              cmdOptions = "")
  f <- app$methods_facts(run)
  expect_true(has(f, "-r 2.1.2"))
  expect_true(has(f, "first 10000000 reads (pairs)"))
  expect_true(has(f, "broadPeak mode (no --narrow_peak)"))
  expect_true(has(f, "BAM and BAI files were deleted"))
  expect_true(has(f, "vst setting"))
  expect_false(has(f, "cmdOptions"))
  ## off / other values
  g <- app$methods_facts(modifyList(run, list(qcMode = FALSE, keepBams = TRUE, peakStyle = "narrow",
                                              varStabilizationMethod = "rlogTransf", cmdOptions = "--skip_fastqc",
                                              pipelineVersion = "2.1.1")))
  expect_false(has(g, "first 10000000 reads") || has(g, "BAM and BAI files were deleted"))
  expect_true(has(g, "--narrow_peak"))
  expect_true(has(g, "--deseq2_vst false"))
  expect_true(has(g, "-r 2.1.1"))
  expect_true(has(g, "--skip_fastqc"))

  ## grouping from the input dataset (the round-8 AtacSeq_ok run had Condition NA everywhere)
  empty <- withInput(run, data.frame(Name = c("N704", "N705"), `Condition [Factor]` = NA, check.names = FALSE))
  expect_true(has(app$methods_facts(empty), "Condition was empty for every sample"))
  part <- withInput(run, data.frame(Name = c("a", "b", "c"), `Condition [Factor]` = c("x", "x", NA), check.names = FALSE))
  expect_true(has(app$methods_facts(part), "sample(s) c"))
  full <- withInput(run, data.frame(Name = c("a", "b"), `Condition [Factor]` = c("x", "x"), check.names = FALSE))
  expect_true(has(app$methods_facts(full), "numbered as replicates"))
  expect_false(has(app$methods_facts(full), "was empty for every sample"))
  nocol <- withInput(run, data.frame(Name = c("a", "b"), check.names = FALSE))
  expect_true(has(app$methods_facts(nocol), "Condition was empty for every sample"))
})

test_that("NfCoreAtacSeq anchors: the values the facts state are in the code", {
  f <- EzAppNfCoreAtacSeq$new()$methods_facts(list(qcMode = TRUE, qcReadsPerSample = 1e7, keepBams = FALSE, peakStyle = "broad"))
  for (v in c("-profile apptainer", "80% of", "--bwa_index", "genes.bed", "first 10000000 reads", "genome, trimgalore, fastqc and igv",
              "replicate", "consensus")) expect_true(has(f, v), label = v)
  cmd <- codeOfFn(ezRun:::buildNfCoreAtacCmd)
  for (v in c("fullGenomeSize * 0.8", "-profile apptainer", "\"-r\"", "param$pipelineVersion",
              "--narrow_peak", "--deseq2_vst false", "--bwa_index", "Sequence/BWAIndex", "genes.bed",
              "param$cmdOptions", "--macs_gsize")) {
    expect_match(cmd, v, fixed = TRUE)
  }
  expect_identical(EzAppNfCoreAtacSeq$new()$appDefaults["pipelineVersion", "DefaultValue"], "2.1.2")
  expect_match(codeOfFn(ezRun:::subsampleAtacFastqs), "head -n", fixed = TRUE)
  expect_match(codeOfFn(ezRun:::subsampleAtacFastqs), "4 * nReads", fixed = TRUE)
  expect_match(codeOfFn(ezRun:::subsampleAtacFastqs), "readCounts[sids] > nReads", fixed = TRUE)
  main <- codeOfFn(ezRun:::ezMethodNfCoreAtacSeq)
  expect_match(main, "c(\"genome\", \"trimgalore\", \"fastqc\", \"igv\")", fixed = TRUE)
  expect_match(main, "_peak/consensus/", fixed = TRUE)
  expect_match(codeOfFn(ezRun:::getAtacSampleSheet), "ezReplicateNumber(groups)", fixed = TRUE)
  expect_match(codeOfFn(ezRun:::getAtacSampleSheet), "groups[isMissing] <- sampleNames[isMissing]", fixed = TRUE)
  expect_match(codeOfFn(ezRun:::cleanupAtacOutFolder), "*.bam(.bai)?$", fixed = TRUE)
  expect_match(codeOfFn(ezRun:::writeNextflowLimits), "resourceLimits", fixed = TRUE)
})

test_that("MACS3 facts: unconditional, gated and carrying the run's values", {
  app <- EzAppMacs3$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 3)
  expect_true(has(f0, "--keep-dup all"))
  expect_true(has(f0, "bedtools getfasta"))
  expect_false(has(f0, "ATAC-seq mode") || has(f0, "ChIPseeker") || has(f0, "alignmentSieve"))

  ## the round-8 run o43054_MACS3_2026-09-17
  atac <- list(mode = "ATAC-seq", paired = TRUE, useControl = FALSE, shiftATAC = TRUE, annotatePeaks = TRUE,
               cmdOptions = "--nomodel --bw 200 --extsize 147", genomeSize = 0, qValue = 0.05,
               removeDuplicates = TRUE)
  f <- app$methods_facts(atac)
  expect_true(has(f, "-q 0.05"))
  expect_true(has(f, "-f BAMPE"))
  expect_true(has(f, "--extsize 147 given in cmdOptions was replaced by --extsize 200"))
  expect_true(has(f, "MAPQ of at least 10"))
  expect_true(has(f, "alignmentSieve --ATACshift"))
  expect_true(has(f, "samtools view -F 1024"))
  expect_true(has(f, "80% of the summed sequence lengths"))
  expect_true(has(f, "bamCoverage"))
  expect_true(has(f, "ChIPseeker"))
  expect_false(has(f, "bdgcmp"))
  ## switched off / other values
  g <- app$methods_facts(modifyList(atac, list(shiftATAC = FALSE, annotatePeaks = FALSE, removeDuplicates = FALSE,
                                               genomeSize = 2.7e9, cmdOptions = "--nomodel --keep-dup 1")))
  expect_false(has(g, "alignmentSieve") || has(g, "ChIPseeker") || has(g, "samtools view -F 1024"))
  expect_false(has(g, "--keep-dup all") || has(g, "80% of the summed"))
  expect_true(has(g, "genomeSize 2.7e+09 was not passed"))
  expect_true(has(g, "--extsize 200 was added"))
  expect_true(has(g, "duplicate reads were not removed"))
  ## ChIP-seq with and without control
  chip <- modifyList(atac, list(mode = "ChIP-seq", paired = FALSE, useControl = TRUE))
  h <- app$methods_facts(chip)
  expect_true(has(h, "bdgcmp -m FE"))
  expect_true(has(h, "no mapping-quality or mitochondrial filter"))
  expect_false(has(h, "bamCoverage") || has(h, "-f BAMPE") || has(h, "MAPQ of at least 10") || has(h, "--extsize 200"))
  expect_true(has(app$methods_facts(modifyList(chip, list(useControl = FALSE))), "bamCoverage"))
  expect_true(has(app$methods_facts(modifyList(chip, list(cmdOptions = "--nomodel"))), "--extsize 147 was added"))
  expect_true(has(app$methods_facts(modifyList(chip, list(cmdOptions = "--broad"))), "broadPeak"))
})

test_that("MACS3 anchors: the values the facts state are in the code", {
  f <- EzAppMacs3$new()$methods_facts(list(mode = "ATAC-seq", paired = TRUE, shiftATAC = TRUE, annotatePeaks = TRUE,
                                           removeDuplicates = TRUE, genomeSize = 0, qValue = 0.05, cmdOptions = ""))
  for (v in c("--keep-dup all", "-f BAMPE", "80% of", "--extsize 200", "MAPQ of at least 10", "chrM", "--ATACshift",
              "-F 1024", "REMOVE_DUPLICATES=true", "--binSize 10", "CPM", "tssRegion -1000 to 1000", "nearestStart",
              "dustyScore", "getfasta")) expect_true(has(f, v), label = v)
  m <- codeOfFn(ezRun:::ezMethodMacs3)
  for (v in c("--keep-dup all", "\"-q\", param$qValue", "-f BAMPE", "gsize * 0.8", "param$genomeSize == 0",
              "--extsize 200", "--extsize 147", "gsub(\"--extsize 147\", \"--extsize 200\", opt)", "-m FE",
              "bedGraphToBigWig", "bedSort", "method = \"deepTools\"", "getfasta", "-name", "grepl(\"broad\", opt)")) {
    expect_match(m, v, fixed = TRUE)
  }
  b <- codeOfFn(ezRun:::atacBamProcess, ezRun:::filteroutBam)
  for (v in c("mapQ = 10", "c(\"M\", \"MT\", \"chrM\", \"chrMT\")", "--ATACshift", "\"-q\", mapQ")) {
    expect_match(b, v, fixed = TRUE)
  }
  expect_match(codeOfFn(ezRun:::removeDuplicatesFromBam), "-F 1024", fixed = TRUE)
  expect_match(codeOfFn(ezRun:::bamHasMarkedDuplicates), "MarkDuplicates|markdup", fixed = TRUE)
  expect_match(codeOfFn(ezRun:::dupBam), "REMOVE_DUPLICATES=", fixed = TRUE)
  bw <- codeOfFn(ezRun:::bam2bw)
  expect_match(bw, "--binSize 10 -of bigwig", fixed = TRUE)
  expect_match(bw, "--normalizeUsing CPM", fixed = TRUE)
  a <- codeOfFn(ezRun:::annotatePeaks)
  for (v in c("tssRegion = c(-1000, 1000)", "output = \"nearestStart\"", "FeatureLocForDistance = \"TSS\"",
              "multiple = FALSE", "dustyScore(seqs)", "decreasing = TRUE")) {
    expect_match(a, v, fixed = TRUE)
  }
})

test_that("PeakCombiner facts: unconditional, gated and carrying the run's values", {
  app <- EzAppPeakCombiner$new()
  f0 <- app$methods_facts(list())
  expect_gte(length(f0), 4)
  expect_true(has(f0, "minMQS 10"))
  expect_true(has(f0, "no normalization"))
  expect_false(has(f0, "longer than 5 characters"))
  f <- app$methods_facts(list(minSamples = 2L, skipExtraChr = TRUE))
  expect_true(has(f, "peak sets of at least 2 samples"))
  expect_true(has(f, "longer than 5 characters"))
  g <- app$methods_facts(list(minSamples = 3L, skipExtraChr = FALSE))
  expect_true(has(g, "at least 3 samples"))
  expect_false(has(g, "longer than 5 characters"))
})

test_that("PeakCombiner anchors: the values the facts state are in the code", {
  f <- EzAppPeakCombiner$new()$methods_facts(list(minSamples = 2L, skipExtraChr = TRUE))
  for (v in c("minMQS 10", "minOverlap 5", "allowMultiOverlap", "countMultiMappingReads", "primaryOnly", "ignoreDup FALSE",
              "strandSpecific 0", "reduce", "5 characters", "peakFraction", "fold enrichment")) expect_true(has(f, v), label = v)
  p <- codeOfFn(ezRun:::ezMethodCombinePeaks)
  for (v in c("support >= param$minSamples", "GenomicRanges::reduce(filtered, with.revmap = TRUE)",
              "nchar(as.character(seqnames(collapsed))) <= 5", "minMQS = 10", "minOverlap = 5",
              "allowMultiOverlap = TRUE", "countMultiMappingReads = TRUE", "primaryOnly = TRUE",
              "ignoreDup = FALSE", "strandSpecific = 0", "requireBothEndsMapped = FALSE",
              "fraction = FALSE", "readExtension5 = 0", "paste0(\"peak_\", seq_along(collapsed))",
              "max(filtered$fold_enrichment[i], na.rm = TRUE)", "mean(filtered$fold_enrichment[i], na.rm = TRUE)")) {
    expect_match(p, v, fixed = TRUE)
  }
  expect_match(p, "countStats\\[, 1\\]\\s*/\\s*countStats\\[, 2\\]")
  expect_match(codeOfFn(ezRun:::import_narrowPeak), "fold_enrichment = \"numeric\"", fixed = TRUE)
})
