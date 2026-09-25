## Guards on the LLM-written Methods text (R/methods-guards.R). Synthetic text only.

guardConfig <- paste(
  "param[['resolution']] = '0.6'", "param[['npcs']] = '20'", "param[['nfeatures']] = '3000'",
  "emptyDrops niters 100000", "window 2000", "threshold 0.5", sep = "\n")
guardAll <- paste(guardConfig, "Seurat_5.1.0", "R_4.6.0", "12,345 cells passed", "85.3% mapped", sep = "\n")
checkNum <- function(text, all = guardAll, n = 24)
  methods_check_numbers(text, guardConfig, all, sample_count = n)

test_that("planted result numbers, versions and parameters are flagged", {
  expect_identical(checkNum("In total 12,345 cells passed QC."), "12,345")
  expect_identical(checkNum("Of the reads, 85.3% of reads mapped."), "85.3%")
  expect_identical(checkNum("Clustering used Seurat 4.9.9."), "4.9.9")
  expect_identical(checkNum("Clusters were found at resolution 0.8."), "0.8")
})

test_that("configured values in other spellings, identifiers and small integers pass", {
  expect_length(checkNum("SCTransform selected 3,000 variable features."), 0)
  expect_length(checkNum("emptyDrops ran with 1e5 iterations."), 0)
  expect_length(checkNum("Windows of 2 kb were used."), 0)
  expect_length(checkNum("The threshold was [0.5], i.e. c(0.5)."), 0)
  expect_length(checkNum("Analyses ran in R version 4.6.0."), 0)
  expect_length(checkNum("R version 4.6.0 was used.", all = paste(guardConfig, "R version 4.6.0")), 0)
  expect_length(checkNum("The neighbour graph used PCs 1 to 20 (1:20, 1-20, 1–20)."), 0)
  expect_length(checkNum("Reads were aligned to GRCh38 from 10x Genomics libraries with log2FC and CD45."), 0)
  expect_length(checkNum("All 24 samples were processed in 3 batches."), 0)
})

test_that("number normalisation covers units, percentages and scientific forms", {
  cfg <- "a 0.05 b 1e-5 c 100000 d 3000000"
  expect_length(methods_check_numbers("5% and 10^-5 and 1E-5, 100k and 3 Mb.", cfg, cfg, 1), 0)
  expect_identical(methods_check_numbers("17% of 2e6", cfg, cfg, 1), c("17%", "2e6"))
  ## a year is a plain number: allowed only when the citation text (config) has it
  expect_identical(methods_check_numbers("Smith et al. (2019)", cfg, cfg, 1), "2019")
  expect_length(methods_check_numbers("Smith et al. (2019)", paste(cfg, "Smith 2019"), cfg, 1), 0)
})

test_that("claim detection: the gated.py selftest cases", {
  kw <- "decoupler|dorothea|progeny"
  expect_true(methodsClaims("Transcription-factor activities were inferred with decoupleR.", kw))
  expect_false(methodsClaims("decoupleR was not run.", kw))
  expect_false(methodsClaims("Clusters were found with Louvain.", kw))
})

test_that("a step described while its parameter was off is flagged, per class", {
  ## class, parameter(s) set off, a sentence claiming the step
  cases <- list(
    list("EzAppScSeurat", list(SingleR = "none"), "Cells were annotated with SingleR."),
    list("EzAppScSeurat", list(Azimuth = "none", AzimuthPanHuman = "true"), "Cells were mapped with Azimuth."),
    list("EzAppScSeuratCombine", list(integrationMethod = "none"), "Samples were integrated with Harmony."),
    list("EzAppScSeuratCombinedLabelClusters", list(tissue = ""), "Groups were scored with AUCell."),
    list("EzAppSpatialSeurat", list(spotClean = FALSE), "Spot swapping was removed with SpotClean."),
    list("EzAppSpatialSeuratSlides", list(batchCorrection = FALSE), "Slides were integrated with CCA anchors."),
    list("EzAppXeniumSeurat", list(rctdFile = "", rctdReference = "None"), "Cells were annotated with RCTD."),
    list("EzAppVisiumHDSeurat", list(rctdFile = "", rctdReference = "None"), "RCTD ran with spacexr."),
    list("EzAppDeseq2", list(useLfcShrink = FALSE), "Fold changes were shrunk with ashr."),
    list("EzAppEdger", list(runGO = FALSE), "GO terms were tested with clusterProfiler enricher."),
    list("EzAppLimma", list(runGO = FALSE), "GSEA was run on the ranked genes."),
    list("EzAppScMultiOmics", list(runWNN = FALSE), "Modalities were combined by WNN."),
    list("EzAppScSeuratCompare", list(replicateGrouping = ""), "Composition was tested with sccomp."),
    list("EzAppCellRanger", list(runVeloCyto = FALSE), "Velocity counts came from velocyto."),
    list("EzAppCellRangerMulti", list(keepBam = FALSE), "BAM files were converted to CRAM."),
    list("EzAppSpaceRanger", list(panelFile = ""), "The panel was passed as --feature-ref."),
    list("EzAppSTAR", list(twopassMode = FALSE), "STAR ran in two-pass mode."),
    list("EzAppBismark", list(generateBigWig = FALSE), "A bigWig file was written."),
    list("EzAppGatkDnaHaplotyper", list(knownSitesAvailable = FALSE), "Qualities were recalibrated with BaseRecalibrator."),
    list("EzAppJoinGenoTypes", list(recalibrateVariants = FALSE, recalibrateInDels = FALSE), "SNPs were filtered by VQSR."),
    list("EzAppMetaPhlAn", list(estimateReadCounts = FALSE), "MetaPhlAn ran with -t rel_ab_w_read_stats."),
    list("EzAppFastqc", list(generate_ai_summary = FALSE, per_section_ai_summaries = FALSE), "A language model summarised the report."),
    list("EzAppSamsa2", list(paired = FALSE), "Pairs were merged with PEAR."),
    list("EzAppCellBender", list(gpu = 0), "CellBender ran with --cuda.")
  )
  expect_setequal(vapply(cases, `[[`, "", 1), names(METHODS_OFFSTEP_RULES))
  for (x in cases) {
    cls <- x[[1]]; off <- x[[2]]; claim <- x[[3]]
    expect_length(methods_check_offsteps(claim, cls, off), 1)
    expect_length(methods_check_offsteps(paste(claim, "Clusters were found."), cls, off), 1)
    on <- lapply(off, function(v) "TRUE")
    expect_length(methods_check_offsteps(claim, cls, on), 0)
    expect_length(methods_check_offsteps(sub("\\.$", " was not used.", claim), cls, off), 0)
  }
  ## a rule on two parameters needs both off; an unknown parameter is not off
  expect_length(methods_check_offsteps("RCTD ran.", "EzAppXeniumSeurat", list(rctdFile = "", rctdReference = "x")), 0)
  expect_length(methods_check_offsteps("Cells were annotated with SingleR.", "EzAppScSeurat", list(name = "x")), 0)
  expect_length(methods_check_offsteps("Cells were annotated with SingleR.", "EzAppScSeurat", list()), 0)
  ## Pan-Human Azimuth on is not a claim about the Azimuth parameter
  expect_length(methods_check_offsteps("Cells were annotated with Pan-Human Azimuth.", "EzAppScSeurat",
                                       list(Azimuth = "none", AzimuthPanHuman = TRUE)), 0)
})
