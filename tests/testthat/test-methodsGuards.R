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
  ## a glued version is skipped whole, not read as its first digits (2.7 of 2.7.11b)
  expect_length(checkNum("Reads were aligned with STAR 2.7.11b and HISAT2 2.2.1a."), 0)
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

test_that("true sentences about an off step are not flagged, and the step's claims still are", {
  ## class, parameters set off, true sentences (guards-7, eval-2, facts-2: archived texts), a claim
  cases <- list(
    list("EzAppScSeurat", list(SingleR = "none", Azimuth = "none"),
         c("SingleR was set to 'none'.", "SingleR annotation was set to \"none\".", "Azimuth and SingleR annotation were both set to none.",
           "Cell-type annotation used Azimuth (set to 'none') and SingleR (set to 'none').", "SingleR = none; Azimuth was switched off."),
         "Cells were annotated with SingleR."),
    list("EzAppSpatialSeurat", list(Azimuth = "none", spotClean = FALSE),
         c("Markers were queried against the Enrichr databases Azimuth_Cell_Types_2021 and PanglaoDB_Augmented_2021.",
           "Markers were queried against the Azimuth Cell Types 2021 and Human Gene Atlas databases.",
           "The Azimuth reference-based label transfer was set to none.", "spotClean = false;"),
         "Cells were mapped to the Azimuth reference."),
    list("EzAppScSeurat", list(tissue = ""), "Markers were queried against the Enrichr library CellMarker_Augmented_2021.",
         "Cells were scored with AUCell against CellMarker gene sets."),
    list("EzAppSTAR", list(twopassMode = FALSE, barcodePattern = ""),
         c("STAR ran in one-pass mode (--twopassMode None).", "Reads were aligned in one-pass mode.",
           "Reads with fewer than 10 UMIs were counted by STARsolo."),
         c("STAR ran in two-pass mode (--twopassMode Basic).", "UMIs were extracted with umi_tools.")),
    list("EzAppXeniumSeurat", list(doSPLIT = FALSE, coocFdr = FALSE),
         c("The data were split by sample before clustering.", "Marker p-values were corrected with the Benjamini-Hochberg method."),
         c("Spatial purification was performed using SPLIT (version 0.1.2).", "The purified counts were re-annotated with RCTD after SPLIT.",
           "Co-occurrence p-values were corrected with the Benjamini-Hochberg method.")),
    list("EzAppSpatialSeuratSlides", list(batchCorrection = FALSE),
         "The reported UMAP and clusters came from the integrated assay.",
         c("Slides were integrated with Seurat CCA anchors.", "Slides were combined with FindIntegrationAnchors and IntegrateData.")),
    list("EzAppGatkDnaHaplotyper", list(markDuplicates = FALSE), "markDuplicates = false;", "Duplicates were marked with MarkDuplicates.")
  )
  for (x in cases) {
    for (t in x[[3]]) expect_identical(methods_check_offsteps(t, x[[1]], x[[2]]), character(0), label = t)
    for (t in x[[4]]) expect_length(methods_check_offsteps(t, x[[1]], x[[2]]), 1)
  }
  ## a false discovery rate is not a negation: the one real arm-A hit (eval-2) stays flagged
  expect_length(methods_check_offsteps("The false discovery rate threshold for over-representation analysis was set to 0.05.",
                                       "EzAppDeseq2", list(runGO = FALSE)), 1)
})

## Stand-in for llm_write_methods that writes `text` to methods.md and exits with `exit`;
## call n writes text[n] (the last one once they run out). Like the deployed CLI
## (llm_write_methods 0.9:320-322) it APPENDS to --output, after a 72-dash rule when the
## file is not empty. It counts its calls and keeps the last task file, so retries can be checked.
withTextWriter <- function(text, code, exit = 0) {
  bin <- tempfile("bin"); dir.create(bin)
  for (i in seq_along(text)) writeLines(text[i], file.path(bin, paste0("body", i, ".txt")))
  writeLines(c(
    "#!/bin/sh",
    paste0("printf '%s\\n' \"$@\" > ", shQuote(file.path(bin, "args_last.txt"))),
    paste0("echo call >> ", shQuote(file.path(bin, "calls.txt"))),
    "while [ $# -gt 0 ]; do [ \"$1\" = --output ] && out=$2; [ \"$1\" = --task-file ] && task=$2; shift; done",
    paste0("n=$(wc -l < ", shQuote(file.path(bin, "calls.txt")), ")"),
    paste0("body=", shQuote(bin), "/body$n.txt; [ -f \"$body\" ] || body=", shQuote(file.path(bin, paste0("body", length(text), ".txt")))),
    "[ -s \"$out\" ] && printf '\\n%s\\n\\n' ------------------------------------------------------------------------ >> \"$out\"",
    "cat \"$body\" >> \"$out\"",
    paste0("cp \"$task\" ", shQuote(file.path(bin, "task_last.txt"))),
    paste0("exit ", exit)
  ), file.path(bin, "llm_write_methods"))
  Sys.chmod(file.path(bin, "llm_write_methods"), "755")
  oldPath <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(bin, oldPath, sep = ":"))
  on.exit(Sys.setenv(PATH = oldPath))
  code(bin)
}
nCalls <- function(bin) length(readLines(file.path(bin, "calls.txt")))

## A result dir with parameters.tsv, one job script and its log, as write_methods reads it.
fakeResultDir <- function(params, ezrun = paste0("ezRun_", utils::packageVersion("ezRun"))) {
  res <- tempfile("res"); dir.create(file.path(res, "scripts"), recursive = TRUE)
  writeLines(paste(names(params), params, sep = "\t"), file.path(res, "parameters.tsv"))
  writeLines(c("#!/bin/bash", paste0("param[['", names(params), "']] = '", params, "'")),
             file.path(res, "scripts", "ScSeurat_S1.sh"))
  writeLines(c("loading ezRun", paste(ezrun, "Seurat_5.1.0"), "Cells after QC: 4,321"),
             file.path(res, "scripts", "ScSeurat_S1.sh_sushiID1_2026-01-01--00-00-00_o.log"))
  res
}
readMd <- function(out) paste(readLines(file.path(out, "methods.md")), collapse = "\n")

test_that("a failing writer still gives methods.md, with the template Declaration", {
  withTextWriter("ignored", exit = 1, function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzApp$new()$write_methods(output_dir = out, analysis_name = "T"),
                   "template fallback \\(llm_write_methods failed with exit code 1\\)")
    md <- readMd(out)
    expect_match(md, "The analysis was run with the FGCZ SUSHI app EzApp.", fixed = TRUE)
    expect_match(md, METHODS_TEMPLATE_DECLARATION, fixed = TRUE)
    expect_no_match(md, "generated by")
  })
})

test_that("a flagged number is retried once with the tokens named, then templated", {
  withTextWriter("In total 12,345 cells passed QC.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzApp$new()$write_methods(output_dir = out, analysis_name = "T"),
                   "guards: 1 numbers, 0 steps")
    expect_equal(nCalls(bin), 2)
    expect_match(paste(readLines(file.path(bin, "task_last.txt")), collapse = "\n"),
                 "not in the run's configuration; remove them or the sentence that states them: 12,345.",
                 fixed = TRUE)
    md <- readMd(out)
    expect_no_match(md, "12,345")
    expect_match(md, METHODS_TEMPLATE_DECLARATION, fixed = TRUE)
  })
})

test_that("a retry returns the second draft only: the appending writer's file is removed first", {
  withTextWriter(c("In total 12,345 cells passed QC.", "Reads were aligned with 3 mismatches."), function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzApp$new()$write_methods(output_dir = out, analysis_name = "T"), "retrying")
    expect_equal(nCalls(bin), 2)
    md <- readMd(out)
    expect_match(md, "Reads were aligned with 3 mismatches.", fixed = TRUE)
    expect_no_match(md, "12,345")
    expect_no_match(md, "-{72}")
    expect_match(md, "generated by", fixed = TRUE)
  })
})

test_that("backticks the writer adds are removed from the Description", {
  withTextWriter("Reads were sorted by name with `samtools sort -n` and counted with ``featureCounts``.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    EzApp$new()$write_methods(output_dir = out, analysis_name = "T")
    md <- readMd(out)
    expect_match(md, "Reads were sorted by name with samtools sort -n and counted with featureCounts.", fixed = TRUE)
    expect_no_match(md, "`", fixed = TRUE)
  })
})

test_that("clean text is kept and the Declaration names LLM_CALLER_MODEL, else the fallback", {
  withTextWriter("Reads were aligned with 3 mismatches.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    withr::with_envvar(c(LLM_CALLER_MODEL = "Model-X"),
                       EzApp$new()$write_methods(output_dir = out, analysis_name = "T"))
    expect_equal(nCalls(bin), 1)
    md <- readMd(out)
    expect_match(md, "Reads were aligned with 3 mismatches.", fixed = TRUE)
    expect_match(md, "generated by Model-X based on", fixed = TRUE)
    withr::with_envvar(c(LLM_CALLER_MODEL = NA),
                       EzApp$new()$write_methods(output_dir = out, analysis_name = "T"))
    expect_match(readMd(out), paste0("generated by ", METHODS_LLM_MODEL_NAME, " based on"), fixed = TRUE)
  })
})

test_that("an off step is retried then templated with the run's parameters, facts and version", {
  res <- fakeResultDir(c(name = "ScSeurat", cores = "8", SingleR = "none", npcs = "20",
                         CyteTypeR.apiKey = "hidden"))
  withTextWriter("Cells were annotated with SingleR using 20 PCs.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzAppScSeurat$new()$write_methods(file.path(res, "scripts"), out, "T",
                                                     example_script = "ScSeurat_S1.sh", sample_count = 3),
                   "guards: 0 numbers, 1 steps")
    expect_equal(nCalls(bin), 2)
    expect_match(paste(readLines(file.path(bin, "task_last.txt")), collapse = "\n"),
                 "do not describe them: singler (SingleR off).", fixed = TRUE)
    md <- readMd(out)
    expect_match(md, paste0("app ScSeurat (ezRun ", utils::packageVersion("ezRun"),
                            ") on 3 samples. Parameters: SingleR = none; npcs = 20."), fixed = TRUE)
    expect_no_match(md, "hidden|cores")
    expect_match(md, "The random seed was set to 38", fixed = TRUE)   # a fact
    expect_match(md, "https://doi.org/10.12688/f1000research.73600.2", fixed = TRUE)  # a candidate
  })
})

test_that("config.csv is read for the example sample only", {
  res <- tempfile("res")
  for (s in c("CTCL_PBMCs3", "PBMCs3", "P7")) {
    dir.create(file.path(res, s), recursive = TRUE)
    writeLines("[gene-expression]", file.path(res, s, "config.csv"))
  }
  pick <- function(script) basename(dirname(methodsConfigCsv(res, script)))
  expect_identical(pick("SingleCell_CTCL_PBMCs3_CTCL_blood_input_114700.sh"), "CTCL_PBMCs3")
  expect_identical(pick("SingleCell_PBMCs3.sh"), "PBMCs3")
  expect_identical(pick("SingleCell_P7_run_1.sh"), "P7")
  expect_identical(pick("SingleCell_unknown.sh"), "CTCL_PBMCs3")   # no match: first only
  expect_length(methodsConfigCsv(res, NULL), 3)
  expect_length(methodsConfigCsv(tempfile(), "x.sh"), 0)
})

test_that("merged with the facts guard: no facts for another ezRun version, empty values unlisted", {
  res <- fakeResultDir(c(name = "ScSeurat", SingleR = "none", tissue = '""'), ezrun = "ezRun_0.0.1")
  withTextWriter("Cells were annotated with SingleR.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    suppressMessages(EzAppScSeurat$new()$write_methods(file.path(res, "scripts"), out, "T",
                                                       example_script = "ScSeurat_S1.sh", sample_count = 1))
    md <- readMd(out)
    expect_match(md, METHODS_TEMPLATE_DECLARATION, fixed = TRUE)
    expect_no_match(md, "The random seed was set to 38")   # facts withheld: run used ezRun 0.0.1
    expect_no_match(md, "tissue =")                          # "" is empty, not a value
  })
})

test_that("numbers the logs hold pass, and the corpus false positives stay unflagged", {
  res <- fakeResultDir(c(name = "ScSeurat", SingleR = "none"))
  withTextWriter("After quality control 4,321 cells were kept.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    EzAppScSeurat$new()$write_methods(file.path(res, "scripts"), out, "T",
                                      example_script = "ScSeurat_S1.sh", sample_count = 1)
    expect_equal(nCalls(bin), 1)
    expect_match(readMd(out), "4,321 cells were kept", fixed = TRUE)
  })
  expect_length(methods_check_offsteps("Markers were queried against the Enrichr library Azimuth_2023.",
                                       "EzAppScSeurat", list(Azimuth = "none")), 0)
  expect_length(methods_check_offsteps("Reads were aligned in single-pass mode (twopassMode = None).",
                                       "EzAppSTAR", list(twopassMode = FALSE)), 0)
  expect_length(methods_check_offsteps("Cells were mapped to the Azimuth reference.",
                                       "EzAppScSeurat", list(Azimuth = "none")), 1)  # positive control
})

test_that("compute resources are asked out once, then their sentences dropped", {
  expect_setequal(methods_check_resources("It used four cores, 8 threads and 12 GB of RAM in dataset mode with 100 GB scratch."),
                  c("four cores", "8 threads", "12 gb of ram", "dataset mode", "scratch"))
  expect_length(methods_check_resources("Clusters were annotated using the FGCZ-hosted language model with 30 PCs."), 0)
  ## the phrasings that reached final texts (guards-2, ronald-1), and their neighbours
  for (x in c("Assembly was run with eight threads.", "It used sixty-four threads.", "It used twenty four cores.",
              "Local memory was set to 100 GB, and 8 local cores were used.", "A memory of 16 GB was requested.",
              "Reads were sorted with samtools sort -n -m 3500M -@ 4.", "BWA ran with --threads 8.",
              "samtools used a memory limit of 2048 MB per thread.", "It ran with 16 parallel threads.",
              "The pipeline completed successfully.", "All jobs finished successfully.", "The job ran to completion."))
    expect_gt(length(methods_check_resources(x)), 0, label = x)
  for (x in c("STAR aligned the reads in one-pass mode.", "Two samples were compared.",
              "The first 4 principal components were used.", "Genes expressed in at least three cells were kept.",
              "Reads were counted with featureCounts -t exon -g gene_id.", "Twenty clusters were found at resolution 0.5.",
              "The core promoter set was used.", "A 5 Mb window was used."))
    expect_length(methods_check_resources(x), 0)
  withTextWriter("Reads were aligned with STAR. The job used 8 threads.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzApp$new()$write_methods(output_dir = out, analysis_name = "T"), "1 resources, 0 vendor steps, 0 style; retrying")
    expect_equal(nCalls(bin), 2)
    expect_match(paste(readLines(file.path(bin, "task_last.txt")), collapse = "\n"), "compute resources", fixed = TRUE)
    expect_match(readMd(out), "generated by", fixed = TRUE)   # delivered, not templated
    expect_match(readMd(out), "Reads were aligned with STAR.", fixed = TRUE)
    expect_no_match(readMd(out), "8 threads")
  })
})

test_that("after the retry only the flagged sentences are dropped; the template only when nothing is left", {
  raw <- "Cells were clustered with Louvain. Cells were annotated with SingleR. SingleR was not run.\n\n## References\nX https://doi.org/1"
  kept <- methodsDropSentences(raw, steps = "singler")
  expect_match(kept, "Cells were clustered with Louvain. SingleR was not run.", fixed = TRUE)  # negated one stays
  expect_match(kept, "https://doi.org/1", fixed = TRUE)
  expect_identical(methodsDropSentences("Version 2016 had 16 samples. It had 160 cells.", numbers = "16"),
                   "It had 160 cells.")
  expect_true(157 %in% methodsConfigValues("/Variants/dbsnp.157.gencode_compatible.vcf"))
  res <- fakeResultDir(c(name = "ScSeurat", SingleR = "none"))
  withTextWriter("Cells were clustered with Louvain. Cells were annotated with SingleR.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzAppScSeurat$new()$write_methods(file.path(res, "scripts"), out, "T",
                                                     example_script = "ScSeurat_S1.sh", sample_count = 1),
                   "dropped the sentences with 0 numbers, 1 steps")
    md <- readMd(out)
    expect_match(md, "Cells were clustered with Louvain.", fixed = TRUE)
    expect_no_match(md, "annotated with SingleR")
    expect_match(md, "generated by", fixed = TRUE)
  })
})

## Martian stage lines as Cell Ranger count (fork0) and multi (fork_<sample>, a stage
## scheduled in two subpipelines and run in one) write them; trimmed from ER13 and CTCL_PBMCs3.
martianLog <- c(
  "2026-08-21 13:58:33 [runtime] (ready)           ID.S.SC_RNA_COUNTER_CS.SC_MULTI_CORE.COUNT_ANALYZER.SC_RNA_ANALYZER.RUN_PCA",
  "2026-08-21 13:58:33 [runtime] (run:local)       ID.S.SC_RNA_COUNTER_CS.SC_MULTI_CORE.COUNT_ANALYZER.SC_RNA_ANALYZER.RUN_PCA.fork0.split",
  "2026-08-21 13:58:33 [runtime] (split_complete)  ID.S.SC_RNA_COUNTER_CS.SC_MULTI_CORE.COUNT_ANALYZER.SC_RNA_ANALYZER.RUN_PCA",
  "2026-09-04 02:26:41 [runtime] (ready)           ID.S.SC_MULTI_CS.SC_MULTI_CORE.COUNT_ANALYZER.SC_RNA_ANALYZER.RUN_UMAP",
  "2026-09-04 02:26:58 [runtime] (ready)           ID.S.SC_MULTI_CS.SC_MULTI_CORE.SAMPLE_ANALYZER.SC_RNA_ANALYZER.RUN_UMAP",
  "2026-09-04 02:26:58 [runtime] (run:local)       ID.S.SC_MULTI_CS.SC_MULTI_CORE.SAMPLE_ANALYZER.SC_RNA_ANALYZER.RUN_UMAP.fork_CTCL_PBMCs3-cellRanger.chnk0.main",
  "2026-08-21 13:59:07 [runtime] (ready)           ID.S.SC_RNA_COUNTER_CS.SC_MULTI_CORE.COUNT_ANALYZER.SC_RNA_ANALYZER.RUN_HIERARCHICAL_CLUSTERING",
  "2026-08-21 13:59:07 [runtime] (ready)           ID.S.SC_MULTI_CS.SC_MULTI_CORE.CELL_ANNOTATE.WRITE_CELL_TYPES_H5")

test_that("vendor stages: ran when any subpipeline ran them, skipped when only scheduled", {
  f <- tempfile(fileext = "_o.log"); writeLines(c("loading ezRun", martianLog), f)
  st <- methodsMartianStages(f)
  expect_identical(st$ran, c("PCA", "UMAP"))
  expect_identical(st$skipped, c("hierarchical clustering", "cell type annotation"))
  expect_null(methodsMartianStages(character(0)))
  g <- tempfile(); writeLines("no martian here", g); expect_null(methodsMartianStages(g))
  expect_identical(methods_check_vendor("Secondary analysis included PCA and hierarchical clustering.", st),
                   "hierarchical clustering")
  expect_identical(methods_check_vendor("The pipeline also performed cell typing.", st), "cell type annotation")
  expect_length(methods_check_vendor("Cell annotation was skipped because the reference is not supported. UMAP was computed.", st), 0)
})

test_that("the run summary states how the samples ran, the reference and the vendor steps", {
  rb <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"
  s <- methodsRunSummary(list(process_mode = "DATASET", refBuild = rb), NULL, 12)
  expect_identical(s, c("The samples were analysed together, in one job.",
                        "Reference: Homo sapiens, GENCODE GRCh38.p14, annotation release 48."))
  expect_match(methodsRunSummary(list(process_mode = "SAMPLE"), NULL, 12),
               "Each of the 12 samples was analysed separately, in its own job", fixed = TRUE)
  expect_identical(methodsRunSummary(list(), data.frame(Species = c("Mus musculus", "Mus musculus")), 1),
                   "Organism: Mus musculus.")
  expect_length(methodsRunSummary(list(process_mode = "SAMPLE"), data.frame(Name = "a"), 1), 0)
  f <- tempfile(fileext = "_o.log"); writeLines(martianLog, f)
  s <- methodsRunSummary(list(), NULL, 1, f)
  expect_identical(s, c("The vendor pipeline ran these steps: PCA, UMAP.",
                        "It did not run these steps, so do not describe them: hierarchical clustering, cell type annotation."))
})

test_that("framework names, files, paths, flags and sweeping [not recorded] are flagged; plain text is not", {
  for (x in c("It was run within the ezRun framework.", "SUSHI launched the job.", "Counts were read from counts.tsv.",
              "Reads in /srv/gstore/projects/p1 were used.", "fastp was run with --length_required 25.",
              "[not recorded: any additional FastQC parameters beyond the adapter list].",
              "Any additional settings used internally by these steps are [not recorded].",
              "The settings for chimeric fragments were [not recorded] beyond the behaviour described above.",
              "Additional assembly parameters were [not recorded].", "The counting parameters were otherwise [not recorded]."))
    expect_gt(length(methods_check_style(x)), 0, label = x)
  for (x in c("Reads were aligned to GRCh38.p14 with STAR 2.7.11b.", "Genes and/or transcripts were counted.",
              "Reads were counted with featureCounts -t exon.", "The minimum UMI count was [not recorded].",
              "Libraries were prepared with the 10x Genomics 5' kit (v2).", "Log2 fold changes -- shrunken -- were used.",
              "Bins were further filtered with scater isOutlier at a threshold [not recorded]."))
    expect_length(methods_check_style(x), 0)
})

test_that("write_methods passes the run summary, drops skipped vendor steps and the framework, adds the ezRun line", {
  res <- tempfile("res"); dir.create(file.path(res, "scripts"), recursive = TRUE)
  writeLines(c("process_mode\tSAMPLE", "refBuild\tHomo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"),
             file.path(res, "parameters.tsv"))
  writeLines("#!/bin/bash", file.path(res, "scripts", "CR_S1.sh"))
  writeLines(c("ezRun_3.23.2", martianLog), file.path(res, "scripts", "CR_S1.sh_sushiID1_2026-01-01--00-00-00_o.log"))
  withTextWriter(c("Reads were aligned with STAR. Secondary analysis included hierarchical clustering. It ran within the ezRun framework.",
                   "Reads were aligned with STAR. Secondary analysis included hierarchical clustering."), function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzApp$new()$write_methods(file.path(res, "scripts"), out, "T", example_script = "CR_S1.sh",
                                             sample_count = 3),
                   "dropped the sentences with 0 numbers, 0 steps, 0 resources, 1 vendor steps")
    expect_equal(nCalls(bin), 2)
    task <- paste(readLines(file.path(bin, "task_last.txt")), collapse = "\n")
    expect_match(task, "did not run these steps in this job; do not describe them: hierarchical clustering.", fixed = TRUE)
    expect_match(task, "Remove the framework names (ezRun, SUSHI)", fixed = TRUE)
    expect_true(file.path(out, "run_summary.txt") %in% readLines(file.path(bin, "args_last.txt")))
    rs <- readLines(file.path(out, "run_summary.txt"))
    expect_true("Each of the 3 samples was analysed separately, in its own job, with the same settings." %in% rs)
    expect_true("It did not run these steps, so do not describe them: hierarchical clustering, cell type annotation." %in% rs)
    md <- readMd(out)
    expect_match(md, "Reads were aligned with STAR.", fixed = TRUE)
    expect_no_match(md, "hierarchical|ezRun framework")
    expect_match(md, "The analysis was run with the FGCZ ezRun package, version 3.23.2 (https://github.com/uzh/ezRun).", fixed = TRUE)
  })
})

test_that("a framework name the retry keeps is delivered, not dropped or templated", {
  withTextWriter("Counts were tested with DESeq2 within the ezRun framework.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzApp$new()$write_methods(output_dir = out, analysis_name = "T"), "1 style; retrying")
    expect_equal(nCalls(bin), 2)
    expect_match(readMd(out), "Counts were tested with DESeq2 within the ezRun framework.", fixed = TRUE)
    expect_match(readMd(out), "generated by", fixed = TRUE)
  })
})

test_that("a sweeping [not recorded] the retry keeps is dropped, the rest of the text delivered", {
  withTextWriter("Reads were assembled with hifiasm. The settings for any additional assembly parameters were [not recorded].", function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzApp$new()$write_methods(output_dir = out, analysis_name = "T"), "dropped the sentences")
    expect_equal(nCalls(bin), 2)
    expect_match(readMd(out), "Reads were assembled with hifiasm.", fixed = TRUE)
    expect_no_match(readMd(out), "not recorded")
    expect_match(readMd(out), "generated by", fixed = TRUE)
  })
})

test_that("a single job's samples: the samples parameter, and for a comparison the two groups only", {
  inp <- data.frame(Name = paste0("S", 1:6), `Condition [Factor]` = c("A", "A", "B", "B", "C", "C"), check.names = FALSE)
  expect_identical(methodsSamplesRun(list(), inp), 6L)
  expect_identical(methodsSamplesRun(list(samples = "S1,S2,S3"), inp), 3L)
  expect_identical(methodsSamplesRun(list(grouping = "Condition", sampleGroup = "A", refGroup = "B"), inp), 4L)
  expect_identical(methodsSamplesRun(list(grouping = "Condition", sampleGroup = "A", refGroup = "B", refGroupBaseline = "",
                                          samples = "S1,S3,S5"), inp), 2L)
  expect_identical(methodsSamplesRun(list(), NULL), 0L)
})

test_that("a [not recorded] marker is not a negation: the claim in its sentence is still checked (round 8)", {
  st <- list(ran = "PCA", skipped = "cell type annotation")
  expect_identical(methods_check_vendor("Cell type annotation was performed with settings [not recorded].", st), "cell type annotation")
  expect_identical(methods_check_vendor("Cell type annotation was also performed within the pipeline, with settings [not recorded].", st), "cell type annotation")
  expect_length(methods_check_vendor("Cell type annotation was not performed.", st), 0)   # a real negation still counts
  expect_length(methods_check_offsteps("Cells were annotated with SingleR, settings [not recorded].", "EzAppScSeurat", list(SingleR = "none")), 1)
  kept <- methodsDropSentences("UMAP was run. Cells were annotated with SingleR [not recorded].", steps = "singler")
  expect_identical(kept, "UMAP was run.")
})

test_that("FastQC: k-mer analysis is always flagged, the app never enables the Kmer Content module (round 8)", {
  p <- list(paired = "false")
  expect_length(methods_check_offsteps("FastQC was run with k-mer analysis enabled at a k-mer length of 7.", "EzAppFastqc", p), 1)
  expect_length(methods_check_offsteps("Reads were assessed with FastQC, with k-mer analysis enabled.", "EzAppFastqc", p), 1)
  expect_length(methods_check_offsteps("The Kmer Content module was disabled by default, so no k-mer analysis was done.", "EzAppFastqc", p), 0)
  expect_length(methods_check_offsteps("A k-mer length of 7 was set.", "EzAppFastqc", p), 0)
})
