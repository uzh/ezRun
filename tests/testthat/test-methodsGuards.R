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

## Stand-in for llm_write_methods that writes `text` as methods.md and exits with `exit`;
## it counts its calls and keeps the last task file, so retries can be checked.
withTextWriter <- function(text, code, exit = 0) {
  bin <- tempfile("bin"); dir.create(bin)
  writeLines(text, file.path(bin, "body.txt"))
  writeLines(c(
    "#!/bin/sh",
    paste0("echo call >> ", shQuote(file.path(bin, "calls.txt"))),
    "while [ $# -gt 0 ]; do [ \"$1\" = --output ] && out=$2; [ \"$1\" = --task-file ] && task=$2; shift; done",
    paste0("cp ", shQuote(file.path(bin, "body.txt")), " \"$out\""),
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

test_that("compute resources are asked out once, and a text that keeps them is still delivered", {
  expect_setequal(methods_check_resources("It used four cores, 8 threads and 12 GB of RAM in dataset mode with 100 GB scratch."),
                  c("8 threads", "12 gb of ram", "dataset mode", "scratch"))
  expect_length(methods_check_resources("Clusters were annotated using the FGCZ-hosted language model with 30 PCs."), 0)
  withTextWriter("Reads were aligned with 8 threads.", function(bin) {
    out <- tempfile("out"); dir.create(out)
    expect_message(EzApp$new()$write_methods(output_dir = out, analysis_name = "T"), "1 resources; retrying")
    expect_equal(nCalls(bin), 2)
    expect_match(paste(readLines(file.path(bin, "task_last.txt")), collapse = "\n"), "compute resources", fixed = TRUE)
    expect_match(readMd(out), "generated by", fixed = TRUE)   # delivered, not templated
  })
})
