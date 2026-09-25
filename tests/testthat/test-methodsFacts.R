## methods_facts(): fixed behaviour of an app's code that its parameter form does
## not show (seed, clustering algorithm, multiple-testing correction, ...). The
## Methods LLM only reads the job script and logs, so without this it has to write
## [not recorded] for these or guess them. methods_description() hands the facts to
## llm_write_methods as one more evidence file.

## Stand-in for the AI/llm_methods_caller binary: records its arguments and writes
## a minimal methods.md, so the plumbing is tested without calling the LLM.
withStubWriter <- function(code) {
  bin <- tempfile("bin"); dir.create(bin)
  argsFile <- file.path(bin, "args.txt")
  writeLines(c(
    "#!/bin/sh",
    paste0("printf '%s\\n' \"$@\" > ", shQuote(argsFile)),
    "while [ $# -gt 0 ]; do [ \"$1\" = --output ] && out=$2; shift; done",
    "printf 'Methods text.\\n' > \"$out\""
  ), file.path(bin, "llm_write_methods"))
  Sys.chmod(file.path(bin, "llm_write_methods"), "755")
  oldPath <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(bin, oldPath, sep = ":"))
  on.exit(Sys.setenv(PATH = oldPath))
  code(argsFile)
}

test_that("the base app has no facts and every drafted app has some", {
  expect_identical(EzApp$new()$methods_facts(), character(0))
  drafted <- list(EzAppScSeurat, EzAppDeseq2, EzAppEdger, EzAppLimma, EzAppHomerDiffPeaks,
                  EzAppCellBender, EzAppScSeuratCombine, EzAppScSeuratCombinedLabelClusters,
                  EzAppScSeuratCompare, EzAppScMultiOmics, EzAppVeloCyto, EzAppSpatialSeurat,
                  EzAppSpatialSeuratSlides, EzAppSpatialSeuratHD, EzAppVisiumHDSeurat, EzAppXeniumSeurat,
                  ## CellRanger and CellRangerMulti are left out: every fact they have is gated
                  ## on TenXLibrary, so they return none without the run's parameters.
                  EzAppCellRangerARC, EzAppSpaceRanger, EzAppXeniumQC, EzAppKraken, EzAppMetaPhlAn,
                  EzAppSamsa2, EzAppHifiasm, EzAppJoinGenoTypes, EzAppSTAR, EzAppBismark, EzAppKallisto,
                  EzAppFeatureCounts, EzAppFastqc, EzAppFastqScreen, EzAppFlash, EzAppCrisprScreenQC,
                  EzAppGatkDnaHaplotyper)
  for (cls in drafted) {
    facts <- cls$new()$methods_facts()
    expect_type(facts, "character")
    expect_gt(length(facts), 0)
    expect_false(any(is.na(facts) | !nzchar(facts)))
  }
})

test_that("methods_description passes app_facts.txt only when there are facts", {
  withStubWriter(function(argsFile) {
    script <- tempfile(fileext = ".sh"); writeLines("echo job", script)

    out <- tempfile("out"); dir.create(out)
    EzAppScSeurat$new()$methods_description(script, character(0), 1, out)
    args <- readLines(argsFile)
    factsFile <- file.path(out, "app_facts.txt")
    expect_true(file.exists(factsFile))
    scriptsAt <- which(args == "--scripts")
    expect_true(factsFile %in% args[-seq_len(scriptsAt)])
    expect_true(all(EzAppScSeurat$new()$methods_facts() %in% readLines(factsFile)))
    expect_match(readLines(factsFile)[1], "EzAppScSeurat")

    out2 <- tempfile("out"); dir.create(out2)
    EzApp$new()$methods_description(script, character(0), 1, out2)
    expect_false(file.exists(file.path(out2, "app_facts.txt")))
    expect_false(any(grepl("app_facts.txt", readLines(argsFile), fixed = TRUE)))
  })
})

test_that("gated facts follow the job's parameters", {
  app <- EzAppScSeurat$new()
  on  <- list(refBuild = "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03",
              computePathwayTFActivity = TRUE, estimateAmbient = TRUE, SingleR = "HumanPrimaryCellAtlasData")
  off <- modifyList(on, list(computePathwayTFActivity = FALSE, estimateAmbient = FALSE, SingleR = "none"))
  expect_true(any(grepl("decoupleR", app$methods_facts(on))))
  expect_false(any(grepl("decoupleR", app$methods_facts(off))))
  expect_false(any(grepl("DecontX", app$methods_facts(off))))
  expect_false(any(grepl("SingleR against", app$methods_facts(off))))
  ## unknown species or no parameters at all: species- and param-gated facts are dropped
  expect_false(any(grepl("decoupleR|cyclone", app$methods_facts(list()))))
  expect_false(any(grepl("cyclone", app$methods_facts(modifyList(on, list(refBuild = "Danio_rerio/Ensembl/GRCz11"))))))
})

test_that("methods_param reads parameters.tsv and fills app defaults", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c("name\tScSeurat", "computePathwayTFActivity\tfalse", "refBuild\tHomo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"), f)
  p <- EzAppScSeurat$new()$methods_param(f)
  expect_identical(p$computePathwayTFActivity, FALSE)
  expect_equal(p$nmad, 3)          # appDefault, not in the file
  expect_identical(EzApp$new()$methods_param(file.path(tempdir(), "absent.tsv")), list())
  ## SUSHI writes empty values as a literal "": they must not count as set
  g <- tempfile(fileext = ".tsv")
  writeLines(c("name\tCellRangerMulti", "secondRef\t\"\"", "controlSeqs\t\"\""), g)
  q <- EzAppCellRangerMulti$new()$methods_param(g)
  expect_false(ezIsSpecified(q$secondRef))
  expect_false(ezIsSpecified(q$controlSeqs))
})

test_that("write_methods demotes headings the model wrote into the description", {
  withStubWriter(function(argsFile) {
    ## stub writes a description that starts with its own "## Methods" heading
    stub <- file.path(dirname(argsFile), "llm_write_methods")
    writeLines(c("#!/bin/sh",
                 "while [ $# -gt 0 ]; do [ \"$1\" = --output ] && out=$2; shift; done",
                 "printf '## Methods\\n\\nReads were aligned.\\n\\n## Alignment\\n\\nWith STAR.\\n' > \"$out\""), stub)
    out <- tempfile("out"); dir.create(out)
    EzApp$new()$write_methods(output_dir = out, analysis_name = "Test")
    md <- readLines(file.path(out, "methods.md"))
    expect_equal(sum(grepl("^## ", md)), 1)          # only the analysis header
    expect_false(any(grepl("Methods$", md)))         # the model's own Methods heading is dropped
    expect_true(any(md == "#### Alignment"))         # its section headings are demoted
  })
})

## Every character constant in a citation() body (references and the strings its gates compare to).
citationLiterals <- function(e) {
  if (is.character(e)) e else if (is.call(e) || is.pairlist(e)) unlist(lapply(as.list(e), citationLiterals))
}

## The generators whose citation() takes param, i.e. offers only the references of steps that ran.
gatedCitationClasses <- function() {
  ns <- asNamespace("ezRun")
  Filter(function(cls) {
    gen <- get(cls, envir = ns)
    if (!inherits(gen, "refObjectGenerator") || cls == "EzAppSCEVANApp") return(FALSE)
    f <- tryCatch(gen$new()$citation, error = function(e) NULL)
    "param" %in% names(formals(f))
  }, ls(ns, pattern = "^EzApp"))
}

test_that("citation(list()) keeps the app's own first reference for every gated citation()", {
  classes <- gatedCitationClasses()
  expect_gt(length(classes), 25)
  for (cls in classes) {
    f <- get(cls, envir = asNamespace("ezRun"))$new()$citation
    b <- body(f)
    first <- b[[length(b)]][[2]]          # first argument of the final c(...)
    expect_type(first, "character")        # entry 1 is unconditional (write_methods always keeps it)
    expect_identical(f(list())[1], first, label = cls)
  }
})

test_that("gated citations follow the job's parameters", {
  has <- function(cits, pattern) any(grepl(pattern, cits))
  hs <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"

  de <- EzAppDeseq2$new()
  expect_false(has(de$citation(list(runGO = FALSE)), "clusterProfiler|Enrichr"))
  expect_true(has(de$citation(list(runGO = TRUE)), "clusterProfiler"))
  expect_false(has(de$citation(list(useLfcShrink = FALSE)), "ashr|new deal"))
  expect_true(has(de$citation(list(useLfcShrink = TRUE)), "new deal"))

  ed <- EzAppEdger$new()
  expect_false(has(ed$citation(list(testMethod = "glm")), "DESeq2|limma powers"))
  expect_true(has(ed$citation(list(testMethod = "deseq2")), "DESeq2"))

  li <- EzAppLimma$new()
  expect_false(has(li$citation(list(modelMethod = "limma-trend")), "voom"))
  expect_true(has(li$citation(list(modelMethod = "voom")), "voom"))

  cr <- EzAppCellRanger$new()
  expect_false(has(cr$citation(list(TenXLibrary = "GEX", runVeloCyto = FALSE)), "RNA velocity"))
  expect_true(has(cr$citation(list(TenXLibrary = "GEX", runVeloCyto = TRUE)), "RNA velocity"))
  expect_false(has(cr$citation(list(TenXLibrary = "GEX", controlSeqs = "")), "Biostrings"))
  expect_true(has(cr$citation(list(TenXLibrary = "GEX", controlSeqs = "ERCC")), "Biostrings"))

  star <- EzAppSTAR$new()
  expect_false(has(star$citation(list(barcodePattern = "", markDuplicates = FALSE)), "UMI-tools|Picard"))
  expect_true(has(star$citation(list(barcodePattern = "NNNNNNNN", markDuplicates = TRUE)), "UMI-tools"))
  expect_true(has(star$citation(list(markDuplicates = TRUE)), "Picard"))

  macs <- EzAppMacs3$new()
  expect_false(has(macs$citation(list(annotatePeaks = FALSE)), "ChIPseeker"))
  expect_true(has(macs$citation(list(annotatePeaks = TRUE)), "ChIPseeker"))
  expect_true(has(macs$citation(list(mode = "ATAC-seq")), "deepTools"))
  expect_false(has(macs$citation(list(mode = "ChIP-seq", useControl = TRUE)), "deepTools"))

  sp <- EzAppSpatialSeurat$new()
  expect_false(has(sp$citation(list(refBuild = hs, spotClean = FALSE, Azimuth = "none")), "SpotClean|Integrated analysis"))
  expect_true(has(sp$citation(list(refBuild = hs, spotClean = TRUE)), "SpotClean"))
  expect_true(has(sp$citation(list(refBuild = hs, Azimuth = "pbmcref")), "Integrated analysis"))

  xe <- EzAppXeniumSeurat$new()
  expect_false(has(xe$citation(list(rctdReference = "None", rctdFile = "", doSPLIT = TRUE)), "Robust decomposition|Bilous"))
  expect_false(has(xe$citation(list(rctdReference = "allen/ref.rds", doSPLIT = FALSE)), "Bilous"))
  expect_true(has(xe$citation(list(rctdReference = "allen/ref.rds", doSPLIT = TRUE)), "Bilous"))

  mo <- EzAppScMultiOmics$new()
  expect_false(has(mo$citation(list(runWNN = FALSE, adtNorm = "CLR")), "Integrated analysis|ADTnorm"))
  expect_true(has(mo$citation(list(runWNN = TRUE, adtNorm = "ADTnorm")), "ADTnorm"))

  mg <- EzAppMageckTest$new()
  expect_false(has(mg$citation(list(runGSEA = FALSE, species = "dre")), "Fast gene set|limma"))
  expect_true(has(mg$citation(list(runGSEA = TRUE, species = "hsa")), "limma"))
})

test_that("every citation() entry ends with exactly one URL, the anchor write_methods matches", {
  for (cls in ls(asNamespace("ezRun"), pattern = "^EzApp")) {
    gen <- get(cls, envir = asNamespace("ezRun"))
    if (!inherits(gen, "refObjectGenerator") || cls == "EzAppSCEVANApp") next  # SCEVAN cannot be instantiated (pre-existing)
    cits <- tryCatch(gen$new()$methods_citations(list()), error = function(e) character(0))
    ## gated entries are absent from citation(list()): check every reference literal in the body too
    f <- tryCatch(gen$new()$citation, error = function(e) NULL)
    if ("param" %in% names(formals(f))) cits <- union(cits, Filter(function(s) nchar(s) > 60, citationLiterals(body(f))))
    for (x in cits) {
      urls <- regmatches(x, gregexpr("https?://\\S+", x))[[1]]
      expect_length(urls, 1)
      expect_true(endsWith(x, urls), label = paste(cls, substr(x, 1, 40)))
    }
  }
})

test_that("write_methods gives facts only to a finished run of this ezRun version", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    sh <- file.path(sd, "job.sh"); writeLines("echo job", sh)
    log <- paste0(sh, "_sushiID1_x_o.log")
    here <- paste0("ezRun_", utils::packageVersion("ezRun"))
    run <- function(lines) {
      writeLines(lines, log); out <- tempfile("out"); dir.create(out)
      EzAppScSeurat$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                                        example_script = "job.sh", sample_count = 1)
      file.exists(file.path(out, "app_facts.txt"))
    }
    expect_true(run(c("other attached packages:", paste0("[1] ", here))))
    expect_false(run(c("[1] ezRun_0.0.1")))                                 # other version
    expect_false(run(c(paste0("[1] ", here), "Error in foo(): bar", "Execution halted")))  # failed
    expect_false(run("no session info"))                                   # version unknown
    other <- file.path(sd, "job2.sh_sushiID2_x_e.log")                    # another sample failed
    writeLines(c("Error: incompatible indices"), other)
    expect_false(run(c("other attached packages:", paste0("[1] ", here))))
    unlink(other)
  })
})

test_that("an empty description gives the template, and the app's own citation is always kept", {
  withStubWriter(function(argsFile) {
    stub <- file.path(dirname(argsFile), "llm_write_methods")
    writeLines(c("#!/bin/sh", "while [ $# -gt 0 ]; do [ \"$1\" = --output ] && out=$2; shift; done",
                 "printf '## Methods\\n\\n## References\\n' > \"$out\""), stub)
    out <- tempfile("out"); dir.create(out)
    expect_message(EzApp$new()$write_methods(output_dir = out, analysis_name = "T"),
                   "template fallback \\(empty description\\)")
    expect_true(any(grepl(METHODS_TEMPLATE_DECLARATION, readLines(file.path(out, "methods.md")), fixed = TRUE)))
    writeLines(c("#!/bin/sh", "while [ $# -gt 0 ]; do [ \"$1\" = --output ] && out=$2; shift; done",
                 "printf 'Counts were tested.\\n\\n## References\\n' > \"$out\""), stub)
    out <- tempfile("out"); dir.create(out)
    EzAppDeseq2$new()$write_methods(output_dir = out, analysis_name = "T")
    expect_true(any(grepl("s13059-014-0550-8", readLines(file.path(out, "methods.md")), fixed = TRUE)))
  })
})

test_that("write_methods keeps the prose when the model writes References first", {
  withStubWriter(function(argsFile) {
    stub <- file.path(dirname(argsFile), "llm_write_methods")
    writeLines(c("#!/bin/sh", "while [ $# -gt 0 ]; do [ \"$1\" = --output ] && out=$2; shift; done",
                 "printf '## References\\nLove et al. https://doi.org/10.1186/s13059-014-0550-8\\n\\nCounts were tested with DESeq2.\\n' > \"$out\""), stub)
    out <- tempfile("out"); dir.create(out)
    EzAppDeseq2$new()$write_methods(output_dir = out, analysis_name = "T")
    md <- readLines(file.path(out, "methods.md"))
    desc <- md[seq(which(md == "### Description") + 1, which(md == "### References") - 1)]
    expect_true(any(grepl("Counts were tested with DESeq2", desc)))
    expect_false(any(grepl("doi.org", desc)))
  })
})

test_that("the facts header names the ezRun version the run used", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    writeLines("[1] ezRun_3.23.1", file.path(sd, "job.sh_sushiID1_x_o.log"))
    out <- tempfile("out"); dir.create(out)
    withr::with_options(list(ezRun.methodsFactsVersion = "^ezRun_3\\.23\\."),
      EzAppScSeurat$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                                        example_script = "job.sh", sample_count = 1))
    expect_match(readLines(file.path(out, "app_facts.txt"))[1], "in ezRun 3.23.1,", fixed = TRUE)
  })
})

test_that("a failed run gets a statement, not a Methods text, and the writer is not called", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    writeLines(c("Error in EzRef(userParam): reference missing", "Execution halted"),
               file.path(sd, "job.sh_sushiID1_x_e.log"))
    out <- tempfile("out"); dir.create(out)
    EzAppScSeurat$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                                      example_script = "job.sh", sample_count = 2)
    md <- paste(readLines(file.path(out, "methods.md")), collapse = "\n")
    expect_match(md, "The analysis did not complete: 1 of 2 job(s)", fixed = TRUE)
    expect_match(md, "Error in EzRef(userParam)", fixed = TRUE)
    expect_false(file.exists(argsFile))                      # llm_write_methods never ran
  })
})

test_that("facts and citations that depend on the input dataset read it", {
  p <- list(refBuild = "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03",
            maxEmptyDropPValue = 1, tissue = "Blood", sctype.enabled = TRUE, sctype.tissue = "auto")
  h5 <- p; attr(h5, "input") <- data.frame(`CountMatrix [Link]` = "a/cellbender_filtered_seurat.h5", check.names = FALSE)
  dirIn <- p; attr(dirIn, "input") <- data.frame(`CountMatrix [Link]` = "a/filtered_feature_bc_matrix", check.names = FALSE)
  app <- EzAppScSeurat$new()
  expect_false(any(grepl("emptyDrops", app$methods_facts(h5))))
  expect_false(any(grepl("EmptyDrops", app$citation(h5))))
  expect_true(any(grepl("removed no cells", app$methods_facts(dirIn))))
  expect_true(any(grepl("EmptyDrops", app$citation(dirIn))))
  expect_true(any(grepl("emptyDrops", app$methods_facts(p))))    # input unknown: stated
  f <- app$methods_facts(p)
  expect_true(any(grepl("for the tissue Blood.", f, fixed = TRUE)))
  expect_true(any(grepl("Immune system (sctype.tissue auto)", f, fixed = TRUE)))
  q <- list(); attr(q, "input") <- data.frame(`Read Count` = c("600000000", "600000000"), check.names = FALSE)
  expect_true(any(grepl("1 billion", EzAppFastqc$new()$methods_facts(q))))
  attr(q, "input")$`Read Count` <- c("100000", "100000")
  expect_false(any(grepl("1 billion", EzAppFastqc$new()$methods_facts(q))))
})

test_that("job logs archived before the sushiID naming are found", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    writeLines("[1] ezRun_0.0.1", file.path(sd, "job.sh_20240329103602628_o.log"))
    out <- tempfile("out"); dir.create(out)
    EzApp$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                              example_script = "job.sh", sample_count = 1)
    expect_true(any(grepl("20240329103602628_o.log", readLines(argsFile), fixed = TRUE)))
  })
})

test_that("a DATASET-mode run takes its sample count from the input dataset", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    writeLines(c("Name\tRead Count", paste0("S", 1:5, "\t100")), file.path(d, "input_dataset.tsv"))
    out <- tempfile("out"); dir.create(out)
    EzApp$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                              example_script = "job.sh", sample_count = 1)
    a <- readLines(argsFile)
    expect_identical(a[which(a == "--sample-count") + 1], "5")
  })
  expect_true(any(grepl("backgroundExpression (4) is not a filter", EzAppDeseq2$new()$methods_facts(list(backgroundExpression = 4)), fixed = TRUE)))
})

test_that("a run of another ezRun version is cited from the parameters it recorded, not today's defaults", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    writeLines("[1] ezRun_3.18.1", file.path(sd, "job.sh_sushiID1_x_o.log"))
    writeLines(c("refBuild\tHomo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03", "SingleR\tnone"),
               file.path(d, "parameters.tsv"))
    out <- tempfile("out"); dir.create(out)
    EzAppScSeurat$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                                      example_script = "job.sh", sample_count = 1)
    cand <- readLines(file.path(out, "citations_candidates.txt"))
    expect_false(any(grepl("Large language model consensus", cand)))   # mLLMCelltype default is on today
    expect_true(any(grepl("Hao, Y. et al. Dictionary learning", cand, fixed = TRUE)))
  })
})
