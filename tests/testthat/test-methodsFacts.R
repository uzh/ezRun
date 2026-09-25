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
                 "printf '## Methods\\n\\nReads were aligned.\\n' > \"$out\""), stub)
    out <- tempfile("out"); dir.create(out)
    EzApp$new()$write_methods(output_dir = out, analysis_name = "Test")
    md <- readLines(file.path(out, "methods.md"))
    expect_equal(sum(grepl("^## ", md)), 1)          # only the analysis header
    expect_true(any(md == "#### Methods"))
  })
})

test_that("every citation() entry ends with exactly one URL, the anchor write_methods matches", {
  for (cls in ls(asNamespace("ezRun"), pattern = "^EzApp")) {
    gen <- get(cls, envir = asNamespace("ezRun"))
    if (!inherits(gen, "refObjectGenerator") || cls == "EzAppSCEVANApp") next  # SCEVAN cannot be instantiated (pre-existing)
    cits <- tryCatch(gen$new()$methods_citations(list()), error = function(e) character(0))
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
