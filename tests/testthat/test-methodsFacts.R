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
    expect_match(readLines(factsFile)[1], "^Fixed behaviour of this app's code")   # the header, then the facts

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

test_that("no fact names a step that the parameter SUSHI sets to off switched off", {
  hs <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"
  ## the off value SUSHI's form offers; otherwise FALSE for a logical default and "" for the rest
  sushiOff <- list(SingleR = "none", Azimuth = "none", integrationMethod = "none", rctdReference = "None")
  hits <- character(0)
  for (cls in names(METHODS_OFFSTEP_RULES)) {
    app <- get(cls)$new()
    def <- suppressMessages(suppressWarnings(ezParam(list(), appDefaults = app$appDefaults)))
    for (p in unique(unlist(strsplit(names(METHODS_OFFSTEP_RULES[[cls]]), "+", fixed = TRUE)))) {
      q <- list(refBuild = hs)
      q[[p]] <- if (!is.null(sushiOff[[p]])) sushiOff[[p]] else if (is.logical(def[[p]])) FALSE else ""
      h <- methods_check_offsteps(suppressWarnings(app$methods_facts(q)), cls, q)
      if (length(h)) hits <- c(hits, paste(cls, h))
    }
  }
  expect_identical(hits, character(0))
  ## positive control: the check does flag a fact that names the switched-off step
  expect_length(methods_check_offsteps("STAR ran in two-pass mode.", "EzAppSTAR", list(twopassMode = FALSE)), 1)
  expect_true(any(grepl("--twopassMode None", EzAppSTAR$new()$methods_facts(list(twopassMode = FALSE)), fixed = TRUE)))
})

test_that("ScSeuratCompare names sccomp only when replicateGrouping is set", {
  app <- EzAppScSeuratCompare$new()
  off <- app$methods_facts(list(replicateGrouping = ""))
  expect_false(any(grepl("sccomp", off)))
  expect_true(any(grepl("seed was set to 38", off)))
  expect_true(any(grepl("sccomp's sampling seed was drawn from", app$methods_facts(list(replicateGrouping = "Replicate")))))
})

test_that("ScMultiOmics states BD Rhapsody and ATAC facts only where they apply", {
  app <- EzAppScMultiOmics$new()
  p <- list(runWNN = TRUE, adtNorm = "ADTnorm", vdjChain = "auto")
  bd <- p; attr(bd, "input") <- data.frame(`SCDataOrigin [Factor]` = "BDRhapsody", check.names = FALSE)
  tenx <- p; attr(tenx, "input") <- data.frame(`CountMatrix [Link]` = "a/filtered_feature_bc_matrix", check.names = FALSE)
  for (q in list(p, tenx)) {
    f <- app$methods_facts(q)
    expect_false(any(grepl("BD Rhapsody input", f)))
    expect_true(any(grepl("upstream ScSeurat object", f)))
  }
  f <- app$methods_facts(bd)
  expect_true(any(grepl("BD Rhapsody input", f)))
  expect_false(any(grepl("upstream ScSeurat object|Antibody Capture features of the count-matrix H5|normalized with ADTnorm", f)))
  ## WNN: LSI components only when ATAC was present
  wnn <- grep("FindMultiModalNeighbors", app$methods_facts(tenx), value = TRUE)
  expect_length(wnn, 1)
  expect_match(wnn, "LSI components 2-30 only when ATAC was present", fixed = TRUE)
})

test_that("VeloCyto and FastqScreen state input-dependent steps from the input dataset", {
  vc <- EzAppVeloCyto$new()
  bd <- list(); attr(bd, "input") <- data.frame(`SCDataOrigin [Factor]` = "BDRhapsody", check.names = FALSE)
  tenx <- list(); attr(tenx, "input") <- data.frame(ResultDir = "a/cellranger", check.names = FALSE)
  expect_false(any(grepl("run10x on the CellRanger output|CellRanger Multi output|CRAM alignments", vc$methods_facts(bd))))
  expect_true(any(grepl("XF:Z:__intergenic", vc$methods_facts(bd))))
  expect_false(any(grepl("BD Rhapsody data", vc$methods_facts(tenx))))
  expect_true(any(grepl("run10x on the CellRanger output", vc$methods_facts(tenx))))
  expect_true(all(c(any(grepl("run10x", vc$methods_facts())), any(grepl("XF:Z", vc$methods_facts())))))   # input unknown

  fs <- EzAppFastqScreen$new()
  virus <- function(q) grep("pathogenic viruses", fs$methods_facts(q), value = TRUE)
  mouse <- list(virusCheck = FALSE); attr(mouse, "input") <- data.frame(`Species` = "Mus musculus (mouse)", check.names = FALSE)
  human <- list(virusCheck = FALSE); attr(human, "input") <- data.frame(`Species` = "Homo sapiens (human)", check.names = FALSE)
  expect_length(virus(mouse), 0)
  expect_length(virus(human), 1)
  expect_no_match(virus(human), "When virusCheck")
  expect_length(virus(list(virusCheck = TRUE)), 1)
  expect_match(virus(list(virusCheck = FALSE)), "always when the dataset's Species is human")   # input unknown
})

test_that("DE facts state GO enrichment only as the code runs it", {
  for (cls in c("EzAppDeseq2", "EzAppEdger", "EzAppLimma")) {
    app <- get(cls)$new()
    f <- app$methods_facts(list(runGO = TRUE, rankMetric = "pValue", testMethod = "glm", deTest = "QL"))
    gsea <- grep("GSEA", f, value = TRUE)
    expect_length(gsea, 1)
    expect_match(gsea, "^When GO annotation was available for the reference", label = cls)
    expect_match(gsea, "ranked by -log10 p-value (rankMetric pValue)", fixed = TRUE, label = cls)
    ## ezGSEA reads param$rankMetric; without it the run stops, so no GO fact
    expect_false(any(grepl("GSEA|over-representation", app$methods_facts(list(runGO = TRUE)))), label = cls)
    ## the backgroundExpression sentence points at the sigThresh sentence where it is
    bg <- grep("presence rule", f); pres <- grep("called present in a sample", f)
    expect_true(grepl("below", f[bg]) && bg < pres || grepl("above", f[bg]) && bg > pres, label = cls)
  }
  ## Limma declares no rankMetric default: a Limma run with runGO gets no GO facts
  typed <- function(app) suppressMessages(suppressWarnings(ezParam(list(runGO = "true"), appDefaults = app$appDefaults)))
  lim <- EzAppLimma$new()
  expect_false(any(grepl("GSEA|over-representation", lim$methods_facts(typed(lim)))))
  de <- EzAppDeseq2$new()          # control: its appDefaults set rankMetric log2Ratio
  expect_true(any(grepl("ranked by log2 ratio (rankMetric log2Ratio)", de$methods_facts(typed(de)), fixed = TRUE)))
})

test_that("ScSeurat says which clustering is reported for a resolution with two decimals", {
  app <- EzAppScSeurat$new()
  expect_true(any(grepl("lowest-resolution clustering was reported instead", app$methods_facts(list(resolution = 0.25)))))
  expect_true(any(grepl("lowest-resolution clustering was reported instead", app$methods_facts(list(resolution = "[0.25, 0.2, 0.4]")))))
  expect_false(any(grepl("lowest-resolution", app$methods_facts(list(resolution = 0.6)))))
  expect_true(any(grepl("plus the resolution parameter rounded to one decimal", app$methods_facts(list()))))
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

## Every reference an app can offer, by class: citation(list()) plus, for a gated
## citation(param), every long literal of its body (gated entries are absent from list()).
allCitationEntries <- function() {
  out <- list()
  for (cls in ls(asNamespace("ezRun"), pattern = "^EzApp")) {
    gen <- get(cls, envir = asNamespace("ezRun"))
    if (!inherits(gen, "refObjectGenerator") || cls == "EzAppSCEVANApp") next  # SCEVAN cannot be instantiated (pre-existing)
    cits <- tryCatch(gen$new()$methods_citations(list()), error = function(e) character(0))
    f <- tryCatch(gen$new()$citation, error = function(e) NULL)
    if ("param" %in% names(formals(f))) cits <- union(cits, Filter(function(s) nchar(s) > 60, citationLiterals(body(f))))
    if (length(cits)) out[[cls]] <- cits
  }
  out
}

test_that("no citation() entry carries an editor's note", {
  entries <- allCitationEntries()
  expect_gt(length(entries), 40)
  for (cls in names(entries)) for (x in entries[[cls]]) {
    ## "[preprint, not peer-reviewed]" is the one bracket meant for the customer
    y <- gsub("[preprint, not peer-reviewed]", "", x, fixed = TRUE)
    expect_false(grepl("[", y, fixed = TRUE) || grepl("could not be|verified", y, ignore.case = TRUE),
                 label = paste(cls, substr(x, 1, 60)))
  }
})

test_that("gated citations offer no reference for a step the app code skips", {
  has <- function(cits, pattern) any(grepl(pattern, cits))
  hs <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"
  mm <- "Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-03"
  dr <- "Danio_rerio/Ensembl/GRCz11/Annotation/Release_110-2023-10-30"
  sc <- "Saccharomyces_cerevisiae/Ensembl/R64/Annotation/Release_110-2023-10-30"
  ## ScSeuratCompare.Rmd: org.Hs/KEGG hsa for Human, org.Mm/mmu for Mouse, no enrichment otherwise
  cmp <- EzAppScSeuratCompare$new()
  expect_false(has(cmp$citation(list(refBuild = hs)), "org.Mm.eg.db"))
  expect_true(has(cmp$citation(list(refBuild = hs)), "org.Hs.eg.db"))
  expect_false(has(cmp$citation(list(refBuild = mm)), "org.Hs.eg.db"))
  expect_false(has(cmp$citation(list(refBuild = dr)), "clusterProfiler|KEGG|org\\.(Hs|Mm)"))
  expect_true(has(cmp$citation(list()), "clusterProfiler"))   # no refBuild: the app infers Human or Mouse
  ## twoGroups.Rmd: Enrichr results only with runGO, doEnrichr (organism, featureLevel gene) and doPrecomputeEnrichr
  on <- list(refBuild = hs, featureLevel = "gene", runGO = TRUE, doPrecomputeEnrichr = TRUE)
  for (app in list(EzAppDeseq2$new(), EzAppEdger$new(), EzAppLimma$new())) {
    expect_true(has(app$citation(on), "Enrichr"))
    expect_false(has(app$citation(modifyList(on, list(doPrecomputeEnrichr = FALSE))), "Enrichr"))
    expect_false(has(app$citation(modifyList(on, list(refBuild = sc))), "Enrichr"))
    expect_false(has(app$citation(modifyList(on, list(featureLevel = "isoform"))), "Enrichr"))
    expect_true(has(app$citation(modifyList(on, list(refBuild = sc))), "clusterProfiler"))
  }
  ## ScMultiOmics: Signac/ensembldb only when the CountMatrix has ATAC files next to it (findATACFiles)
  rna <- file.path(tempfile("mo"), "filtered_feature_bc_matrix"); dir.create(rna, recursive = TRUE)
  atac <- file.path(tempfile("mo"), "filtered_feature_bc_matrix"); dir.create(atac, recursive = TRUE)
  file.create(file.path(dirname(atac), c("atac_fragments.tsv.gz", "atac_peaks.bed")))
  mo <- EzAppScMultiOmics$new()
  pr <- list(refBuild = hs, dataRoot = ""); attr(pr, "input") <- data.frame(`CountMatrix [Link]` = rna, check.names = FALSE)
  pa <- pr; attr(pa, "input") <- data.frame(`CountMatrix [Link]` = atac, check.names = FALSE)
  expect_false(has(mo$citation(pr), "Signac|ensembldb"))
  expect_true(has(mo$citation(pa), "Signac"))
  expect_true(has(mo$citation(pa), "ensembldb"))
  expect_true(has(mo$citation(list(refBuild = hs)), "Signac"))  # input unknown: kept
})

## write_methods always keeps entry 1, so it must be a tool that runs whatever the parameters.
test_that("the always-kept first citation is a tool that ran, even with the optional step off", {
  dna <- EzAppDnaBamStats$new()$methods_citations(list(runQualimap = FALSE))
  expect_match(dna[1], "Rsamtools")                       # getBamMultiMatching, every sample
  expect_false(any(grepl("Qualimap", dna)))
  expect_true(any(grepl("Qualimap", EzAppDnaBamStats$new()$methods_citations(list(runQualimap = TRUE)))))
  fl <- EzAppFlash$new()$methods_citations(list(skipFlash = TRUE))
  expect_match(fl[1], "fastp")                            # ezMethodFastpTrim, every sample
  expect_false(any(grepl("FLASH", fl)))
  expect_true(any(grepl("FLASH", EzAppFlash$new()$methods_citations(list(skipFlash = FALSE)))))
})

test_that("citations name the registered author and the paper of the step that ran", {
  has <- function(cits, pattern) any(grepl(pattern, cits))
  hs <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"
  ## package DOIs: first author as registered (DataCite for 10.18129, Crossref for 10.32614), 2026-09-28
  creator <- c("10.18129/B9.bioc.Biostrings" = "Pagès", "10.18129/B9.bioc.celldex" = "Aran",
               "10.18129/B9.bioc.GO.db" = "Carlson", "10.18129/B9.bioc.org.Hs.eg.db" = "Carlson",
               "10.18129/B9.bioc.org.Mm.eg.db" = "Carlson", "10.18129/B9.bioc.rhdf5" = "Fischer",
               "10.18129/B9.bioc.Rsamtools" = "Morgan", "10.18129/B9.bioc.seqLogo" = "Bembom",
               "10.32614/CRAN.package.metap" = "Dewey")
  entries <- unlist(allCitationEntries(), use.names = FALSE)
  pkg <- grep("doi.org/10\\.(18129|32614)/", entries, value = TRUE)
  expect_gt(length(pkg), 8)
  for (x in pkg) {
    doi <- sub(".*doi.org/", "", x)
    expect_true(doi %in% names(creator), label = doi)
    expect_true(startsWith(x, creator[doi]), label = substr(x, 1, 50))
  }
  ## CellBender reads/writes 10x files with DropletUtils but never runs emptyDrops
  expect_false(has(EzAppCellBender$new()$methods_citations(list()), "EmptyDrops"))
  ## the nanopore paper only for ONT input
  hi <- EzAppHifiasm$new()
  expect_false(has(hi$methods_citations(list(inputType = "HiFi")), "nanopore"))
  expect_true(has(hi$methods_citations(list(inputType = "ONT")), "nanopore"))
  ## ScSeuratCombine: the paper of the integration that ran
  co <- EzAppScSeuratCombine$new()
  expect_false(has(co$methods_citations(list(integrationMethod = "Harmony")), "Comprehensive Integration"))
  expect_true(has(co$methods_citations(list(integrationMethod = "Harmony")), "Harmony"))
  expect_false(has(co$methods_citations(list(integrationMethod = "CCA")), "Harmony"))
  expect_true(has(co$methods_citations(list(integrationMethod = "RPCA")), "Comprehensive Integration"))
  expect_false(has(co$methods_citations(list(integrationMethod = "none")), "Harmony|Comprehensive Integration"))
  expect_false(has(co$methods_citations(list(refBuild = hs, computePathwayTFActivity = FALSE, enrichrDatabase = "")),
                   "decoupleR|Enrichr"))
  expect_true(has(co$methods_citations(list(refBuild = hs, computePathwayTFActivity = TRUE, enrichrDatabase = "x")),
                  "decoupleR"))
})

## Software with no paper and no registered DOI (checked 2026-09-28 on Crossref,
## DataCite, Zenodo and the CRAN DOI prefix): the entry ends in its web page instead.
citationUrlAllowList <- c(
  "https://www.10xgenomics.com/support/software/cell-ranger",
  "https://www.10xgenomics.com/support/software/cell-ranger-arc",
  "https://www.10xgenomics.com/support/software/space-ranger",
  "https://broadinstitute.github.io/picard/",
  "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
  "https://github.com/lh3/seqtk",
  "https://github.com/uzh/ezRun",
  "https://github.com/10XGenomics/loupeR",
  "https://github.com/satijalab/seurat-wrappers",
  "https://github.com/p-gueguen/rctd-py"
)

test_that("every citation() entry ends with exactly one URL, the anchor write_methods matches", {
  entries <- allCitationEntries()
  for (cls in names(entries)) {
    cits <- entries[[cls]]
    for (x in cits) {
      urls <- regmatches(x, gregexpr("https?://\\S+", x))[[1]]
      expect_length(urls, 1)
      expect_true(endsWith(x, urls), label = paste(cls, substr(x, 1, 40)))
      ## a DOI where the tool has one; a web page only for allow-listed software
      expect_true(grepl("^https://doi\\.org/10\\.[0-9]+/\\S+$", urls[1]) || urls[1] %in% citationUrlAllowList,
                  label = paste(cls, urls[1]))
    }
  }
})

test_that("write_methods gives facts only to a finished run of this ezRun version", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    sh <- file.path(sd, "job.sh"); writeLines("echo job", sh)
    log <- paste0(sh, "_sushiID1_x_o.log")
    here <- paste0("ezRun_", utils::packageVersion("ezRun"))
    ## R writes "Execution halted" to stderr, i.e. the job's _e.log (as in the archived runs)
    run <- function(lines, err = "done") {
      writeLines(lines, log); writeLines(err, sub("_o\\.log$", "_e.log", log))
      out <- tempfile("out"); dir.create(out)
      EzAppScSeurat$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                                        example_script = "job.sh", sample_count = 1)
      file.exists(file.path(out, "app_facts.txt"))
    }
    expect_true(run(c("other attached packages:", paste0("[1] ", here))))
    expect_false(run(c("[1] ezRun_0.0.1")))                                 # other version
    expect_false(run(paste0("[1] ", here), err = c("Error in foo(): bar", "Execution halted")))  # failed
    expect_false(run("no session info"))                                   # version unknown
    other <- file.path(sd, "job2.sh_sushiID2_x_e.log")                    # another sample failed
    writeLines(c("Error: incompatible indices", "Execution halted"), other)
    expect_false(run(c("other attached packages:", paste0("[1] ", here))))
    unlink(other)
  })
})

test_that("facts are withheld when any job's latest attempt failed, and kept after a successful rerun", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh")); writeLines("echo job", file.path(sd, "job2.sh"))
    writeLines(c("other attached packages:", paste0("[1] ezRun_", utils::packageVersion("ezRun"))),
               file.path(sd, "job.sh_sushiID1_2025-10-14--16-50-51_o.log"))
    writeLines("done", file.path(sd, "job.sh_sushiID1_2025-10-14--16-50-51_e.log"))
    other <- file.path(sd, "job2.sh_sushiID2_2025-10-14--14-41-50_e.log")
    facts <- function() {
      out <- tempfile("out"); dir.create(out)
      suppressMessages(EzAppScSeurat$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                                                         example_script = "job.sh", sample_count = 2))
      file.exists(file.path(out, "app_facts.txt"))
    }
    writeLines("done", other)
    expect_true(facts())                                                  # positive control
    writeLines(c("Error: incompatible indices", "Execution halted"), other)  # another job failed
    expect_false(facts())
    writeLines("done", other)                                            # the example job was rerun
    writeLines(c("Error: boom", "Execution halted"), file.path(sd, "job.sh_sushiID1_2025-10-14--14-41-53_e.log"))
    expect_true(facts())
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

test_that("the facts header says what the facts are without an ezRun version the writer would repeat", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    writeLines("[1] ezRun_3.23.1", file.path(sd, "job.sh_sushiID1_x_o.log"))
    out <- tempfile("out"); dir.create(out)
    withr::with_options(list(ezRun.methodsFactsVersion = "^ezRun_3\\.23\\."),
      EzAppScSeurat$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                                        example_script = "job.sh", sample_count = 1))
    head <- readLines(file.path(out, "app_facts.txt"))[1]
    expect_identical(head, paste("Fixed behaviour of this app's code, read from its source and the run's parameters.",
                                 "It may name options and code values; describe them in words."))
    expect_no_match(head, "ezRun|3\\.23")
  })
})

test_that("a failed run gets a statement, not a Methods text, and the writer is not called", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    writeLines(c("Error in EzRef(userParam) :", "  reference missing in /srv/GT/reference/Genes/genes.gtf", "Execution halted"),
               file.path(sd, "job.sh_sushiID1_x_e.log"))
    out <- tempfile("out"); dir.create(out)
    EzAppScSeurat$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                                      example_script = "job.sh", sample_count = 2)
    md <- paste(readLines(file.path(out, "methods.md")), collapse = "\n")
    expect_match(md, "The analysis did not complete: its job stopped with an error", fixed = TRUE)
    expect_match(md, "Error in EzRef(userParam) : reference missing in genes.gtf", fixed = TRUE)
    expect_no_match(md, "/srv/GT")
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
    writeLines("samples\tS1,S3", file.path(d, "parameters.tsv"))      # a run on 2 of the 5 rows
    EzApp$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                              example_script = "job.sh", sample_count = 1)
    a <- readLines(argsFile)
    expect_identical(a[which(a == "--sample-count") + 1], "2")
    EzApp$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                              example_script = "job.sh", sample_count = 3)   # SAMPLE mode: SUSHI's count
    a <- readLines(argsFile)
    expect_identical(a[which(a == "--sample-count") + 1], "3")
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

test_that("a caught error is not a failure; one halted job of several gives a note", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    writeLines(c("Error in autoEstCont(sc): caught", "done"), file.path(sd, "job.sh_sushiID1_x_e.log"))
    writeLines(c("Error in plot_layout(): boom", "Execution halted"), file.path(sd, "job2.sh_sushiID2_x_e.log"))
    writeLines("done", file.path(sd, "job3.sh_sushiID3_x_e.log"))
    writeLines("done", file.path(sd, "methods_dataset_1.sh_sushiID9_x_e.log"))   # the Methods job itself: not a job of the run
    out <- tempfile("out"); dir.create(out)
    EzApp$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                              example_script = "job.sh", sample_count = 3)
    md <- paste(readLines(file.path(out, "methods.md")), collapse = "\n")
    expect_no_match(md, "did not complete")
    expect_match(md, "Note: 1 of 3 jobs of this run stopped with an error (first: \"Error in plot_layout(): boom\")", fixed = TRUE)
  })
})

test_that("a job rerun after a halted attempt counts once, by its latest attempt", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    halted <- c("Error in plot_layout(): boom", "Execution halted")
    early <- file.path(sd, "job.sh_sushiID1_2025-10-14--14-41-53_e.log")
    late <- file.path(sd, "job.sh_sushiID1_2025-10-14--16-50-51_e.log")   # resubmitted: same sushiID
    md <- function() {
      out <- tempfile("out"); dir.create(out)
      EzApp$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                                example_script = "job.sh", sample_count = 1)
      paste(readLines(file.path(out, "methods.md")), collapse = "\n")
    }
    writeLines(halted, early); writeLines("done", late)                  # the rerun succeeded
    x <- md()
    expect_no_match(x, "did not complete")
    expect_no_match(x, "Note:")
    writeLines("done", early); writeLines(halted, late)                  # the rerun halted
    expect_match(md(), "The analysis did not complete: its job stopped", fixed = TRUE)
    writeLines("echo job", file.path(sd, "job2.sh"))                     # one of two jobs failed
    writeLines(c("done"), file.path(sd, "job2.sh_sushiID2_2025-10-14--14-41-50_e.log"))
    expect_match(md(), "Note: 1 of 2 jobs of this run stopped", fixed = TRUE)
  })
  ## logs archived before the sushiID naming group the same way
  expect_identical(methodsLatestJobLogs(c("s/x.sh_20220615164150845_e.log", "s/x.sh_20220614164150845_e.log",
                                          "s/y.sh_sushiID3_2026-01-01--00-00-00_e.log")),
                   c("s/x.sh_20220615164150845_e.log", "s/y.sh_sushiID3_2026-01-01--00-00-00_e.log"))
})

test_that("the writer reads only the latest attempt of the example job", {
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    for (st in c("2025-10-14--14-41-53", "2025-10-14--16-50-51")) for (k in c("o", "e"))
      writeLines(if (st < "2025-10-14--16") c("Error in f(): boom", "Execution halted") else "done",
                 file.path(sd, sprintf("job.sh_sushiID1_%s_%s.log", st, k)))
    out <- tempfile("out"); dir.create(out)
    EzApp$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                              example_script = "job.sh", sample_count = 1)
    a <- readLines(argsFile)
    expect_true(any(grepl("16-50-51_e.log", a, fixed = TRUE)))
    expect_true(any(grepl("16-50-51_o.log", a, fixed = TRUE)))
    expect_false(any(grepl("14-41-53", a, fixed = TRUE)))   # the superseded, halted attempt
  })
})

test_that("a job SLURM killed, or whose traced script stopped before the copy to gStore, failed", {
  ## fixture logs shaped like the archived ones: SUSHI's `set -eux` trace starts with
  ## "+ umask 0002" and a finished job's _e.log ends with the footer's "+ rm -rf <scratch>"
  traced <- c("+ set -o pipefail", "+ umask 0002", "+ SCRATCH_DIR=/scratch/r_S1_temp1", "+ R --vanilla --slave",
              "INFO (ezRun) [2026-01-01 00:00:00] Starting\tEzAppX")
  done <- c(traced, "+ g-req -w copy S1 /srv/gstore/projects/p1/r", "OK", "+ cd /scratch", "+ rm -rf /scratch/r_S1_temp1", "")
  timeLimit <- c(traced, "slurmstepd: error: *** JOB 123 ON fgcz-c-1 CANCELLED AT 2026-01-01T02:00:00 DUE TO TIME LIMIT ***")
  oom <- c("Loading ezRun", "Killed", paste("slurmstepd-fgcz-h-1: error: Detected 8 oom-kill event(s) in step 898.batch cgroup.",
                                            "Some of your processes may have been killed by the cgroup out-of-memory handler."))
  cancelled <- c(traced, "slurmstepd-fgcz-h-151: error: *** JOB 220315 ON fgcz-h-151 CANCELLED AT 2026-07-21T13:53:04 ***")
  ## the analysis finished and only the trailing g-req copy failed: SUSHI marks it FAILED, but the
  ## outputs exist (routine for Cell Ranger jobs), so it is written up, not reported as not completed
  copyFailed <- c(traced, "+ g-req -w copy S1 /srv/gstore/projects/p1/r", "Error while processing request: Destination path already exists!")
  cut <- c(traced, "+ samtools index S1.bam", "bash: line 12: samtools: command not found")   # stopped before the copy
  old <- c("Loading ezRun", "done")                                                      # before the trace: no marker
  f <- function(lines) { x <- tempfile(fileext = "_e.log"); writeLines(lines, x); methodsJobFailed(x) }
  expect_false(f(done))
  expect_false(f(old))
  expect_true(f(timeLimit))
  expect_true(f(oom))
  expect_true(f(cancelled))
  expect_true(f(cut))
  expect_false(f(copyFailed))
  expect_true(f(c(done[-length(done)], "Execution halted")))
  withStubWriter(function(argsFile) {
    d <- tempfile("res"); sd <- file.path(d, "scripts"); dir.create(sd, recursive = TRUE)
    writeLines("echo job", file.path(sd, "job.sh"))
    writeLines(c(traced, "Error in autoEstCont(sc): caught", timeLimit[length(timeLimit)]),
               file.path(sd, "job.sh_sushiID1_2026-01-01--00-00-00_e.log"))
    out <- tempfile("out"); dir.create(out)
    EzApp$new()$write_methods(gstore_script_dir = sd, output_dir = out, analysis_name = "T",
                              example_script = "job.sh", sample_count = 1)
    md <- paste(readLines(file.path(out, "methods.md")), collapse = "\n")
    expect_match(md, "The analysis did not complete: its job stopped with an error (first: \"SLURM cancelled the job due to time limit\")", fixed = TRUE)
    expect_false(file.exists(argsFile))
  })
})

test_that("FastQC cites ShortRead only when it subsampled (or when the input is unknown)", {
  q <- list(); attr(q, "input") <- data.frame(`Read Count` = "100000", check.names = FALSE)
  expect_false(any(grepl("ShortRead", EzAppFastqc$new()$citation(q))))
  attr(q, "input")$`Read Count` <- "2000000000"
  expect_true(any(grepl("ShortRead", EzAppFastqc$new()$citation(q))))
  expect_true(any(grepl("ShortRead", EzAppFastqc$new()$citation(list()))))
  ## a known input without a Read Count column: sum(NULL) = 0, ezMethodFastQC never subsamples
  attr(q, "input") <- data.frame(Name = "S1", `Read1 [File]` = "S1.fastq.gz", check.names = FALSE)
  expect_false(any(grepl("ShortRead", EzAppFastqc$new()$citation(q))))
})

test_that("round-5 gaps: marker pre-filters, mLLMCelltype tissue, MAD fields, gene summing, cell calling", {
  sc <- EzAppScSeurat$new()
  has <- function(fs, x) any(grepl(x, fs, fixed = TRUE))
  expect_true(has(sc$methods_facts(list()), "at least a fraction 0.1 of either group (min.pct) with an average log2 fold change of at least 0.25"))
  expect_true(has(sc$methods_facts(list(min.pct = 0.2, logfc.threshold = 0.5)), "fraction 0.2 of either group (min.pct) with an average log2 fold change of at least 0.5"))
  expect_true(has(sc$methods_facts(list(mLLMCelltype = TRUE, mLLMCelltype.tissue = "PBMC")), "with the tissue context 'PBMC'."))
  expect_true(has(sc$methods_facts(list(mLLMCelltype = TRUE, mLLMCelltype.tissue = "auto", sctype.tissue = "Immune system")),
                  "with the tissue context 'Immune system' (sctype.tissue"))
  expect_true(has(sc$methods_facts(list(mLLMCelltype = TRUE, mLLMCelltype.tissue = "auto")), "without a tissue context"))
  ## p40607a: perc_riboprot fixed, the other three empty -> only those three by MAD
  fs <- sc$methods_facts(list(nUMI = "", ngenes = "", perc_mito = "", perc_riboprot = 70, nmad = 3))
  expect_true(has(fs, "Cells were removed at fixed thresholds: perc_riboprot = 70"))
  expect_true(has(fs, "For nUMI, ngenes, perc_mito, cells more than 3 median absolute deviations"))
  for (app in list(EzAppDeseq2$new(), EzAppEdger$new()))
    expect_true(has(app$methods_facts(list()), "counts given per transcript, as from kallisto, were first summed per gene over all of its transcripts; gene-level counts, as from featureCounts, were used as they are"))
  expect_true(has(EzAppDeseq2$new()$methods_facts(list()), "Only the samples of the sample group and the reference group were used in the DESeq2 fit"))
  cr <- EzAppCellRanger$new()
  expect_true(has(cr$methods_facts(list(TenXLibrary = "GEX", expectedCells = "")), "so Cell Ranger called cells automatically"))
  expect_false(has(cr$methods_facts(list(TenXLibrary = "GEX", expectedCells = 5000)), "called cells automatically"))
  crm <- EzAppCellRangerMulti$new()
  fs <- crm$methods_facts(list(TenXLibrary = "GEX,VDJ-T", expectedCells = "", chemistry = "auto"))
  expect_true(has(fs, "no expect-cells in config.csv"))
  expect_true(has(fs, "Chemistry was left on auto, so Cell Ranger detected it"))
  expect_false(has(crm$methods_facts(list(TenXLibrary = "GEX", expectedCells = 3000, chemistry = "SC5P-R2")), "left on auto"))
})
