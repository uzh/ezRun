## methods_facts() states literal values of the app code (seeds, algorithms, thresholds).
## Each row ties such a fact to the code it describes: the fact must still be emitted for
## the given parameters, and the named function(s) must still contain the value. When one
## side changes without the other, this test fails. The value is looked for only inside the
## named top-level definition(s), so another occurrence elsewhere in the file (a second
## set.seed(38), say) cannot keep a drifted anchor passing.

pkgRoot <- test_path("..", "..")
skip_if_not(dir.exists(file.path(pkgRoot, "R")), "package sources not available")

## Line numbers of the top-level definition `fn <- ...` / `fn = ...` in x, up to the line
## before the next top-level definition. An error when fn is not defined exactly once.
defWindow <- function(x, fn) {
  defs <- grep("^[A-Za-z._][A-Za-z0-9._]*\\s*(<-|=)", x)
  s <- defs[grepl(paste0("^", gsub(".", "\\.", fn, fixed = TRUE), "\\s*(<-|=)"), x[defs])]
  if (length(s) != 1) stop("no single top-level definition of ", fn)
  s:(c(defs[defs > s], length(x) + 1)[1] - 1)
}

## Code of `fn` in `file` (the whole file when fn is NA) without comments and without the
## methods_facts()/citation() bodies, so an anchor can only match code, not the fact text
## or a comment quoting the value.
codeOf <- function(file, fn = NA, x = readLines(file.path(pkgRoot, file), warn = FALSE)) {
  if (!is.na(fn)) x <- x[defWindow(x, fn)]
  drop <- integer(0)
  for (s in grep("(methods_facts|citation) = function\\(", x)) {
    depth <- 0
    for (i in s:length(x)) {
      depth <- depth + lengths(gregexpr("\\{", x[i])) * grepl("\\{", x[i]) -
        lengths(gregexpr("\\}", x[i])) * grepl("\\}", x[i])
      if (i > s && depth <= 0) break
    }
    drop <- c(drop, s:i)
  }
  if (length(drop)) x <- x[-drop]
  paste(sub("#.*$", "", x), collapse = "\n")
}

## A code regex starting with "!" must NOT match (e.g. an option the app never passes).
## code may hold several functions' code: the value must be in each of them.
anchorOk <- function(factRegex, facts, codeRegex, code) {
  neg <- startsWith(codeRegex, "!")
  inCode <- grepl(sub("^!", "", codeRegex), code, perl = TRUE)
  any(grepl(factRegex, facts, perl = TRUE)) && (if (neg) !any(inCode) else all(inCode))
}

human <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"
rctd <- list(rctdFile = "", rctdReference = "someRef", doSPLIT = TRUE, splitMode = "neighborhood",
             coocFdr = TRUE)
anchors <- list(
  ## class, param, fact regex, source file, code regex, function(s) holding the code (NA: whole file)
  list("EzAppScSeurat", list(), "seed was set to 38", "R/app-ScSeurat.R", "set\\.seed\\(38\\)", "ezMethodScSeurat"),
  list("EzAppScSeurat", list(), "seed\\.use = 38", "R/seuratUtils.R", "seed = 38\\s*\\)[\\s\\S]*seed\\.use = seed", "seuratStandardSCTPreprocessing"),
  list("EzAppScSeurat", list(), "niters = 1e5", "R/app-ScSeurat.R", "emptyDrops\\([^)]*niters = 1e5", "ezMethodScSeurat"),
  list("EzAppScSeurat", list(), "clusters = TRUE", "R/app-ScSeurat.R", "scDblFinder\\([^)]*clusters = TRUE", "addCellQcToSeurat"),
  list("EzAppScSeurat", list(), "vst\\.flavor = v2", "R/seuratUtils.R", "vst\\.flavor = \"v2\"", "seuratStandardSCTPreprocessing"),
  list("EzAppScSeurat", list(), "FindClusters algorithm 1", "R/seuratUtils.R",
       "FindClusters\\(\\s*object = scData,\\s*resolution = myResolutions,\\s*verbose = FALSE\\s*\\)", "seuratStandardWorkflow"),
  list("EzAppScSeurat", list(), "resolutions 0\\.2, 0\\.4, 0\\.6, 0\\.8 and 1\\.0", "R/seuratUtils.R",
       "seq\\(from = 0\\.2, to = 1, by = 0\\.2\\)", "seuratStandardWorkflow"),
  list("EzAppScSeurat", list(), "only\\.pos = TRUE\\); p-values were Bonferroni", "R/seuratUtils.R", "only\\.pos = TRUE", "getSeuratMarkers"),
  list("EzAppScSeurat", list(refBuild = human, computePathwayTFActivity = TRUE), "times = 100, minsize = 5",
       "R/app-ScSeurat.R", "times = 100,\\s*minsize = 5", c("computeTFActivityAnalysis", "computePathwayActivityAnalysis")),
  list("EzAppScSeurat", list(mLLMCelltype = TRUE), "temperature 0, seed 42", "R/app-ScSeurat.R",
       "temperature = 0,[^\\n]*\\n\\s*seed = 42L", "registerFgczVllmProvider"),
  list("EzAppScSeurat", list(CyteTypeR = TRUE), "p below 0\\.05 and log2 fold change above 0\\.5", "R/app-ScSeurat.R",
       "p_val_adj < 0\\.05, avg_log2FC > 0\\.5", "ezMethodScSeurat"),
  list("EzAppSpatialSeurat", list(), "seed was set to 38", "R/app-SpatialSeurat.R", "set\\.seed\\(38\\)", "ezMethodSpatialSeurat"),
  list("EzAppSpatialSeurat", list(), "r\\.metric = 5", "R/seuratUtils.R", "r\\.metric = 5", "spatialMarkers"),
  list("EzAppSpatialSeurat", list(), "below 0\\.01 \\(pvalue_allMarkers", "R/app-SpatialSeurat.R",
       "pvalue_allMarkers = ezFrame\\([^)]*DefaultValue = 0\\.01", "EzAppSpatialSeurat"),
  list("EzAppSpatialSeuratHD", list(), "seed was set to 38", "R/app-SpatialSeuratHD.R", "set\\.seed\\(38\\)", "ezMethodSpatialSeuratHD"),
  list("EzAppScSeuratCombine", list(), "seed was set to 38", "R/app-ScSeuratCombine.R", "set\\.seed\\(38\\)", "ezMethodScSeuratCombine"),
  list("EzAppScSeuratCombine", list(integrationMethod = "RPCA"), "k\\.anchor = 20", "R/seuratUtils.R", "k\\.anchor = 20", "cellClustWithCorrection"),
  list("EzAppXeniumSeurat", rctd, "seed was set to 42 immediately before RCTD", "R/app-XeniumSeurat.R", "set\\.seed\\(42\\)", "ezMethodXeniumSeurat"),
  list("EzAppXeniumSeurat", list(), "min\\.pct 0\\.25 and logfc\\.threshold 0\\.25", "R/app-XeniumSeurat.R",
       "min\\.pct = 0\\.25,\\s*logfc\\.threshold = 0\\.25", "ezMethodXeniumSeurat"),
  list("EzAppXeniumSeurat", list(), "first 30 \\(fixed in the code\\)", "R/app-XeniumSeurat.R",
       "FindNeighbors\\(sdata, dims = 1:30", "ezMethodXeniumSeurat"),
  list("EzAppXeniumSeurat", rctd, "at most 10,000 cells per cell type", "R/app-XeniumSeurat.R",
       "spacexr::Reference\\(ref_counts, ref_celltypes\\)", "ezMethodXeniumSeurat"),
  list("EzAppXeniumSeurat", rctd, "Benjamini-Hochberg-adjusted", "inst/templates/XeniumSeurat.Rmd",
       "p\\.adjust\\(P\\[ut\\], method = \"BH\"\\)", NA),
  list("EzAppXeniumSeurat", rctd, "fewer than 20 cells dropped, a random subsample of 100,000 cells when larger \\(seed 42\\)",
       "R/xeniumCooccurrence.R", "min_cells = 20,\\s*max_cells = 1e5, seed = 42", "computeCelltypeCooccurrence"),
  list("EzAppXeniumSeurat", rctd, "20 nearest spatial neighbours \\(edges over 15 um pruned\\)", "R/app-XeniumSeurat.R",
       "rad_pruning = 15, k_knn = 20", "ezMethodXeniumSeurat"),
  list("EzAppVisiumHDSeurat", list(), "seed was set to 38", "R/app-VisiumHDSeurat.R", "set\\.seed\\(38\\)", "ezMethodVisiumHDSeurat"),
  list("EzAppVisiumHDSeurat", list(), "fewer than 50,000 bins", "R/app-VisiumHDSeurat.R", "nrow\\(scData@meta\\.data\\) < 50000", "ezMethodVisiumHDSeurat"),
  list("EzAppVisiumHDSeurat", list(), "PCA computed 80 components", "R/app-VisiumHDSeurat.R", "RunPCA\\(scData, npcs = 80\\)", "ezMethodVisiumHDSeurat"),
  list("EzAppVisiumHDSeurat", list(), "k_geom = 30 \\(fixed", "R/app-VisiumHDSeurat.R", "k_geom = 30", "ezMethodVisiumHDSeurat"),
  list("EzAppVisiumHDSeurat", list(), "the first 12 \\(fixed", "R/app-VisiumHDSeurat.R", "dims = 1:12", "ezMethodVisiumHDSeurat"),
  list("EzAppVisiumHDSeurat", list(), "at resolution 0\\.5 \\(nicheResolution\\)", "R/app-VisiumHDSeurat.R",
       "nicheResolution = ezFrame\\([^)]*DefaultValue = 0\\.5", "EzAppVisiumHDSeurat"),
  list("EzAppDeseq2", list(), "minReplicatesForReplace = Inf", "R/twoGroups.R", "minReplicatesForReplace = Inf", "runDeseq2"),
  list("EzAppDeseq2", list(), "Benjamini-Hochberg adjustment of the DESeq2 Wald p-values", "R/twoGroups.R",
       "p\\.adjust\\(pValue\\[useProbe\\], method = \"fdr\"\\)", "twoGroupCountComparison"),
  list("EzAppDeseq2", list(useLfcShrink = TRUE), "lfcShrink using the ashr method", "R/twoGroups.R",
       "lfcShrink\\(dds[\\s\\S]{0,200}?type=\"ashr\"", "runDeseq2"),
  list("EzAppDeseq2", list(), "sigThresh \\(ezRun default 10\\)", "inst/extdata/EZ_PARAM_DEFAULTS.txt",
       "(?m)^sigThresh\\tnumeric\\t10\\t", NA),
  list("EzAppEdger", list(testMethod = "glm", deTest = "QL"), "dispersion was fixed at 0\\.1", "R/twoGroups.R",
       "common\\.dispersion <- 0\\.1", c("runEdger", "runGlm")),
  list("EzAppEdger", list(testMethod = "glm", deTest = "QL"), "glmQLFit, glmQLFTest", "R/twoGroups.R", "glmQLFit\\(", "runGlm"),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "150 training epochs", "R/app-cellBender.R", "!--epochs", NA),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "\\(--fpr\\) of 0\\.01", "R/app-cellBender.R", "!--fpr", NA),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "--cpu-threads \\(gpu 0\\)", "R/app-cellBender.R",
       "'--cpu-threads', param\\$cores", "ezMethodCellBender"),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 1), "--cuda \\(gpu above 0\\)", "R/app-cellBender.R",
       "param\\$gpu > 0\\) \\{\\s*cmd <- paste\\(cmd, \"--cuda\"\\)", "ezMethodCellBender"),
  list("EzAppSTAR", list(), "infer_experiment\\.py on 1,000,000 sampled reads", "R/app-mapping.R", "\"-s 1000000\"", "ezMethodSTAR"),
  list("EzAppFastqc", structure(list(), input = data.frame(`Read Count` = "2000000000", check.names = FALSE)), "ShortRead FastqSampler, seed 123", "R/fastqIO.R",
       "subsampleFastqFile <- function\\([^)]*seed = 123L\\)", "subsampleFastqFile")
)

test_that("every fact that states a literal code value still matches the code", {
  for (a in anchors) {
    cls <- a[[1]]; file <- a[[4]]; fns <- a[[6]]
    facts <- get(cls)$new()$methods_facts(a[[2]])
    expect_true(any(grepl(a[[3]], facts, perl = TRUE)), label = paste(cls, "fact:", a[[3]]))
    ## a positive anchor in R code names the function(s) it describes; NA (whole file) is
    ## for absence checks and non-R files
    neg <- startsWith(a[[5]], "!")
    if (!neg && grepl("^R/", file)) expect_false(anyNA(fns), label = paste(cls, a[[5]], "names no function"))
    code <- vapply(fns, function(f) codeOf(file, f), "")
    expect_true(anchorOk(a[[3]], facts, a[[5]], code),
                label = paste(cls, "code:", file, paste(fns, collapse = "+"), a[[5]]))
  }
  expect_true(all(c("EzAppScSeurat", "EzAppSpatialSeurat", "EzAppXeniumSeurat", "EzAppVisiumHDSeurat",
                    "EzAppDeseq2", "EzAppEdger", "EzAppCellBender") %in% vapply(anchors, `[[`, "", 1)))
})

test_that("the anchor check fails when the code value drifts (positive control)", {
  facts <- EzAppScSeurat$new()$methods_facts(list())
  code <- codeOf("R/app-ScSeurat.R")
  expect_true(anchorOk("seed was set to 38", facts, "set\\.seed\\(38\\)", code))
  expect_false(anchorOk("seed was set to 38", facts, "set\\.seed\\(38\\)", gsub("set.seed(38)", "set.seed(39)", code, fixed = TRUE)))
  cb <- codeOf("R/app-cellBender.R")
  cbFacts <- EzAppCellBender$new()$methods_facts(list(cmdOptions = "", gpu = 0))
  expect_true(anchorOk("150 training epochs", cbFacts, "!--epochs", cb))
  expect_false(anchorOk("150 training epochs", cbFacts, "!--epochs", paste(cb, "cmd <- paste(cmd, '--epochs 300')")))
  ## the facts text itself is not code: stripping it hides "niters = 1e5" in the fact string
  expect_false(grepl("DropletUtils emptyDrops (niters", code, fixed = TRUE))
})

test_that("an anchor checks the function the fact describes, not the whole file (positive control)", {
  x <- readLines(file.path(pkgRoot, "R/app-ScSeurat.R"), warn = FALSE)
  facts <- EzAppScSeurat$new()$methods_facts(list())
  ## change only the set.seed(38) of ezMethodScSeurat; the one in addCellQcToSeurat stays
  win <- defWindow(x, "ezMethodScSeurat")
  mut <- x; mut[win] <- gsub("set.seed(38)", "set.seed(39)", mut[win], fixed = TRUE)
  expect_false(identical(mut, x))
  expect_true(grepl("set.seed(38)", codeOf("R/app-ScSeurat.R", x = mut), fixed = TRUE))
  ## the whole-file check (the old one) still passes on the mutant; the scoped one fails
  expect_true(anchorOk("seed was set to 38", facts, "set\\.seed\\(38\\)", codeOf("R/app-ScSeurat.R", x = mut)))
  expect_false(anchorOk("seed was set to 38", facts, "set\\.seed\\(38\\)", codeOf("R/app-ScSeurat.R", "ezMethodScSeurat", x = mut)))
  expect_true(anchorOk("seed was set to 38", facts, "set\\.seed\\(38\\)", codeOf("R/app-ScSeurat.R", "ezMethodScSeurat", x = x)))
  ## a named function that does not exist is an error, not a silent pass
  expect_error(codeOf("R/app-ScSeurat.R", "noSuchFunction"))
})

test_that("every off-step rule names a parameter the app code or its report reads", {
  files <- c(file.path("R", list.files(file.path(pkgRoot, "R"), pattern = "\\.R$")),
             file.path("inst/templates", list.files(file.path(pkgRoot, "inst/templates"), pattern = "\\.Rmd$")))
  code <- paste(vapply(files, codeOf, ""), collapse = "\n")
  for (cls in names(METHODS_OFFSTEP_RULES)) {
    expect_true(exists(cls), label = cls)
    for (p in unlist(strsplit(names(METHODS_OFFSTEP_RULES[[cls]]), "+", fixed = TRUE))) {
      q <- gsub(".", "\\.", p, fixed = TRUE)
      read <- paste0("param\\$", q, "(?![A-Za-z0-9._])|param\\[\\[[\"']", q, "[\"']\\]\\]")
      expect_true(grepl(read, code, perl = TRUE), label = paste(cls, p))
    }
  }
})
