## methods_facts() states literal values of the app code (seeds, algorithms, thresholds).
## Each row ties such a fact to the code line it describes: the fact must still be emitted
## for the given parameters, and the code must still contain the value. When one side
## changes without the other, this test fails.

pkgRoot <- test_path("..", "..")
skip_if_not(dir.exists(file.path(pkgRoot, "R")), "package sources not available")

## Source code without comments and without the methods_facts()/citation() bodies, so an
## anchor can only match code, not the fact text or a comment quoting the value.
codeOf <- function(file) {
  x <- readLines(file.path(pkgRoot, file), warn = FALSE)
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
anchorOk <- function(factRegex, facts, codeRegex, code) {
  neg <- startsWith(codeRegex, "!")
  inCode <- grepl(sub("^!", "", codeRegex), code, perl = TRUE)
  any(grepl(factRegex, facts, perl = TRUE)) && (if (neg) !inCode else inCode)
}

human <- "Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03"
rctd <- list(rctdFile = "", rctdReference = "someRef", doSPLIT = TRUE, splitMode = "neighborhood",
             coocFdr = TRUE)
anchors <- list(
  ## class, param, fact regex, source file, code regex
  list("EzAppScSeurat", list(), "seed was set to 38", "R/app-ScSeurat.R", "set\\.seed\\(38\\)"),
  list("EzAppScSeurat", list(), "seed\\.use = 38", "R/seuratUtils.R", "seed\\.use = 38"),
  list("EzAppScSeurat", list(), "niters = 1e5", "R/app-ScSeurat.R", "emptyDrops\\([^)]*niters = 1e5"),
  list("EzAppScSeurat", list(), "clusters = TRUE", "R/app-ScSeurat.R", "scDblFinder\\([^)]*clusters = TRUE"),
  list("EzAppScSeurat", list(), "vst\\.flavor = v2", "R/seuratUtils.R", "vst\\.flavor = \"v2\""),
  list("EzAppScSeurat", list(), "FindClusters algorithm 1", "R/seuratUtils.R",
       "FindClusters\\(\\s*object = scData,\\s*resolution = myResolutions,\\s*verbose = FALSE\\s*\\)"),
  list("EzAppScSeurat", list(), "resolutions 0\\.2, 0\\.4, 0\\.6, 0\\.8 and 1\\.0", "R/seuratUtils.R",
       "seq\\(from = 0\\.2, to = 1, by = 0\\.2\\)"),
  list("EzAppScSeurat", list(), "only\\.pos = TRUE\\); p-values were Bonferroni", "R/seuratUtils.R", "only\\.pos = TRUE"),
  list("EzAppScSeurat", list(refBuild = human, computePathwayTFActivity = TRUE), "times = 100, minsize = 5",
       "R/app-ScSeurat.R", "times = 100,\\s*minsize = 5"),
  list("EzAppScSeurat", list(mLLMCelltype = TRUE), "temperature 0, seed 42", "R/app-ScSeurat.R",
       "temperature = 0,[^\\n]*\\n\\s*seed = 42L"),
  list("EzAppScSeurat", list(CyteTypeR = TRUE), "p below 0\\.05 and log2 fold change above 0\\.5", "R/app-ScSeurat.R",
       "p_val_adj < 0\\.05, avg_log2FC > 0\\.5"),
  list("EzAppSpatialSeurat", list(), "seed was set to 38", "R/app-SpatialSeurat.R", "set\\.seed\\(38\\)"),
  list("EzAppSpatialSeurat", list(), "r\\.metric = 5", "R/seuratUtils.R", "r\\.metric = 5"),
  list("EzAppSpatialSeurat", list(), "below 0\\.01 \\(pvalue_allMarkers", "R/app-SpatialSeurat.R",
       "pvalue_allMarkers = ezFrame\\([^)]*DefaultValue = 0\\.01"),
  list("EzAppSpatialSeuratHD", list(), "seed was set to 38", "R/app-SpatialSeuratHD.R", "set\\.seed\\(38\\)"),
  list("EzAppScSeuratCombine", list(), "seed was set to 38", "R/app-ScSeuratCombine.R", "set\\.seed\\(38\\)"),
  list("EzAppScSeuratCombine", list(integrationMethod = "RPCA"), "k\\.anchor = 20", "R/seuratUtils.R", "k\\.anchor = 20"),
  list("EzAppXeniumSeurat", rctd, "seed was set to 42 immediately before RCTD", "R/app-XeniumSeurat.R", "set\\.seed\\(42\\)"),
  list("EzAppXeniumSeurat", list(), "min\\.pct 0\\.25 and logfc\\.threshold 0\\.25", "R/app-XeniumSeurat.R",
       "min\\.pct = 0\\.25,\\s*logfc\\.threshold = 0\\.25"),
  list("EzAppXeniumSeurat", list(), "first 30 \\(fixed in the code\\)", "R/app-XeniumSeurat.R",
       "FindNeighbors\\(sdata, dims = 1:30"),
  list("EzAppXeniumSeurat", rctd, "at most 10,000 cells per cell type", "R/app-XeniumSeurat.R",
       "spacexr::Reference\\(ref_counts, ref_celltypes\\)"),
  list("EzAppXeniumSeurat", rctd, "Benjamini-Hochberg-adjusted", "inst/templates/XeniumSeurat.Rmd",
       "p\\.adjust\\(P\\[ut\\], method = \"BH\"\\)"),
  list("EzAppXeniumSeurat", rctd, "fewer than 20 cells dropped, a random subsample of 100,000 cells when larger \\(seed 42\\)",
       "R/xeniumCooccurrence.R", "min_cells = 20,\\s*max_cells = 1e5, seed = 42"),
  list("EzAppXeniumSeurat", rctd, "20 nearest spatial neighbours \\(edges over 15 um pruned\\)", "R/app-XeniumSeurat.R",
       "rad_pruning = 15, k_knn = 20"),
  list("EzAppVisiumHDSeurat", list(), "seed was set to 38", "R/app-VisiumHDSeurat.R", "set\\.seed\\(38\\)"),
  list("EzAppVisiumHDSeurat", list(), "fewer than 50,000 bins", "R/app-VisiumHDSeurat.R", "nrow\\(scData@meta\\.data\\) < 50000"),
  list("EzAppVisiumHDSeurat", list(), "PCA computed 80 components", "R/app-VisiumHDSeurat.R", "RunPCA\\(scData, npcs = 80\\)"),
  list("EzAppVisiumHDSeurat", list(), "k_geom = 30 \\(fixed", "R/app-VisiumHDSeurat.R", "k_geom = 30"),
  list("EzAppVisiumHDSeurat", list(), "the first 12 \\(fixed", "R/app-VisiumHDSeurat.R", "dims = 1:12"),
  list("EzAppVisiumHDSeurat", list(), "at resolution 0\\.5 \\(nicheResolution\\)", "R/app-VisiumHDSeurat.R",
       "nicheResolution = ezFrame\\([^)]*DefaultValue = 0\\.5"),
  list("EzAppDeseq2", list(), "minReplicatesForReplace = Inf", "R/twoGroups.R", "minReplicatesForReplace = Inf"),
  list("EzAppDeseq2", list(), "Benjamini-Hochberg adjustment of the DESeq2 Wald p-values", "R/twoGroups.R",
       "p\\.adjust\\(pValue\\[useProbe\\], method = \"fdr\"\\)"),
  list("EzAppDeseq2", list(useLfcShrink = TRUE), "lfcShrink using the ashr method", "R/twoGroups.R",
       "lfcShrink\\(dds[\\s\\S]{0,200}?type=\"ashr\""),
  list("EzAppDeseq2", list(), "sigThresh \\(ezRun default 10\\)", "inst/extdata/EZ_PARAM_DEFAULTS.txt",
       "(?m)^sigThresh\\tnumeric\\t10\\t"),
  list("EzAppEdger", list(testMethod = "glm", deTest = "QL"), "dispersion was fixed at 0\\.1", "R/twoGroups.R",
       "common\\.dispersion <- 0\\.1"),
  list("EzAppEdger", list(testMethod = "glm", deTest = "QL"), "glmQLFit, glmQLFTest", "R/twoGroups.R", "glmQLFit\\("),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "150 training epochs", "R/app-cellBender.R", "!--epochs"),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "\\(--fpr\\) of 0\\.01", "R/app-cellBender.R", "!--fpr"),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "--cpu-threads \\(gpu 0\\)", "R/app-cellBender.R",
       "'--cpu-threads', param\\$cores"),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 1), "--cuda \\(gpu above 0\\)", "R/app-cellBender.R",
       "param\\$gpu > 0\\) \\{\\s*cmd <- paste\\(cmd, \"--cuda\"\\)"),
  list("EzAppSTAR", list(), "infer_experiment\\.py on 1,000,000 sampled reads", "R/app-mapping.R", "\"-s 1000000\""),
  list("EzAppFastqc", structure(list(), input = data.frame(`Read Count` = "2000000000", check.names = FALSE)), "ShortRead FastqSampler, seed 123", "R/fastqIO.R",
       "subsampleFastqFile <- function\\([^)]*seed = 123L\\)")
)

test_that("every fact that states a literal code value still matches the code", {
  codes <- list()
  for (a in anchors) {
    cls <- a[[1]]; file <- a[[4]]
    facts <- get(cls)$new()$methods_facts(a[[2]])
    expect_true(any(grepl(a[[3]], facts, perl = TRUE)), label = paste(cls, "fact:", a[[3]]))
    if (is.null(codes[[file]])) codes[[file]] <- codeOf(file)
    neg <- startsWith(a[[5]], "!")
    expect_identical(grepl(sub("^!", "", a[[5]]), codes[[file]], perl = TRUE), !neg,
                     label = paste(cls, "code:", file, a[[5]]))
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
