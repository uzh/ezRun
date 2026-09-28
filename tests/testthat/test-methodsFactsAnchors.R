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
mouse <- "Mus_musculus/GENCODE/GRCm39/Annotation/Release_M37-2025-07-01"
bdIn <- structure(list(runWNN = TRUE), input = data.frame(`SCDataOrigin [Factor]` = "BDRhapsody", check.names = FALSE))
atac <- list(refBuild = human)
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
  list("EzAppScSeurat", list(), "only\\.pos = TRUE\\), testing genes", "R/seuratUtils.R", "only\\.pos = TRUE", "getSeuratMarkers"),
  list("EzAppScSeurat", list(), "at least a fraction 0\\.1 of either group \\(min\\.pct\\)", "R/seuratUtils.R",
       "min\\.pct = ifelse\\(ezIsSpecified\\(param\\$min\\.pct\\), param\\$min\\.pct, 0\\.1\\)", "getSeuratMarkers"),
  list("EzAppScSeurat", list(), "at least 0\\.25 \\(logfc\\.threshold\\)", "R/seuratUtils.R",
       "param\\$logfc\\.threshold,\\s*0\\.25\\s*\\)", "getSeuratMarkers"),
  list("EzAppCellRanger", list(TenXLibrary = "GEX", expectedCells = ""), "no --expect-cells", "R/app-cellRanger.R",
       "if \\(ezIsSpecified\\(param\\$expectedCells\\)\\) \\{\\s*paste0\\(\"--expect-cells=\"", "ezMethodCellRanger"),
  list("EzAppCellRangerMulti", list(TenXLibrary = "GEX", expectedCells = ""), "no expect-cells in config\\.csv", "R/app-cellRangerMulti.R",
       "if \\(ezIsSpecified\\(param\\$expectedCells\\)\\) \\{[\\s\\S]{0,80}\"expect-cells,%s\"", "buildMultiConfigFile"),
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
       "subsampleFastqFile <- function\\([^)]*seed = 123L\\)", "subsampleFastqFile"),
  ## -- added for the apps of the evaluation matrix (ScSeurat, ScMultiOmics, CellBender, CellRanger,
  ## CellRangerMulti, DESeq2, edgeR, limma, FastQC, hifiasm, STAR with the fastp helper, kallisto)
  ## ScSeurat
  list("EzAppScSeurat", list(refBuild = human, enrichrDatabase = "Tabula_Sapiens"), "adjusted p-value below 0\\.001 and more than 3 overlapping genes were kept, the top 5",
       "R/app-ScSeurat.R", "overlapGeneCutOff = 3,\\s*adjPvalueCutOff = 0\\.001,\\s*reportTopN = 5[\\s\\S]*Adjusted\\.P\\.value < adjPvalueCutOff[\\s\\S]*OverlapGenesN > overlapGeneCutOff[\\s\\S]*head\\(reportTopN\\)",
       "querySignificantClusterAnnotationEnrichR"),
  list("EzAppScSeurat", list(refBuild = human, tissue = "Blood"), "CellMarker 2\\.0 \\(2023-09-27 release\\) gene sets of at least 3 genes",
       "R/scTools.R", "CellMarker_2\\.0-2023-09-27/Cell_marker_All\\.txt", "createCellMarker2_GeneSets"),
  list("EzAppScSeurat", list(refBuild = human, tissue = "Blood"), "gene sets of at least 3 genes", "R/scTools.R", "minGsSize = 3", "cellsLabelsWithAUC"),
  list("EzAppScSeurat", list(sctype.enabled = TRUE, sctype.tissue = "auto"), "ScTypeDB_full marker database fetched from the sc-type GitHub repository at run time, for the tissue Immune system",
       "R/app-ScSeurat.R", "tissue_type <- \"Immune system\"[\\s\\S]*sc-type/master/ScTypeDB_full\\.xlsx", "ezMethodScSeurat"),
  list("EzAppScSeurat", list(mLLMCelltype = TRUE), "from their 10 top markers by average log2 fold change", "R/app-ScSeurat.R",
       "nTopMarkers <- 10\\s*topMarkers <- markers %>%\\s*group_by\\(cluster\\) %>%\\s*slice_max\\(n = nTopMarkers, order_by = avg_log2FC\\)", "ezMethodScSeurat"),
  list("EzAppScSeurat", list(CyteTypeR = TRUE), "clusters with fewer than 5 such markers were not submitted", "R/app-ScSeurat.R",
       "dplyr::filter\\(n >= 5\\)", "ezMethodScSeurat"),
  list("EzAppScSeurat", list(AzimuthPanHuman = TRUE, refBuild = human), "covered at least 50% of the barcodes", "R/app-ScSeurat.R",
       "length\\(shared\\) >= 0\\.5 \\* ncol\\(scData\\)", "ezMethodScSeurat"),
  list("EzAppScSeurat", list(), "\\^MT-, \\^RPS/\\^RPL and \\^HB\\[\\^P\\] \\(case-insensitive\\)", "R/app-ScSeurat.R",
       "\"\\(\\?i\\)\\^MT-\"[\\s\\S]*\"\\(\\?i\\)\\^RPS\\|\\^RPL\"[\\s\\S]*\"\\(\\?i\\)\\^HB\\[\\^P\\]\"", "addCellQcToSeurat"),
  list("EzAppScSeurat", list(), "run once without mitochondrial and ribosomal-protein genes and once on all genes, keeping the larger p-value", "R/app-ScSeurat.R",
       "emptyDrops\\(\\s*rawCts\\[!featInfo\\$isMito & !featInfo\\$isRiboprot, \\][\\s\\S]*emptyDrops\\(rawCts,[\\s\\S]*pmin\\(", "ezMethodScSeurat"),
  list("EzAppScSeurat", list(), "if scDblFinder failed twice, no doublet filtering was applied", "R/app-ScSeurat.R",
       "doubletsInfo <- callScDblFinder\\(\\)\\s*if \\(is\\.null\\(doubletsInfo\\)\\)[\\s\\S]*doubletsInfo <- callScDblFinder\\(nfeatures", "addCellQcToSeurat"),
  list("EzAppScSeurat", list(refBuild = human), "scran cyclone", "R/scTools.R", "method = c\\(\"cyclone\", \"seurat\"\\)", "addCellCycleToSeurat"),
  list("EzAppScSeurat", list(), "a LogNormalize layer, scale factor 10000", "R/seuratUtils.R",
       "normalization\\.method = \"LogNormalize\",\\s*scale\\.factor = 10000", "seuratStandardSCTPreprocessing"),
  list("EzAppScSeurat", list(SCT.regress.CellCycle = TRUE), "cell-cycle score difference \\(S minus G2M\\)", "R/seuratUtils.R",
       "vars\\.to\\.regress <- c\\(\"CC\\.Difference\"\\)", "getSeuratVarsToRegress"),
  list("EzAppScSeurat", list(refBuild = human, SingleR = "HumanPrimaryCellAtlasData"), "fine labels \\(label\\.fine\\)", "R/scTools.R",
       "labels = ref\\$label\\.fine", "cellsLabelsWithSingleR"),
  list("EzAppScSeurat", list(estimateAmbient = TRUE), "DecontX using the cluster labels", "R/sc-estimateAmbient.R",
       "decontX\\(sce, z = sce\\$clusters\\)", "addAmbientEstimateToSeurat"),
  ## ScMultiOmics
  list("EzAppScMultiOmics", list(adtNorm = "CLR"), "CLR-normalized per cell \\(Seurat NormalizeData, margin = 2\\)", "R/multiOmicsUtils.R",
       "normalization\\.method = \"CLR\", margin = 2", "processADT"),
  list("EzAppScMultiOmics", list(adtNorm = "ADTnorm"), "exclude_zeroes = TRUE", "R/multiOmicsUtils.R", "exclude_zeroes = TRUE", "processADT"),
  list("EzAppScMultiOmics", list(), "Non-finite normalized or scaled ADT values were set to 0", "R/multiOmicsUtils.R",
       "norm_mat\\[!is\\.finite\\(norm_mat\\)\\] <- 0[\\s\\S]*sd_mat\\[!is\\.finite\\(sd_mat\\)\\] <- 0", "processADT"),
  list("EzAppScMultiOmics", list(), "exact PCA \\(approx = FALSE\\)[\\s\\S]*30 neighbours", "R/multiOmicsUtils.R",
       "approx = FALSE[\\s\\S]*n_neighbors <- min\\(30L", "processADT"),
  list("EzAppScMultiOmics", bdIn, "30 PCs, a k = 20 shared-nearest-neighbour graph on those 30 PCs and Louvain clustering at a fixed resolution of 0\\.5",
       "R/multiOmicsUtils.R", "RunPCA\\(obj, npcs = 30[\\s\\S]*FindNeighbors\\(obj, dims = 1:30[\\s\\S]*FindClusters\\(obj, resolution = 0\\.5", "loadBDRhapsody"),
  list("EzAppScMultiOmics", list(), "min\\.cells = 0\\), Signac FindTopFeatures used min\\.cutoff q5", "R/multiOmicsUtils.R",
       "min\\.cells = 0[\\s\\S]*FindTopFeatures\\(obj, min\\.cutoff = \"q5\"\\)", "processATAC"),
  list("EzAppScMultiOmics", list(), "ATAC UMAP used LSI components 2 to 30", "R/multiOmicsUtils.R",
       "n_lsi <- min\\(30L[\\s\\S]*dims = 2:n_lsi", "processATAC"),
  list("EzAppScMultiOmics", atac, "EnsDb\\.Hsapiens\\.v86", "R/multiOmicsUtils.R", "\"EnsDb\\.Hsapiens\\.v86\"", "getATACAnnotation"),
  list("EzAppScMultiOmics", list(refBuild = mouse), "EnsDb\\.Mmusculus\\.v79", "R/multiOmicsUtils.R", "\"EnsDb\\.Mmusculus\\.v79\"", "getATACAnnotation"),
  list("EzAppScMultiOmics", atac, "scale factor equal to the median total gene-activity count", "R/multiOmicsUtils.R",
       "scale\\.factor = median\\(obj\\$nCount_GeneActivity\\)", "processATAC"),
  list("EzAppScMultiOmics", list(runWNN = TRUE), "RNA PCs 1-20 \\(fixed, independent of upstream npcs\\), ADT PCs 1 to at most 18, and ATAC LSI components 2-30",
       "R/multiOmicsUtils.R", "dims_rna = 1:20[\\s\\S]*seq_len\\(min\\(18L,[\\s\\S]*2:min\\(30L,", "runWNN"),
  list("EzAppScMultiOmics", list(runWNN = TRUE), "SLM algorithm \\(FindClusters algorithm = 3", "R/multiOmicsUtils.R",
       "FindClusters\\(obj, graph\\.name = \"wsnn\", algorithm = 3", "runWNN"),
  list("EzAppScMultiOmics", list(vdjChain = "auto"), "removeNA, removeMulti and filterMulti all FALSE", "R/multiOmicsUtils.R",
       "removeNA = FALSE, removeMulti = FALSE,\\s*filterMulti = FALSE", "processVDJ"),
  list("EzAppScMultiOmics", list(), "Single \\(1 cell\\), Small \\(2-5\\), Medium \\(6-20\\), Large \\(21-100\\) and Hyperexpanded \\(101-500\\)", "R/multiOmicsUtils.R",
       "proportion = FALSE,\\s*cloneSize = c\\(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500\\)", "processVDJ"),
  ## CellBender
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "the full model \\(ambient RNA plus barcode swapping\\)", "R/app-cellBender.R", "!--model", NA),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "a learning rate of 1e-4", "R/app-cellBender.R", "!--learning-rate", NA),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "a 64-dimensional latent space", "R/app-cellBender.R", "!--z-dim", NA),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "estimated by CellBender's own heuristic", "R/app-cellBender.R", "!--expected-cells|--total-droplets-included", "ezMethodCellBender"),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "MCKP estimator", "R/app-cellBender.R", "!--estimator", NA),
  ## CellRanger
  list("EzAppCellRanger", list(TenXLibrary = "GEX", includeIntrons = TRUE), "--include-introns=true", "R/app-cellRanger.R", "\"--include-introns=true\"", "ezMethodCellRanger"),
  list("EzAppCellRanger", list(TenXLibrary = "GEX", chemistry = "auto"), "--chemistry=auto", "R/app-cellRanger.R", "paste0\\(\"--chemistry=\", param\\$chemistry\\)", "ezMethodCellRanger"),
  list("EzAppCellRanger", list(nReads = 1000), "seqtk sample \\(seed 42, two-pass mode\\)", "R/app-cellRanger.R", "\"seqtk sample -s 42 -2\"", "subsample"),
  list("EzAppCellRanger", list(TenXLibrary = "GEX", keepAlignment = TRUE), "converted to CRAM with samtools", "R/app-cellRanger.R",
       "doCramConversion <- !ezIsSpecified\\(param\\$controlSeqs\\) &&\\s*!ezIsSpecified\\(param\\$secondRef\\)[\\s\\S]*'-C'", "ezMethodCellRanger"),
  list("EzAppCellRanger", list(TenXLibrary = "GEX", keepAlignment = TRUE, bamStats = TRUE), "more than 3 bases of TSO[\\s\\S]*more than 20 million alignments per GB", "R/app-cellRanger.R",
       "nAlign / ram > 20e6[\\s\\S]*x > 3, cb", "computeBamStatsSC"),
  list("EzAppCellRanger", list(TenXLibrary = "GEX", CellRangerVersion = "Aligner/CellRanger/10.1.0", refBuild = human), "reference\\.json declares the genome as GRCh38", "R/app-cellRanger.R",
       "crVersion < numeric_version\\(\"10\\.1\\.0\"\\)[\\s\\S]*ref\\$genomes <- \"GRCh38\"", "cellRangerAnnotatableRef"),
  list("EzAppCellRanger", list(TenXLibrary = "VDJ"), "cellranger mkvdjref", "R/app-cellRanger.R", "\"cellranger mkvdjref\"", "getCellRangerVDJReference"),
  list("EzAppCellRanger", list(TenXLibrary = "GEX"), "built by ezRun with mkref", "R/app-cellRanger.R", "\"cellranger mkref\"", "getCellRangerGEXReference"),
  list("EzAppCellRanger", list(runVeloCyto = TRUE), "velocyto run10x", "R/app-cellRanger.R", "'velocyto run10x'", "ezMethodCellRanger"),
  ## CellRangerMulti
  list("EzAppCellRangerMulti", list(nReads = 1000), "seqtk sample \\(seed 42, two-pass mode\\)", "R/app-cellRanger.R", "\"seqtk sample -s 42 -2\"", "subsample"),
  list("EzAppCellRangerMulti", list(TenXLibrary = "GEX"), "built by ezRun with mkref", "R/app-cellRanger.R", "\"cellranger mkref\"", "getCellRangerGEXReference"),
  list("EzAppCellRangerMulti", list(TenXLibrary = "fixedRNA", probesetFile = "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv", chemistry = "auto"),
       "SFRP \\(singleplex Flex v1\\)", "R/app-cellRangerMulti.R", "\"SFRP\"", "buildMultiConfigFile"),
  list("EzAppCellRangerMulti", list(TenXLibrary = "fixedRNA,Multiplexing", probesetFile = "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv", chemistry = "auto"),
       "MFRP \\(multiplexed Flex v1\\)", "R/app-cellRangerMulti.R", "\"MFRP\"", "buildMultiConfigFile"),
  ## DESeq2 / edgeR / limma
  list("EzAppDeseq2", list(), "median-ratio method on the present genes only \\(controlGenes\\)", "R/twoGroups.R",
       "estimateSizeFactors\\(dds, controlGenes = isPresent\\)", "runDeseq2"),
  list("EzAppDeseq2", list(grouping2 = "Batch"), "~ grouping \\+ grouping2 \\(additive, no interaction\\)", "R/twoGroups.R", "design = ~ grouping \\+ grouping2", "runDeseq2"),
  list("EzAppDeseq2", list(), "results\\(\\) contrast", "R/twoGroups.R", "contrast = c\\(\"grouping\", sampleGroup, refGroup\\)", "runDeseq2"),
  list("EzAppDeseq2", list(runGO = TRUE, rankMetric = "log2Ratio"), "at least 3 genes per term", "R/go-analysis.R",
       "minGenesOverlap <- 3[\\s\\S]*tempTable\\$Count >= minGenesOverlap", "ezEnricher"),
  list("EzAppDeseq2", list(runGO = TRUE, rankMetric = "log2Ratio"), "p-value at or below pValThreshGO and absolute log2 ratio above log2RatioThreshGO", "R/go-analysis.R",
       "pValue <= param\\$pValThreshGO[\\s\\S]*log2Ratio > param\\$log2RatioThreshGO[\\s\\S]*log2Ratio < -param\\$log2RatioThreshGO", "compileEnrichmentInput"),
  list("EzAppDeseq2", list(runGO = TRUE, rankMetric = "pValue"), "ranked by -log10 p-value \\(rankMetric pValue\\)", "R/go-analysis.R",
       "rankMetric == \"pValue\"\\) \\{\\s*df\\$rank_metric <- -log10\\(", "ezGSEA"),
  list("EzAppEdger", list(testMethod = "glm"), "no-intercept design ~0 \\+ group\\.", "R/twoGroups.R", "model\\.matrix\\(~ 0 \\+ groupFactor\\)", "runGlm"),
  list("EzAppEdger", list(testMethod = "glm"), "glm contrast sample group minus reference group", "R/twoGroups.R", "contrastsIndices <- c\\(-1, 1, rep\\(0, ncol\\(fitGlm\\) - 2\\)\\)", "runGlm"),
  list("EzAppEdger", list(testMethod = "glm", deTest = "LR"), "glmFit and glmLRT", "R/twoGroups.R", "glmFit\\([\\s\\S]*glmLRT\\(", "runGlm"),
  list("EzAppEdger", list(testMethod = "glm", deTest = "LR", robust = TRUE), "estimateGLMRobustDisp\\) was used", "R/twoGroups.R", "estimateGLMRobustDisp\\(cds, design\\)", "runGlm"),
  list("EzAppEdger", list(testMethod = "exactTest"), "edgeR exact test", "R/twoGroups.R", "exactTest\\(", "runEdger"),
  list("EzAppEdger", list(testMethod = "glm"), "calcNormFactors using the normMethod method", "R/twoGroups.R", "calcNormFactors\\(cds, method = normMethod\\)", c("runEdger", "runGlm")),
  list("EzAppLimma", list(), "TMM normalization factors \\(edgeR calcNormFactors default\\)", "R/twoGroups.R", "calcNormFactors\\(cds\\)", "runLimma"),
  list("EzAppLimma", list(), "design ~ group with the reference group as baseline", "R/twoGroups.R", "model\\.matrix\\(~groupFactor\\)", "runLimma"),
  list("EzAppLimma", list(grouping2 = "Batch"), "consensus within-block correlation from duplicateCorrelation", "R/twoGroups.R",
       "duplicateCorrelation\\([\\s\\S]*correlation = corfit\\$consensus", "runLimma"),
  list("EzAppLimma", list(modelMethod = "limma-trend"), "eBayes\\(trend = TRUE\\)", "R/twoGroups.R", "eBayes\\(fit, trend = TRUE\\)", "runLimma"),
  list("EzAppLimma", list(modelMethod = "limma-trend"), "log2 CPM with cpm \\(TMM-normalized library sizes, prior\\.count = priorCount\\)", "R/twoGroups.R",
       "cpm\\(cds, log = TRUE, prior\\.count = priorCount\\)", "runLimma"),
  list("EzAppLimma", list(modelMethod = "voom"), "transformed with voom", "R/twoGroups.R", "voom\\(cds, design, plot = FALSE\\)", "runLimma"),
  ## FastQC
  list("EzAppFastqc", list(), "FGCZ adapter list for the adapter content module \\(-a\\)", "R/app-fastQC.R", "\"-a\",\\s*FASTQC_ADAPTERS", "ezMethodFastQC"),
  list("EzAppFastqc", structure(list(), input = data.frame(`Read Count` = "2000000000", check.names = FALSE)), "more than 1 billion reads, FastQC was run on a random subsample of 1,000,000 reads",
       "R/app-fastQC.R", "sum\\(dataset\\$`Read Count`\\) > 1e9\\) \\{\\s*input <- ezMethodSubsampleFastq\\(input = input, param = param, n = 1e6\\)", "ezMethodFastQC"),
  list("EzAppFastqc", list(), "combined into one report with MultiQC", "R/app-fastQC.R", "multiqc", "ezMethodFastQC"),
  list("EzAppFastqc", list(generate_ai_summary = TRUE), "--ai-summary-full", "R/app-fastQC.R", "--ai-summary-full", "ezMethodFastQC"),
  ## Hifiasm
  list("EzAppHifiasm", list(), "--primary mode[\\s\\S]*--n-hap set to the ploidy", "R/app-hifiasm.R", "paste\\(\"--n-hap\", param\\$ploidy\\),\\s*\"--primary\"", "ezMethodHifiasm"),
  list("EzAppHifiasm", list(inputType = "ONT"), "\\(--ont\\)", "R/app-hifiasm.R", "\"--ont\"", "ezMethodHifiasm"),
  list("EzAppHifiasm", list(), "segment \\(S\\) lines of hifiasm's p_ctg\\.gfa", "R/app-hifiasm.R", "p_ctg\\.gfa", "ezMethodHifiasm"),
  ## STAR (+ fastp helper)
  list("EzAppSTAR", list(), "--outSAMattributes All appended to cmdOptions \\(unless cmdOptions already set --outSAMattributes\\)", "R/app-mapping.R",
       "if \\(!str_detect\\(param\\$cmdOptions, \"outSAMattributes\"\\)\\) \\{\\s*param\\$cmdOptions <- str_c\\(\\s*param\\$cmdOptions,\\s*\"--outSAMattributes All\"", "ezMethodSTAR"),
  list("EzAppSTAR", list(twopassMode = TRUE), "--twopassMode Basic", "R/app-mapping.R", "\"--twopassMode\",\\s*if_else\\(param\\$twopassMode, \"Basic\", \"None\"\\)", "ezMethodSTAR"),
  list("EzAppSTAR", list(twopassMode = FALSE), "--twopassMode None", "R/app-mapping.R", "\"--twopassMode\",\\s*if_else\\(param\\$twopassMode, \"Basic\", \"None\"\\)", "ezMethodSTAR"),
  list("EzAppSTAR", list(markDuplicates = TRUE), "REMOVE_DUPLICATES=false, OPTICAL_DUPLICATE_PIXEL_DISTANCE set to dupDistance, ezRun default 2500", "inst/extdata/EZ_PARAM_DEFAULTS.txt",
       "(?m)^dupDistance\\tnumeric\\t2500\\t", NA),
  list("EzAppSTAR", list(markDuplicates = TRUE), "REMOVE_DUPLICATES=false, OPTICAL_DUPLICATE_PIXEL_DISTANCE set to dupDistance", "R/bamUtils.R",
       "paste0\\(\"OPTICAL_DUPLICATE_PIXEL_DISTANCE=\", dupDistance\\)[\\s\\S]*\"REMOVE_DUPLICATES=\",\\s*if_else\\(operation == \"mark\", \"false\", \"true\"\\)", "dupBam"),
  list("EzAppSTAR", list(markDuplicates = TRUE), "flagged, not removed", "R/app-mapping.R",
       "dupBam\\([\\s\\S]{0,120}?operation = \"mark\",[\\s\\S]{0,60}?dupDistance = param\\$dupDistance", "ezMethodSTAR"),
  list("EzAppSTAR", list(barcodePattern = "NNNNNNNN"), "umi_tools dedup \\(default directional method", "R/app-mapping.R", "'umi_tools dedup --temp-dir=\\. --verbose=0'", "ezMethodSTAR"),
  list("EzAppSTAR", list(secondRef = "x.fa"), "STAR --genomeFastaFiles", "R/app-mapping.R", "--genomeFastaFiles", "ezMethodSTAR"),
  list("EzAppSTAR", list(trimAdapter = TRUE), "allIllumina-forTrimmomatic-20160202\\.fa", "inst/extdata/EZ_GLOBAL_VARIABLES.txt",
       "TRIMMOMATIC_ADAPTERS=\"[^\"]*/allIllumina-forTrimmomatic-20160202\\.fa\"", NA),
  list("EzAppSTAR", list(trimAdapter = FALSE), "--disable_quality_filtering\\) because average_qual was not set", "R/app-trim.R",
       "if \\(ezIsSpecified\\(param\\$average_qual\\)\\) \\{\\s*paste\\(\"--average_qual\", param\\$average_qual\\)\\s*\\} else \\{\\s*\"--disable_quality_filtering\"", "ezMethodFastpTrim"),
  list("EzAppSTAR", list(poly_x_min_len = 10), "fastp --trim_poly_x", "R/app-trim.R", "paste\\(\"--trim_poly_x\", \"--poly_x_min_len\", param\\$poly_x_min_len\\)", "ezMethodFastpTrim"),
  ## Kallisto
  list("EzAppKallisto", list(strandMode = "sense"), "--fr-stranded", "R/app-kallisto.R", "\"sense\" = \"--fr-stranded\"", "ezMethodKallisto"),
  list("EzAppKallisto", list(strandMode = "antisense"), "--rf-stranded", "R/app-kallisto.R", "\"antisense\" = \"--rf-stranded\"", "ezMethodKallisto"),
  list("EzAppKallisto", list(paired = FALSE), "kallisto --single; a fragment-length or sd of 0 was replaced by a mean fragment length of 180 and a standard deviation of 50", "R/app-kallisto.R",
       "param\\$\"fragment-length\" = 180[\\s\\S]*param\\$sd = 50[\\s\\S]*iftrue\\(param\\$paired, \"\", \"--single\"\\)", "ezMethodKallisto"),
  list("EzAppKallisto", list(gpu = 1), "bootstrap samples was forced to 0", "R/app-kallisto.R", "param\\$gpu > 0\\) \\{\\s*param\\$\"bootstrap-samples\" = 0", "ezMethodKallisto"),
  list("EzAppScSeurat", list(Azimuth = "pbmcref"), "Azimuth RunAzimuth to the Azimuth reference named in Azimuth, using the RNA counts", "R/seuratUtils.R",
       "Azimuth::RunAzimuth\\(scData, param\\$Azimuth, assay = \"RNA\"\\)", "getSeuratMarkersAndAnnotate"),
  list("EzAppScMultiOmics", list(), "Spearman correlation of per-cluster average ADT and matched-gene RNA values", "inst/templates/_scMultiOmics_adt.Rmd",
       "stats::cor\\(x, y, method = \"spearman\"\\)", NA),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "the log then says \"Dropping Peaks \\(ATAC\\) for CellBender\"", "R/app-cellBender.R",
       "ezLog\\(\"Dropping Peaks \\(ATAC\\) for CellBender", "dropPeaksFromH5"),
  list("EzAppCellRangerMulti", list(TenXLibrary = "fixedRNA", probesetFile = "Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"),
       "match a gene of the Cell Ranger reference \\(star/geneInfo\\.tab\\)", "R/app-cellRangerMulti.R", "file\\.path\\(refDir, 'star', 'geneInfo\\.tab'\\)", "buildMultiConfigFile"),
  list("EzAppDeseq2", list(runGO = FALSE), "The DESeq2 model used the design ~ grouping\\.", "R/twoGroups.R", "design = ~grouping", "runDeseq2"),
  list("EzAppDeseq2", list(runGO = FALSE), "Cook's distance filtering of p-values was not applied \\(cooksCutoff false", "R/twoGroups.R",
       "cooksCutoff = ezIsSpecified\\(param\\$cooksCutoff\\) && param\\$cooksCutoff", "twoGroupCountComparison"),
  list("EzAppDeseq2", list(runGO = TRUE, rankMetric = "log2Ratio"), "ranked by log2 ratio \\(rankMetric log2Ratio\\)", "R/go-analysis.R",
       "rankMetric == \"log2Ratio\"\\) \\{\\s*df\\$rank_metric <- df\\$log2Ratio", "ezGSEA"),
  list("EzAppDeseq2", list(), "added to the normalized counts before they are log2-transformed", "inst/templates/twoGroups.Rmd",
       "log2\\(assays\\(se\\)\\$xNorm \\+ param\\$backgroundExpression\\)", NA),
  list("EzAppEdger", list(), "added to the normalized counts before they are log2-transformed", "inst/templates/twoGroups.Rmd",
       "log2\\(assays\\(se\\)\\$xNorm \\+ param\\$backgroundExpression\\)", NA),
  list("EzAppLimma", list(), "added to the normalized counts before they are log2-transformed", "inst/templates/twoGroups.Rmd",
       "log2\\(assays\\(se\\)\\$xNorm \\+ param\\$backgroundExpression\\)", NA),
  list("EzAppEdger", list(), "sigThresh \\(ezRun default 10\\)", "inst/extdata/EZ_PARAM_DEFAULTS.txt", "(?m)^sigThresh\\tnumeric\\t10\\t", NA),
  list("EzAppLimma", list(), "sigThresh \\(ezRun default 10\\)", "inst/extdata/EZ_PARAM_DEFAULTS.txt", "(?m)^sigThresh\\tnumeric\\t10\\t", NA),
  list("EzAppEdger", list(testMethod = "glm"), "prior count equal to backgroundExpression \\(prior\\.count\\)", "R/twoGroups.R",
       "exactTest = runEdger\\([\\s\\S]*priorCount = param\\$backgroundExpression[\\s\\S]*glm = runGlm\\([\\s\\S]*priorCount = param\\$backgroundExpression", "twoGroupCountComparison"),
  list("EzAppEdger", list(testMethod = "glm"), "prior count equal to backgroundExpression \\(prior\\.count\\)", "R/twoGroups.R",
       "prior\\.count = priorCount", c("runEdger", "runGlm")),
  list("EzAppSTAR", list(), "computed with dupRadar on the delivered BAM", "R/app-mapping.R", "getDupRateFromBam\\(", "ezMethodSTAR"),
  list("EzAppSTAR", list(), "computed with dupRadar on the delivered BAM", "R/app-RnaBamStats.R", "analyzeDuprates\\(", "getDupRateFromBam"),
  list("EzAppKallisto", list(), "The kallisto index \\(default k-mer length 31\\)", "R/app-kallisto.R", "cmdTemplate = \"kallisto index -i %s\\.idx %s\"", "getKallistoReference"),
  list("EzAppScSeurat", list(), "Cells with zero UMIs were always removed", "R/app-ScSeurat.R",
       "scData\\$qc\\.lib <- scData\\$qc\\.lib \\| scData@meta\\.data\\[, att_nCounts\\] == 0", "addCellQcToSeurat"),
  list("EzAppScSeurat", list(), "Genes with no counts in the filtered matrix were dropped; after cell filtering, genes were additionally required to have at least geneMinUMI UMIs in at least the cellsFraction",
       "R/app-ScSeurat.R", "cts\\[rowSums2\\(cts > 0\\) > 0, \\][\\s\\S]*num\\.cells <- param\\$cellsFraction \\* ncol\\(scData\\)[\\s\\S]*>= param\\$geneMinUMI[\\s\\S]*cellsPerGene >= num\\.cells",
       "ezMethodScSeurat"),
  list("EzAppScMultiOmics", list(runWNN = TRUE), "FindAllMarkers \\(Wilcoxon test, only\\.pos = TRUE, min\\.pct = 0\\.25, logfc\\.threshold = 0\\.25\\)", "inst/templates/_scMultiOmics_wnn.Rmd",
       "FindAllMarkers\\(scMultiData_wm, assay = assay_name, only\\.pos = TRUE,\\s*min\\.pct = 0\\.25, logfc\\.threshold = 0\\.25", NA),
  list("EzAppScMultiOmics", list(vdjChain = "BCR"), "combineBCR on IGH CDR3 nucleotide sequences with length-normalized Levenshtein similarity at bcrSimilarityThreshold", "R/multiOmicsUtils.R",
       "combineBCR\\(contigs, samples = sampleName,\\s*threshold = bcrSimilarityThreshold", "processVDJ"),
  list("EzAppScMultiOmics", list(cloneCallTCR = "strict"), "Clone identity in combineExpression followed cloneCallTCR \\(strict\\) for both TCR and BCR", "R/multiOmicsUtils.R",
       "effectiveClone <- cloneCallTCR[\\s\\S]*cloneCall = effectiveClone", "processVDJ"),
  list("EzAppCellRangerMulti", list(keepBam = FALSE), "The per-sample BAM files were deleted after the run", "R/app-cellRangerMulti.R",
       "if \\(!param\\$keepBam\\) \\{[\\s\\S]{0,200}?-name \"\\*_alignments\\.bam\\*\" -type f -delete", "ezMethodCellRangerMulti"),
  list("EzAppEdger", list(), "calcNormFactors using the normMethod method on all genes\\.", "R/twoGroups.R", "calcNormFactors\\(cds, method = normMethod\\)", c("runEdger", "runGlm")),
  list("EzAppEdger", list(testMethod = "exactTest"), "calcNormFactors using the normMethod method on all genes, on all samples of the dataset", "R/twoGroups.R",
       "cds <- DGEList\\(counts = x, group = grouping\\)\\s*cds <- calcNormFactors\\(cds, method = normMethod\\)", "runEdger"),
  list("EzAppEdger", list(), "prior count equal to backgroundExpression \\(prior\\.count\\)", "R/twoGroups.R", "prior\\.count = priorCount", c("runEdger", "runGlm")),
  list("EzAppEdger", list(runGO = TRUE, rankMetric = "log2Ratio"), "at least 3 genes per term", "R/go-analysis.R",
       "minGenesOverlap <- 3[\\s\\S]*tempTable\\$Count >= minGenesOverlap", "ezEnricher"),
  list("EzAppEdger", list(runGO = TRUE, rankMetric = "log2Ratio"), "ranked by log2 ratio \\(rankMetric log2Ratio\\)", "R/go-analysis.R",
       "rankMetric == \"log2Ratio\"\\) \\{\\s*df\\$rank_metric <- df\\$log2Ratio", "ezGSEA"),
  list("EzAppKallisto", list(strandMode = "both"), "without strand information \\(no kallisto strand option\\)", "R/app-kallisto.R", "\"both\" = \"\"", "ezMethodKallisto"),
  list("EzAppKallisto", list(trimAdapter = TRUE), "allIllumina-forTrimmomatic-20160202\\.fa", "inst/extdata/EZ_GLOBAL_VARIABLES.txt",
       "TRIMMOMATIC_ADAPTERS=\"[^\"]*/allIllumina-forTrimmomatic-20160202\\.fa\"", NA),
  list("EzAppKallisto", list(trimAdapter = FALSE), "--disable_quality_filtering\\) because average_qual was not set", "R/app-trim.R",
       "if \\(ezIsSpecified\\(param\\$average_qual\\)\\) \\{\\s*paste\\(\"--average_qual\", param\\$average_qual\\)\\s*\\} else \\{\\s*\"--disable_quality_filtering\"", "ezMethodFastpTrim"),
  list("EzAppSTAR", list(trimAdapter = FALSE), "Reads shorter than 15 bases after trimming were discarded \\(fastp default length_required, not set by the wrapper\\)", "R/app-trim.R",
       "if \\(ezIsSpecified\\(param\\$length_required\\)\\) \\{\\s*paste\\(\"--length_required\", param\\$length_required\\)", "ezMethodFastpTrim"),
  list("EzAppKallisto", list(trimAdapter = FALSE), "Reads shorter than 15 bases after trimming were discarded \\(fastp default length_required, not set by the wrapper\\)", "R/app-trim.R",
       "if \\(ezIsSpecified\\(param\\$length_required\\)\\) \\{\\s*paste\\(\"--length_required\", param\\$length_required\\)", "ezMethodFastpTrim"),
  list("EzAppScSeurat", list(), "PCA computed 50 components \\(Seurat default\\) and the first npcs were used for the neighbour graph, UMAP and t-SNE", "R/seuratUtils.R",
       "RunPCA\\(object = scData, verbose = FALSE\\)[\\s\\S]*RunTSNE\\([\\s\\S]*dims = 1:param\\$npcs[\\s\\S]*RunUMAP\\(object = scData, reduction = reduction, dims = 1:param\\$npcs\\)",
       "seuratStandardWorkflow"),
  list("EzAppScSeurat", list(), "UMAP used uwot with the cosine metric, 30 neighbours and seed 42 \\(RunUMAP defaults\\)", "R/seuratUtils.R",
       "RunUMAP\\(object = scData, reduction = reduction, dims = 1:param\\$npcs\\)", "seuratStandardWorkflow"),
  list("EzAppScMultiOmics", list(), "taken unchanged from the upstream ScSeurat object \\(Azimuth, SingleR fine-label and cellxgene results saved next to it were re-attached\\)",
       "R/app-ScMultiOmics.R", "qs2::qs_read\\(scDataPath[\\s\\S]*attachUpstreamAnnotations\\(obj, dirname\\(scDataPath\\)\\)", "ezMethodScMultiOmics"),
  list("EzAppScMultiOmics", list(), "ADT counts were the Antibody Capture features of the count-matrix H5", "R/multiOmicsUtils.R",
       "adt <- cts\\[\\[\"Antibody Capture\"\\]\\][\\s\\S]*hashtags <- findCellRangerMultiHashtags\\(h5path\\)", "readADTCounts"),
  list("EzAppScMultiOmics", list(), "removeNA, removeMulti and filterMulti all FALSE", "R/multiOmicsUtils.R",
       "removeNA = FALSE, removeMulti = FALSE,\\s*filterMulti = FALSE", "processVDJ"),
  list("EzAppScMultiOmics", list(vdjChain = "BCR"), "removeNA, removeMulti and filterMulti all FALSE", "R/multiOmicsUtils.R",
       "combineBCR\\(contigs, samples = sampleName,[\\s\\S]{0,120}?removeNA = FALSE, removeMulti = FALSE,\\s*filterMulti = FALSE", "processVDJ"),
  list("EzAppScMultiOmics", list(vdjChain = "TCR"), "removeNA, removeMulti and filterMulti all FALSE", "R/multiOmicsUtils.R",
       "combineTCR\\(contigs, samples = sampleName,[\\s\\S]{0,80}?removeNA = FALSE, removeMulti = FALSE,\\s*filterMulti = FALSE", "processVDJ"),
  list("EzAppScMultiOmics", list(), "Clone identity in combineExpression followed cloneCallTCR for both TCR and BCR", "R/multiOmicsUtils.R",
       "effectiveClone <- cloneCallTCR[\\s\\S]*cloneCall = effectiveClone", "processVDJ"),
  list("EzAppSTAR", list(), "fastp polyG tail trimming was not set by the wrapper", "R/app-trim.R", "!trim_poly_g|poly_g_min_len", "ezMethodFastpTrim"),
  list("EzAppKallisto", list(), "fastp polyG tail trimming was not set by the wrapper", "R/app-trim.R", "!trim_poly_g|poly_g_min_len", "ezMethodFastpTrim"),
  list("EzAppCellBender", list(cmdOptions = "", gpu = 0), "without an UnfilteredCountMatrix column, the library-level raw matrix of the pool, or the sample's raw matrix when the pool has none",
       "R/app-cellBender.R", "\"UnfilteredCountMatrix\" %in% input\\$colNames[\\s\\S]*\"multi/count/raw_feature_bc_matrix\"[\\s\\S]*if \\(!file\\.exists\\(cmDir\\)\\)[\\s\\S]*\"count/sample_raw_feature_bc_matrix\"",
       "ezMethodCellBender")
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
                    "EzAppDeseq2", "EzAppEdger", "EzAppCellBender", "EzAppScMultiOmics", "EzAppCellRanger",
                    "EzAppCellRangerMulti", "EzAppFastqc", "EzAppHifiasm", "EzAppSTAR", "EzAppLimma",
                    "EzAppKallisto") %in% vapply(anchors, `[[`, "", 1)))
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
