## Guards on the LLM-written Methods text, run by write_methods() on the Description:
## numbers the run's configuration does not contain, steps whose parameter was off,
## and a template Description for when the writer fails or keeps failing a guard.

## A standalone number: not glued to letters (GRCh38, CD45, 10x, log2FC, p28409), except
## a unit suffix or %. Thousands commas, scientific forms (1e5, 10^-5) and ranges
## (1-20, 1:20, "1 to 20": each endpoint is its own token) are read as numbers.
METHODS_NUMBER_CORE <- paste0("(?:\\d{1,3}(?:,\\d{3})+(?:\\.\\d+)?|\\d+(?:\\.\\d+)*)",
                              "(?:[eE][-+]?\\d+|\\^[-+]?\\d+)?")
METHODS_NUMBER_TOKEN <- paste0("(?<![A-Za-z0-9_.])(?:[vV](?=\\d))?", METHODS_NUMBER_CORE,
                               "(?:%|\\s?(?:kb|Kb|bp|Mb|Gb)(?![A-Za-z0-9_])|(?:k|K|M)(?![A-Za-z0-9_]))?",
                               "(?![A-Za-z0-9_])")
METHODS_UNIT_SCALE <- c(kb = 1e3, Kb = 1e3, bp = 1, Mb = 1e6, Gb = 1e9, k = 1e3, K = 1e3, M = 1e6)

## The values a number token can stand for: "5%" is 5 or 0.05, "2 kb" is 2000.
methodsNumberValues <- function(tok) {
  x <- sub("^[vV]", "", tok)
  if (endsWith(x, "%")) {
    v <- methodsNumberValues(sub("%$", "", x))
    return(c(v, v / 100))
  }
  unit <- regmatches(x, regexpr("(kb|Kb|bp|Mb|Gb|k|K|M)$", x))
  scale <- if (length(unit)) METHODS_UNIT_SCALE[[unit]] else 1
  x <- gsub(",", "", trimws(sub("(kb|Kb|bp|Mb|Gb|k|K|M)$", "", x)))
  v <- if (grepl("^", x, fixed = TRUE)) {
    p <- as.numeric(strsplit(x, "^", fixed = TRUE)[[1]])
    p[1]^p[2]
  } else suppressWarnings(as.numeric(x))
  v * scale
}

## Every numeric value in the configuration text; loose on purpose (Release_48 gives 48,
## "2kb" gives 2 and 2000), so a configured value is never flagged for its spelling.
methodsConfigValues <- function(config_text) {
  toks <- regmatches(config_text, gregexpr(paste0("(?<![0-9.])", METHODS_NUMBER_CORE,
                                                  "(?:%|\\s?(?:kb|Kb|bp|Mb|Gb|k|K|M))?"),
                                           config_text, perl = TRUE))[[1]]
  vals <- unlist(lapply(toks, function(t) c(methodsNumberValues(t),
                                            methodsNumberValues(sub("[^0-9]+$", "", t)))))
  unique(vals[is.finite(vals)])
}

## Numbers in the Description that the run's configuration does not contain. config_text is
## what configured the run (job script, parameters, app defaults, config.csv, app facts,
## citation candidates); all_text adds the logs and is only used for versions, because logs
## hold results: nearly every number an LLM writes appears somewhere in them.
methods_check_numbers <- function(description, config_text, all_text, sample_count, allow_small = 0:10) {
  text <- paste(description, collapse = "\n")
  m <- gregexpr(METHODS_NUMBER_TOKEN, text, perl = TRUE)[[1]]
  if (m[1] == -1) return(character(0))
  toks <- regmatches(text, list(m))[[1]]
  config <- methodsConfigValues(paste(config_text, collapse = "\n"))
  all <- paste(all_text, collapse = "\n")
  flagged <- vapply(seq_along(toks), function(i) {
    tok <- toks[i]
    core <- sub("^[vV]", "", tok)
    before <- substr(text, max(1, m[i] - 8), m[i] - 1)
    ## x.y.z, v-prefixed, or after "version" / "R": a version, looked up as a string
    if (grepl("^\\d+(\\.\\d+){2,}$", core) || grepl("^[vV]\\d", tok) ||
        grepl("(version|\\bR)[ -]$", before)) {
      return(!grepl(paste0("(?<![0-9])", gsub(".", "\\.", core, fixed = TRUE), "(?![0-9])"),
                    all, perl = TRUE))
    }
    v <- methodsNumberValues(tok)
    v <- v[is.finite(v)]
    if (!length(v)) return(FALSE)
    if (any(v == round(v) & v %in% c(allow_small, sample_count))) return(FALSE)
    !any(vapply(v, function(x) any(abs(config - x) <= 1e-9 * pmax(abs(config), abs(x))), logical(1)))
  }, logical(1))
  unique(toks[flagged])
}

## A parameter value meaning "step off". "None" is rctdReference's off value.
METHODS_OFF_VALUE <- "^(false|FALSE|False|0|none|None|NONE|)$"
METHODS_NEGATION <- "\\b(not|no|without|neither|nor|disabled|skipped|omitted)\\b"

## Per app class: parameter -> lower-case regex of the tool that step runs. When the job's
## value is off, a Description sentence naming the tool, not negated, is flagged.
## "a+b" means the step is off only when both are off. Ported from the facts A/B
## gated.py and extended from the params each methods_facts()/citation() gates on.
METHODS_OFFSTEP_RULES <- list(
  ## annotation, ambient-RNA and pathway steps (facts + citation gates)
  EzAppScSeurat = c(computePathwayTFActivity = "decoupler|dorothea|progeny", SingleR = "singler",
                    estimateAmbient = "decontx|soupx", enrichrDatabase = "enrichr",
                    tissue = "aucell|cellmarker", Azimuth = "(?<!pan-human )azimuth",
                    AzimuthPanHuman = "pan-human azimuth", sctype.enabled = "\\bsctype\\b|sc-type",
                    mLLMCelltype = "mllmcelltype", CyteTypeR = "cytetype",
                    SCT.regress.CellCycle = "regress\\w*[^.;]*cell.cycle|cell.cycle[^.;]*regress"),
  ## integration method and annotation steps
  EzAppScSeuratCombine = c(computePathwayTFActivity = "decoupler|dorothea|progeny", SingleR = "singler",
                           enrichrDatabase = "enrichr", tissue = "aucell|cellmarker",
                           integrationMethod = "harmony|\\bcca\\b|\\brpca\\b|reciprocal pca",
                           SCT.regress.CellCycle = "regress\\w*[^.;]*cell.cycle|cell.cycle[^.;]*regress"),
  ## relabelling only; annotation steps
  EzAppScSeuratCombinedLabelClusters = c(computePathwayTFActivity = "decoupler|dorothea|progeny",
                                         SingleR = "singler", enrichrDatabase = "enrichr",
                                         tissue = "aucell|cellmarker"),
  ## SpotClean, Azimuth and Enrichr
  EzAppSpatialSeurat = c(spotClean = "spotclean", Azimuth = "azimuth", enrichrDatabase = "enrichr",
                         SCT.regress.CellCycle = "regress\\w*[^.;]*cell.cycle|cell.cycle[^.;]*regress"),
  ## batch correction by CCA
  EzAppSpatialSeuratSlides = c(batchCorrection = "\\bcca\\b|integrat",
                               SCT.regress.CellCycle = "regress\\w*[^.;]*cell.cycle|cell.cycle[^.;]*regress"),
  ## SPLIT, co-occurrence FDR, RCTD
  EzAppXeniumSeurat = c(doSPLIT = "\\bsplit\\b(?! into)", coocFdr = "benjamini",
                        "rctdFile+rctdReference" = "\\brctd\\b|spacexr"),
  ## RCTD
  EzAppVisiumHDSeurat = c("rctdFile+rctdReference" = "\\brctd\\b|spacexr|rctd-py"),
  ## GO enrichment, LFC shrinkage
  EzAppDeseq2 = c(runGO = "enricher|\\bgsea\\b|over-representation", useLfcShrink = "\\bashr\\b|lfcshrink"),
  EzAppEdger = c(runGO = "enricher|\\bgsea\\b|over-representation"),
  EzAppLimma = c(runGO = "enricher|\\bgsea\\b|over-representation"),
  ## WNN and TCR clonalCluster
  EzAppScMultiOmics = c(runWNN = "\\bwnn\\b|weighted nearest", tcrSimilarityMerge = "clonalcluster"),
  ## sccomp needs replicateGrouping; the DESeq2 test only runs in pseudobulk mode
  EzAppScSeuratCompare = c(replicateGrouping = "sccomp", pseudoBulkMode = "deseq2|aggregateexpression"),
  ## velocyto, CRAM, FASTQ subsampling
  EzAppCellRanger = c(runVeloCyto = "velocyto", keepAlignment = "\\bcram\\b", nReads = "seqtk|subsampl"),
  EzAppCellRangerMulti = c(keepBam = "\\bcram\\b", nReads = "seqtk|subsampl",
                           customProbesFile = "customprobesfile|custom probe"),
  ## probe set, protein panel, CRAM, TIFF splitting
  EzAppSpaceRanger = c(probesetFile = "probe.set", panelFile = "feature-ref|antibody capture|protein panel",
                       keepAlignment = "\\bcram\\b", splitTif = "tiffsplit"),
  ## two-pass mapping, UMI extraction
  EzAppSTAR = c(twopassMode = "two-pass|twopass", barcodePattern = "umi_tools|\\bumis?\\b"),
  ## deduplication, bigWig
  EzAppBismark = c(deduplicate = "deduplicat", generateBigWig = "bigwig"),
  ## read groups, duplicate marking, BQSR
  EzAppGatkDnaHaplotyper = c(addReadGroup = "addorreplacereadgroups", markDuplicates = "markduplicates",
                             knownSitesAvailable = "baserecalibrator|applybqsr|\\bbqsr\\b|recalibrat"),
  ## VQSR, SnpEff
  EzAppJoinGenoTypes = c("recalibrateVariants+recalibrateInDels" = "\\bvqsr\\b|variantrecalibrator",
                         snpEffDB = "snpeff"),
  ## read-count estimation mode
  EzAppMetaPhlAn = c(estimateReadCounts = "rel_ab_w_read_stats"),
  ## AI summaries in the MultiQC report
  EzAppFastqc = c("generate_ai_summary+per_section_ai_summaries" = "language model|ai-generated|ai summar"),
  ## PEAR merging of pairs
  EzAppSamsa2 = c(paired = "\\bpear\\b"),
  ## GPU
  EzAppCellBender = c(gpu = "cuda")
)

methodsParamOff <- function(param, names) {
  all(vapply(names, function(p) {
    v <- param[[p]]
    !is.null(v) && grepl(METHODS_OFF_VALUE, paste(as.character(unlist(v)), collapse = ","))
  }, logical(1)))
}

## TRUE when some sentence names the tool and not every such sentence is negated.
methodsClaims <- function(text, keyword) {
  sents <- strsplit(tolower(paste(text, collapse = "\n")), "(?<=[.;])\\s+", perl = TRUE)[[1]]
  hit <- sents[grepl(keyword, sents, perl = TRUE)]
  length(hit) > 0 && !all(grepl(METHODS_NEGATION, hit, perl = TRUE))
}

## "param:tool" for each step described in the Description whose parameter was off.
methods_check_offsteps <- function(description, class_name, param) {
  rules <- METHODS_OFFSTEP_RULES[[class_name]]
  if (is.null(rules) || length(param) == 0) return(character(0))
  text <- tolower(paste(description, collapse = "\n"))
  out <- character(0)
  for (p in names(rules)) {
    if (methodsParamOff(param, strsplit(p, "+", fixed = TRUE)[[1]]) && methodsClaims(text, rules[[p]]))
      out <- c(out, paste0(p, ":", regmatches(text, regexpr(rules[[p]], text, perl = TRUE))))
  }
  out
}
