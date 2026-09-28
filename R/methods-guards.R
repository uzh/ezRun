## Guards on the LLM-written Methods text, run by write_methods() on the Description:
## numbers the run's configuration does not contain, steps whose parameter was off,
## and a template Description for when the writer fails or keeps failing a guard.

## A standalone number: not glued to letters (GRCh38, CD45, 10x, log2FC, p28409), except
## a unit suffix or %; nor a partial match of a longer one (2.7 of STAR 2.7.11b). Thousands commas, scientific forms (1e5, 10^-5) and ranges
## (1-20, 1:20, "1 to 20": each endpoint is its own token) are read as numbers.
METHODS_NUMBER_CORE <- paste0("(?:\\d{1,3}(?:,\\d{3})+(?:\\.\\d+)?|\\d+(?:\\.\\d+)*)",
                              "(?:[eE][-+]?\\d+|\\^[-+]?\\d+)?")
METHODS_NUMBER_TOKEN <- paste0("(?<![A-Za-z0-9_.])(?:[vV](?=\\d))?", METHODS_NUMBER_CORE,
                               "(?:%|\\s?(?:kb|Kb|bp|Mb|Gb)(?![A-Za-z0-9_])|(?:k|K|M)(?![A-Za-z0-9_]))?",
                               "(?![A-Za-z0-9_]|[.,]\\d)")
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
  ## digit runs inside identifiers too (dbsnp.157.gencode_compatible gives 157)
  runs <- as.numeric(regmatches(config_text, gregexpr("[0-9]+", config_text))[[1]])
  unique(c(as.numeric(vals[is.finite(vals)]), runs))
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
## "none"/"= false"/"switched off": the tool named with its off value ("SingleR was set to
## 'none'") is a true sentence, not a claim; "false" only as a value, so a false discovery
## rate still reads as a claim.
METHODS_NEGATION <- paste0("\\b(not|no|none|without|neither|nor|disabled|skipped|omitted|single-pass)\\b|",
                           "(=|\\bset to|\\bwas|\\bwere)\\s*['\"]?false\\b|\\b(switched|turned|set) off\\b")

## Per app class: parameter -> lower-case regex of the tool that step runs. When the job's
## value is off, a Description sentence naming the tool, not negated, is flagged.
## "a+b" means the step is off only when both are off. Ported from the facts A/B
## gated.py and extended from the params each methods_facts()/citation() gates on.
## Tool names inside Enrichr library names (Azimuth_Cell_Types_2021, CellMarker_2024) and
## generic words (split, integrated, UMI, Benjamini-Hochberg) are not claims.
METHODS_OFFSTEP_RULES <- list(
  ## annotation, ambient-RNA and pathway steps (facts + citation gates)
  EzAppScSeurat = c(computePathwayTFActivity = "decoupler|dorothea|progeny", SingleR = "singler",
                    estimateAmbient = "decontx|soupx", enrichrDatabase = "enrichr",
                    tissue = "aucell|(?<!_)cellmarker(?!_)", Azimuth = "(?<!pan-human )(?<!_)azimuth(?![_\\w]|[ _]cell[ _]types)",
                    AzimuthPanHuman = "pan-human azimuth", sctype.enabled = "\\bsctype\\b|sc-type",
                    mLLMCelltype = "mllmcelltype", CyteTypeR = "cytetype",
                    SCT.regress.CellCycle = "regress\\w*[^.;]*cell.cycle|cell.cycle[^.;]*regress"),
  ## integration method and annotation steps
  EzAppScSeuratCombine = c(computePathwayTFActivity = "decoupler|dorothea|progeny", SingleR = "singler",
                           enrichrDatabase = "enrichr", tissue = "aucell|(?<!_)cellmarker(?!_)",
                           integrationMethod = "harmony|\\bcca\\b|\\brpca\\b|reciprocal pca",
                           SCT.regress.CellCycle = "regress\\w*[^.;]*cell.cycle|cell.cycle[^.;]*regress"),
  ## relabelling only; annotation steps
  EzAppScSeuratCombinedLabelClusters = c(computePathwayTFActivity = "decoupler|dorothea|progeny",
                                         SingleR = "singler", enrichrDatabase = "enrichr",
                                         tissue = "aucell|(?<!_)cellmarker(?!_)"),
  ## SpotClean, Azimuth and Enrichr
  EzAppSpatialSeurat = c(spotClean = "spotclean", Azimuth = "(?<!_)azimuth(?![_\\w]|[ _]cell[ _]types)", enrichrDatabase = "enrichr",
                         SCT.regress.CellCycle = "regress\\w*[^.;]*cell.cycle|cell.cycle[^.;]*regress"),
  ## batch correction by CCA
  EzAppSpatialSeuratSlides = c(batchCorrection = "\\bcca\\b|canonical correlation|integration anchors|findintegrationanchors|integratedata|integratelayers|batch.correct",
                               SCT.regress.CellCycle = "regress\\w*[^.;]*cell.cycle|cell.cycle[^.;]*regress"),
  ## SPLIT, co-occurrence FDR, RCTD
  EzAppXeniumSeurat = c(doSPLIT = "spatial purification|split-purified|\\bsplit (\\(|version|algorithm|method|spatial|purif)|(using|with|after|via) (the )?split\\b",
                        coocFdr = "co-?occurrence[^.;]*benjamini|benjamini[^.;]*co-?occurrence",
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
  EzAppSTAR = c(twopassMode = "two-pass|twopassmode\\W+basic", barcodePattern = "umi_tools|umi-tools|umitools"),
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
## raw with every Description sentence removed that states one of `numbers`, names one of the
## step regexes `steps` without negating it, or states one of the `resources` phrases.
## Headings and URL (reference) lines are left alone.
methodsDropSentences <- function(raw, numbers = character(0), steps = character(0), resources = character(0)) {
  esc <- function(t) gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", t, perl = TRUE)
  numRe <- if (length(numbers)) paste0("(?<![0-9.,])(", paste(esc(numbers), collapse = "|"), ")(?![0-9])")
  bad <- function(x) {
    lx <- tolower(x)
    (!is.null(numRe) && grepl(numRe, x, perl = TRUE)) ||
      any(vapply(tolower(resources), grepl, logical(1), x = lx, fixed = TRUE)) ||
      (any(vapply(steps, grepl, logical(1), x = lx, perl = TRUE)) && !grepl(METHODS_NEGATION, lx, perl = TRUE))
  }
  lines <- strsplit(raw, "\n", fixed = TRUE)[[1]]
  prose <- !grepl("^#|https?://", lines) & nzchar(trimws(lines))
  lines[prose] <- vapply(lines[prose], function(l) {
    sents <- strsplit(l, "(?<=[.;])\\s+", perl = TRUE)[[1]]
    paste(sents[!vapply(sents, bad, logical(1))], collapse = " ")
  }, "", USE.NAMES = FALSE)
  paste(lines, collapse = "\n")
}

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

## Compute resources, job modes and run outcome, which the task already forbids and the
## writer still wrote into about 1 text in 10 ("four cores and 12 GB of RAM", "8 threads",
## "dataset mode"). Numbers may be spelled out (one to sixty-four) and have a word before
## the unit ("8 local cores"); memory may come before its amount ("memory was set to
## 100 GB"), in MB or M; thread flags may be quoted ("-@ 4", "--threads 8").
METHODS_NUMBER_WORDS <- paste0("(?:(?:twenty|thirty|forty|fifty|sixty)(?:[- ](?:one|two|three|four|five|six|seven|eight|nine))?|",
                               "ten|eleven|twelve|thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|nineteen|",
                               "one|two|three|four|five|six|seven|eight|nine)")
METHODS_RESOURCE_PATTERN <- paste0(
  "\\b(?:\\d+|", METHODS_NUMBER_WORDS, ")(?:[\\s-]+[a-z]+)?[\\s-]*(?:cores?|cpus?|threads?|workers)\\b|",
  "\\b\\d+(?:\\.\\d+)?\\s*[gm]i?b?\\b[^.;]{0,20}\\b(?:ram|memory)\\b|",
  "\\b(?:ram|memory)\\b[^.;]{0,30}?\\b\\d+(?:\\.\\d+)?\\s*[gm]i?b?\\b|",
  "(?<!\\S)-@\\s*\\d+|--(?:threads|runthreadn|num[-_]?threads|nthreads|cores|cpus|local(?:mem|cores))\\b|",
  "\\b(?:completed|finished) successfully\\b|\\bran to completion\\b|",
  "\\bscratch\\b|\\b(?:dataset|sample) mode\\b|process_mode")
methods_check_resources <- function(description) {
  text <- tolower(paste(description, collapse = "\n"))
  unique(regmatches(text, gregexpr(METHODS_RESOURCE_PATTERN, text, perl = TRUE))[[1]])
}

## What the prompt forbids and the writer still wrote in round 5 (3 of 12 texts named ezRun,
## 1 a file): the framework, file names, paths, command-line options, and a [not recorded]
## that covers "any additional" settings rather than one a reader needs.
METHODS_STYLE_PATTERN <- paste0(
  "\\bezrun\\b|\\bsushi\\b|",
  "(?<![\\w.])[\\w-]+\\.(?:tsv|csv|h5|h5ad|qs2|rds|html|bam|cram|fastq|fq|gz|txt|json|mtx|loom|bed|gtf|fa|fasta)\\b|",
  "(?<![\\w:/])/[\\w.-]+/[\\w./-]+|(?<![\\w-])--[a-z][\\w-]+|",
  "\\[not recorded[^]]*\\b(?:any (?:other|additional|further)|beyond)\\b[^]]*\\]|\\[not recorded\\][^.;]*\\bbeyond\\b|",
  "\\b(?:any (?:other|additional|further)|beyond)\\b[^.;]*\\[not recorded\\]")
methods_check_style <- function(description) {
  text <- tolower(paste(description, collapse = "\n"))
  unique(regmatches(text, gregexpr(METHODS_STYLE_PATTERN, text, perl = TRUE))[[1]])
}

## 10x martian pipelines (Cell Ranger, Space Ranger) log each stage "(ready)" when scheduled
## and "(run:local)" when it runs; a disabled stage is scheduled and never run. The writer
## described vendor steps from memory (hierarchical clustering, cell typing that never ran).
## Steps a reader would name: the stage that produces each (it ran when that stage ran in any
## subpipeline) and the words that describe it.
METHODS_MARTIAN_STEPS <- list(
  "PCA" = c("RUN_PCA", "\\bpca\\b|principal component"),
  "UMAP" = c("RUN_UMAP", "\\bumap\\b"),
  "t-SNE" = c("RUN_TSNE", "\\bt-?sne\\b"),
  "graph-based clustering" = c("RUN_GRAPH_CLUSTERING", "graph-based clustering"),
  "k-means clustering" = c("RUN_KMEANS", "k-means"),
  "hierarchical clustering" = c("RUN_HIERARCHICAL_CLUSTERING", "hierarchical clustering"),
  "differential expression between clusters" = c("RUN_DIFFERENTIAL_EXPRESSION", "differential expression[^.;]*clusters"),
  "cell type annotation" = c("WRITE_CELL_TYPES_H5", "cell[- ]typing|cell[- ]type (?:annotation|assignment|calls?)|(?:annotated|assigned)[^.;]*cell types"),
  "differential expression between cell types" = c("TIDY_CELLTYPE_DIFFEXP", "differential[^.;]*between (?:the )?cell types"),
  "clonotype grouping" = c("RUN_ENCLONE", "enclone|clonotype"),
  "chemistry batch correction" = c("CORRECT_CHEMISTRY_BATCH", "chemistry batch|batch[- ]correct"))

## list(ran, skipped) of METHODS_MARTIAN_STEPS names, from the _o.logs; NULL without martian lines.
methodsMartianStages <- function(o_logs) {
  l <- unlist(lapply(Filter(file.exists, as.character(o_logs)), readLines, warn = FALSE))
  m <- regmatches(l, regexec("\\[runtime\\] \\((ready|run:\\w+)\\)\\s+(\\S+)", l))
  m <- do.call(rbind, m[lengths(m) == 3])
  if (is.null(m)) return(NULL)
  stage <- sub("^.*\\.", "", sub("\\.fork.*$", "", m[, 3]))
  ran <- unique(stage[m[, 2] != "ready"])
  st <- vapply(METHODS_MARTIAN_STEPS, `[`, "", 1)
  list(ran = names(st)[st %in% ran], skipped = names(st)[st %in% setdiff(stage, ran)])
}

## Vendor steps the Description names although the pipeline did not run them.
methods_check_vendor <- function(description, stages) {
  if (!length(stages$skipped)) return(character(0))
  text <- tolower(paste(description, collapse = "\n"))
  Filter(function(s) methodsClaims(text, METHODS_MARTIAN_STEPS[[s]][2]), stages$skipped)
}

## The samples a single (DATASET) job ran on: the "samples" parameter, else the input rows;
## for a two-group comparison only the rows of the compared groups (and baselines), which is
## all a DE fit uses (a review found "all 43 samples" for a DESeq2 fit on 25).
methodsSamplesRun <- function(user_param, input = NULL) {
  s <- user_param$samples
  rows <- if (ezIsSpecified(s)) strsplit(s, ",")[[1]] else input$Name
  groups <- unlist(user_param[c("sampleGroup", "refGroup", "sampleGroupBaseline", "refGroupBaseline")])
  groups <- groups[nzchar(groups)]
  col <- intersect(c(user_param$grouping, paste(user_param$grouping, "[Factor]")), names(input))
  if (length(groups) >= 2 && length(col) && "Name" %in% names(input)) {
    inGroups <- input$Name[input[[col[1]]] %in% groups]
    return(length(intersect(if (length(rows)) rows else input$Name, inGroups)))
  }
  if (length(rows)) length(rows) else NROW(input)
}

## Facts about the run itself, read from its record rather than from the code, so they hold
## for any ezRun version: how the samples were run (a review found "applied to all 12
## samples" for one DATASET job), organism and reference, and the vendor steps that ran.
methodsRunSummary <- function(user_param, input = NULL, sample_count = 1, o_logs = character(0)) {
  n <- sample_count %||% 1
  mode <- user_param$process_mode
  rb <- strsplit(user_param$refBuild %||% "", "/", fixed = TRUE)[[1]]
  sp <- if ("Species" %in% names(input)) unique(input$Species[nzchar(input$Species)])
  stages <- methodsMartianStages(o_logs)
  ## DATASET: no count, an app may use a subset of the input rows (a DE fit keeps two groups;
  ## a review found "all 43 samples" for a fit on 25)
  c(if (identical(mode, "DATASET")) "The samples were analysed together, in one job."
    else if (identical(mode, "SAMPLE") && n > 1)
      paste0("Each of the ", n, " samples was analysed separately, in its own job, with the same settings."),
    if (length(rb) >= 3)
      paste0("Reference: ", gsub("_", " ", rb[1]), ", ", rb[2], " ", rb[3],
             if (length(rb) >= 5 && grepl("^Release_", rb[5])) paste0(", annotation release ", sub("^Release_([^-]+).*$", "\\1", rb[5])), ".")
    else if (length(sp) == 1) paste0("Organism: ", sp, "."),
    if (length(stages$ran)) paste0("The vendor pipeline ran these steps: ", paste(stages$ran, collapse = ", "), "."),
    if (length(stages$skipped)) paste0("It did not run these steps, so do not describe them: ", paste(stages$skipped, collapse = ", "), "."))
}

## Keys left out of the template's parameter list: scheduler, bookkeeping and credentials.
METHODS_TEMPLATE_SKIP_PARAMS <- c("cores", "ram", "scratch", "partition", "process_mode", "samples",
                                  "name", "mail", "adminMail", "sushi_app", "sushiApp", "specialOptions",
                                  "dataRoot", "resultDir", "Rversion", "isLastJob", "inputDatasetName",
                                  "projectId", "node", "nodes", "gpu_feature", "appName")
METHODS_TEMPLATE_DECLARATION <- paste(
  "This description is template text assembled by ezRun from the run's parameters and",
  "the app's recorded behaviour; no language model was used.")

## The untyped key/value pairs of <resultDir>/parameters.tsv, i.e. what the job was given.
methodsParamTable <- function(param_file) {
  if (is.null(param_file) || !file.exists(param_file)) return(list())
  tab <- utils::read.delim(param_file, header = FALSE, colClasses = "character",
                           quote = "", comment.char = "")
  stats::setNames(as.list(sub('^"(.*)"$', "\\1", tab[[2]])), tab[[1]])
}

## Fallback Description when the writer fails or its text fails the guards twice:
## app, ezRun version (from the logs' ezRun_x.y.z), sample count, parameters, app facts.
methods_template <- function(class_name, param, facts, citations, sample_count, log_paths) {
  ver <- character(0)
  for (f in Filter(file.exists, as.character(log_paths))) {
    l <- readLines(f, warn = FALSE)
    ver <- regmatches(l, regexpr("ezRun_\\d+(\\.\\d+)+", l))
    if (length(ver)) break
  }
  keep <- vapply(param, function(v) is.atomic(v) && any(nzchar(as.character(v))), logical(1)) &
    !(names(param) %in% METHODS_TEMPLATE_SKIP_PARAMS) &
    !grepl("apikey|password|token|secret", names(param), ignore.case = TRUE)
  kv <- paste(names(param)[keep], vapply(param[keep], paste, "", collapse = ", "), sep = " = ")
  description <- paste0("The analysis was run with the FGCZ SUSHI app ", sub("^EzApp(.)", "\\1", class_name),
                        if (length(ver)) paste0(" (ezRun ", sub("ezRun_", "", ver[1]), ")"),
                        if (isTRUE(sample_count > 1)) paste0(" on ", sample_count, " samples"), ".",
                        if (length(kv)) paste0(" Parameters: ", paste(kv, collapse = "; "), "."))
  if (length(facts)) description <- paste0(description, "\n\n", paste(facts, collapse = " "))
  list(description = description,
       references = if (length(citations)) paste(citations, collapse = "\n") else "pending")
}

## The Description of a writer response. The model may put its References block first or
## last; it echoes candidate entries, each ending in a URL, and the prose carries no URLs.
## So with a "## References" header the Description is every other non-URL line.
methodsDescriptionPart <- function(raw) {
  lines <- strsplit(raw, "\n", fixed = TRUE)[[1]]
  if (!any(grepl("^## References", lines))) return(raw)
  ## the header, URL (reference) lines and a bare "None" the model writes for no references
  drop <- grepl("^## References", lines) | grepl("https?://", lines) | grepl("^\\s*(none|n/a)\\.?\\s*$", lines, ignore.case = TRUE)
  trimws(paste(lines[!drop], collapse = "\n"))
}

methodsParamText <- function(param) {
  vapply(param, function(v) tryCatch(paste(as.character(unlist(v)), collapse = " "),
                                     error = function(e) ""), "")
}

## write_methods() calls this instead of methods_description(): the writer's text if it
## passes the number and off-step guards (after at most one retry told what failed),
## otherwise list(template = methods_template(...)). A methods_description() override
## without extra_task (static text) is returned unchecked.
methodsGuardedWrite <- function(app, script_paths, log_paths, sample_count, output_dir,
                                param = list(), user_param = list(), use_facts = TRUE) {
  if (!"extra_task" %in% names(formals(app$methods_description)))
    return(list(raw = app$methods_description(script_paths, log_paths, sample_count, output_dir, param)))
  cls <- class(app)[1]
  facts <- if (use_facts) app$methods_facts(param) else character(0)
  citations <- app$methods_citations(param)
  readAll <- function(paths) unlist(lapply(Filter(file.exists, as.character(paths)), readLines, warn = FALSE))
  config <- c(readAll(script_paths), methodsParamText(param), methodsParamText(user_param),
              as.character(unlist(app$appDefaults$DefaultValue)), facts, citations)
  all <- c(config, readAll(log_paths))
  stages <- methodsMartianStages(log_paths[grepl("_o\\.log$", log_paths)])
  template <- function(reason) {
    message("write_methods: template fallback (", reason, ")")
    list(template = methods_template(cls, user_param, facts, citations, sample_count, log_paths))
  }
  extra <- NULL
  for (attempt in 1:2) {
    raw <- tryCatch(app$methods_description(script_paths, log_paths, sample_count, output_dir,
                                            param, use_facts = use_facts, extra_task = extra),
                    error = function(e) e)
    if (inherits(raw, "error")) return(template(conditionMessage(raw)))
    description <- methodsDescriptionPart(raw)
    ## headings alone (the model's "## Methods") are no Description either
    if (!nzchar(trimws(gsub("(^|\n)#+[^\n]*", "", description)))) return(template("empty description"))
    ## Plain numbers are checked against the logs too: on 140 archived runs, matching the
    ## configuration only flagged 57 texts, all but one for values the logs do hold (sample
    ## counts, tool defaults on command lines, reference sizes), and the retry deleted them.
    numbers <- methods_check_numbers(description, all, all, sample_count)
    steps <- methods_check_offsteps(description, cls, param)
    resources <- methods_check_resources(description)
    vendor <- methods_check_vendor(description, stages)
    style <- methods_check_style(description)
    ## style is asked out once, not dropped (except a sweeping [not recorded]): its sentences carry the tool and its version
    ## ("DESeq2 1.52.0 within the ezRun framework"), so it is left when the retry keeps it
    ## a sweeping [not recorded] sentence says nothing, so it is dropped like the others
    notrec <- grep("not recorded", style, value = TRUE)
    if (!length(numbers) && !length(steps) && !length(resources) && !length(vendor) && !length(notrec) &&
        (!length(style) || attempt == 2))
      return(list(raw = raw))
    message("write_methods: guards flagged ", length(numbers), " numbers, ", length(steps), " steps, ",
            length(resources), " resources, ", length(vendor), " vendor steps, ", length(style), " style",
            if (attempt == 1) "; retrying")
    extra <- paste(c(
      if (length(numbers)) paste0("These values are not in the run's configuration; remove them or the ",
                                  "sentence that states them: ", paste(numbers, collapse = ", "), "."),
      if (length(steps)) paste0("These steps were not run in this job; do not describe them: ",
                                paste(sub("^(.*):(.*)$", "\\2 (\\1 off)", steps), collapse = ", "), "."),
      if (length(resources)) paste0("Remove the compute resources, job settings and statements that the ",
                                    "run completed, which are not part of the method: ", paste(resources, collapse = ", "), "."),
      if (length(vendor)) paste0("The vendor pipeline did not run these steps in this job; do not describe them: ",
                                 paste(vendor, collapse = ", "), "."),
      if (length(style)) paste0("Remove the framework names (ezRun, SUSHI), file names, paths and command-line ",
                                "options, and any [not recorded] that is not about one specific setting: ",
                                paste(style, collapse = ", "), ".")),
      collapse = "\n")
  }
  ## Still flagged after the retry: drop just those sentences, so one bad sentence does not
  ## cost the whole text; the template only when nothing is left.
  kept <- methodsDropSentences(raw, numbers, c(METHODS_OFFSTEP_RULES[[cls]][sub(":.*$", "", steps)],
                                               vapply(METHODS_MARTIAN_STEPS[vendor], `[`, "", 2)),
                               c(resources, notrec))
  if (nzchar(trimws(gsub("(^|\n)#+[^\n]*", "", methodsDescriptionPart(kept))))) {
    message(sprintf("write_methods: dropped the sentences with %d numbers, %d steps, %d resources, %d vendor steps",
                    length(numbers), length(steps), length(resources), length(vendor)))
    return(list(raw = kept))
  }
  template(sprintf("guards: %d numbers, %d steps", length(numbers), length(steps)))
}

## CellRanger(-Multi) writes <resultDir>/<sample>/config.csv. SUSHI names a SAMPLE-mode job
## script <category>_<sample>.sh or <category>_<sample>_<dataset>_<id>.sh, so the example
## sample's config is the one whose directory name sits between underscores in the script
## name; the longest such name wins ("CTCL_PBMCs3" over "PBMCs3"). No match: the first one
## only. Without an example script (DATASET mode) every config.csv is read, as the scripts are.
methodsConfigCsv <- function(result_dir, example_script = NULL) {
  csv <- sort(Sys.glob(file.path(result_dir, "*", "config.csv")))
  if (is.null(example_script) || length(csv) <= 1) return(csv)
  samples <- basename(dirname(csv))
  script <- paste0("_", sub("\\.sh$", "_", basename(example_script)))
  hit <- vapply(samples, function(s) grepl(paste0("_", s, "_"), script, fixed = TRUE), logical(1))
  if (!any(hit)) return(csv[1])
  csv[hit][which.max(nchar(samples[hit]))]
}

## The _e.log of each job's latest attempt. SUSHI names a job log
## <script>.sh_sushiID<n>_<YYYY-MM-DD--HH-MM-SS>_e.log (before the sushiID naming:
## <script>.sh_<stamp>_e.log); a resubmitted job keeps its script and sushiID and gets a
## newer stamp, so counting log files counted a job that was rerun successfully as stopped.
methodsLatestJobLogs <- function(logs) {
  if (!length(logs)) return(character(0))
  b <- basename(logs)
  script <- sub("\\.sh_.*$", ".sh", b)
  stamp <- sub("^sushiID\\d+_", "", sub("_e\\.log$", "", substring(b, nchar(script) + 2)))
  o <- order(script, stamp, decreasing = TRUE)
  logs[sort(o[!duplicated(script[o])])]
}

## SLURM's own lines in a job's _e.log when it killed the job: a cancel (time limit, scancel)
## or an out-of-memory kill. Neither leaves "Execution halted".
METHODS_SLURM_KILL <- "\\*\\*\\* (JOB|STEP) \\S+ ON \\S+ CANCELLED AT|Detected \\d+ oom[-_]kill event"

## A job failed when R stopped ("Execution halted"; an "Error in" line alone can be a caught
## error, SoupX autoEstCont in a ScSeurat run that delivered), when SLURM killed it, or when
## its traced script did not reach its end. SUSHI job scripts run under `set -eux` (the
## trace opens with "+ umask 0002"), so the footer's last command, "+ rm -rf <scratch dir>",
## is traced only when every command before it succeeded: on p2220 + p28409 it ends all 535
## traced _e.logs without a halt or a kill message and none of the 18 others (13 stopped
## in the g-req copy, which does not count as failed, see below). Untraced (older) logs have
## no such marker, so only the first two apply.
methodsJobFailed <- function(e_log) {
  l <- readLines(e_log, warn = FALSE)
  if (any(grepl("^Execution halted", l)) || any(grepl(METHODS_SLURM_KILL, l))) return(TRUE)
  traced <- any(grepl("^\\+ umask ", utils::head(l, 5)))
  ## The footer copies the results with g-req before that rm -rf. A trace that reached the copy
  ## means the analysis itself finished; a failed copy (SUSHI marks the job FAILED, routine for
  ## Cell Ranger with complete outputs) is a delivery problem, not a run that did not complete.
  traced && !any(grepl("^\\+ rm -rf ", utils::tail(l[nzchar(trimws(l))], 3))) && !any(grepl("^\\+ g-req ", l))
}

## Why the jobs of these _e.logs failed, for the failed-run statement and note: SLURM's
## kill, else the first "Error" line (with the next line when it ends in a colon, as
## "Error in ezSystem(cmd) :" does), with paths reduced to file names.
methodsJobError <- function(e_logs) {
  l <- unlist(lapply(e_logs, readLines, warn = FALSE))
  k <- grep(METHODS_SLURM_KILL, l, value = TRUE)[1]
  i <- grep("^Error", l)[1]
  e <- if (!is.na(k)) {
    if (grepl("oom", k)) "SLURM killed the job for exceeding its memory" else
      paste0("SLURM cancelled the job", if (grepl("DUE TO", k)) paste0(" due to ", tolower(sub(".*DUE TO (.*?) *\\*+.*$", "\\1", k))))
  } else if (!is.na(i)) {
    paste(trimws(l[i:(i + (grepl(":\\s*$", l[i]) && i < length(l)))]), collapse = " ")
  } else if (any(grepl("^Execution halted", l))) "Execution halted" else "the job ended before its last step"
  substr(gsub("(/[^ /]+)+/([^ /]*)", "\\2", e), 1, 200)   # no paths in a Methods file
}
