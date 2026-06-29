###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.character(a) && length(a) == 1 && nchar(a) == 0) return(b)
  a
}

## ------------------------------------------------------------------
## Read one MetaPhlAn-4 profile file. Returns a data.frame whose
## column names come from the last #-prefixed header line (e.g.
## "#clade_name<TAB>clade_taxid<TAB>relative_abundance<TAB>coverage<TAB>
##  estimated_number_of_reads_from_the_clade").
## Older / default-mode profiles (no `-t rel_ab_w_read_stats`) carry
## "#clade_name<TAB>NCBI_tax_id<TAB>relative_abundance<TAB>
##  additional_species" — the `estimated_number_of_reads_from_the_clade`
## column will simply be absent, which the caller handles.
## ------------------------------------------------------------------
.read_metaphlan_profile <- function(path) {
  con <- file(path, "r"); on.exit(close(con))
  lastComment <- NULL
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#")) {
      lastComment <- line
    } else {
      break
    }
  }
  if (is.null(lastComment))
    stop("MetaPhlAn profile '", path, "' has no header line.")
  cols <- strsplit(sub("^#", "", lastComment), "\t", fixed = TRUE)[[1]]
  cols <- trimws(cols)
  df <- tryCatch(
    read.delim(path, comment.char = "#", header = FALSE,
               stringsAsFactors = FALSE, check.names = FALSE,
               quote = "", na.strings = c("", "NA")),
    error = function(e) data.frame()
  )
  if (nrow(df) == 0) {
    return(setNames(data.frame(matrix(ncol = length(cols), nrow = 0)),
                    cols))
  }
  if (ncol(df) < length(cols))
    df[, (ncol(df) + 1):length(cols)] <- NA
  if (ncol(df) > length(cols))
    df <- df[, seq_along(cols), drop = FALSE]
  colnames(df) <- cols
  df
}

## ------------------------------------------------------------------
## Parse a MetaPhlAn lineage string like
##   "k__Bacteria|p__Firmicutes|c__Clostridia|...|s__Faecalibacterium_prausnitzii"
## into the 7-rank vector that observation_metadata() will surface.
## Returns a named character vector with NA for missing ranks.
## ------------------------------------------------------------------
.parse_metaphlan_lineage <- function(lineage) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order",
             "Family", "Genus", "Species")
  prefMap <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order",
               f = "Family", g = "Genus", s = "Species")
  out <- setNames(rep(NA_character_, 7), ranks)
  if (is.na(lineage) || !nzchar(lineage)) return(out)
  for (tok in strsplit(lineage, "|", fixed = TRUE)[[1]]) {
    pre <- substr(tok, 1, 1)
    if (substr(tok, 2, 3) != "__") next
    rk <- prefMap[pre]
    if (is.na(rk)) next
    out[rk] <- sub("^[a-z]__", "", tok)
  }
  out
}

## ------------------------------------------------------------------
## Build the bracken_table.biom equivalent from MetaPhlAn profiles.
## - profilePaths    : character vector of file paths
## - sampleNames     : same length, in the same order
## - needsCounts     : if TRUE, require estimated_number_of_reads_from_the_clade
##                     and abort with a clear message otherwise; values used
##                     as-is (already integer-ish). If FALSE, use the
##                     relative_abundance column scaled to integer pseudo-
##                     counts (round(relab * 1e4)).
## - metadataFile    : sample-metadata.tsv to embed into the BIOM
## - outBiomFile     : output path
## - methodLabelStr  : included in the abort message for count methods
## Only species-resolved rows (lineage ending in "|s__...") are kept; the
## Rmd's tax_glom then rolls up to Genus, mirroring how kraken-biom's
## --min S --max D table behaves.
## ------------------------------------------------------------------
build_biom_from_metaphlan <- function(profilePaths, sampleNames,
                                      needsCounts, metadataFile, outBiomFile,
                                      methodLabelStr = "") {
  if (!requireNamespace("biomformat", quietly = TRUE))
    stop("Package 'biomformat' is required to load MetaPhlAn profiles. ",
         "Install it from Bioconductor.")

  countsCol  <- "estimated_number_of_reads_from_the_clade"
  relabCol   <- "relative_abundance"
  valueCol   <- if (needsCounts) countsCol else relabCol

  perSample <- lapply(profilePaths, .read_metaphlan_profile)

  if (needsCounts) {
    missingCols <- vapply(perSample, function(d) !(countsCol %in% colnames(d)),
                          logical(1))
    if (any(missingCols)) {
      bad <- sampleNames[missingCols]
      stop("DiffShot ", methodLabelStr,
           " requires integer counts, but the following MetaPhlAn profile(s) ",
           "do NOT carry the '", countsCol, "' column: ",
           paste(bad, collapse = ", "),
           ".\nThis means MetaPhlAn was run WITHOUT '-t rel_ab_w_read_stats'. ",
           "Re-run the MetaPhlAn SUSHI app with 'Estimation of read counts ",
           "mapped to clade' = true (the default), or pick the relative-",
           "abundance-based DiffShot apps (DiffShotMaAsLin3 or DiffShotLEfSe) ",
           "instead, which accept any MetaPhlAn profile.")
    }
  }

  ## Keep only species-resolved rows. MetaPhlAn writes the full lineage in
  ## clade_name; the strain rank "|t__..." is rare (only some DBs) and is
  ## treated as below-species — drop those too so species counts don't
  ## double-count their own strains.
  isSpeciesRow <- function(name)
    grepl("\\|s__[^|]+$", name, perl = FALSE)

  perSampleVals <- lapply(seq_along(perSample), function(i) {
    d <- perSample[[i]]
    if (!("clade_name" %in% colnames(d)) || !(valueCol %in% colnames(d)))
      return(setNames(numeric(0), character(0)))
    keep <- isSpeciesRow(d[["clade_name"]]) &
            !is.na(suppressWarnings(as.numeric(d[[valueCol]])))
    setNames(as.numeric(d[[valueCol]][keep]), d[["clade_name"]][keep])
  })

  allTaxa <- sort(unique(unlist(lapply(perSampleVals, names),
                                use.names = FALSE)))
  if (length(allTaxa) == 0)
    stop("No species-level rows found across the supplied MetaPhlAn profiles.")

  mat <- matrix(0, nrow = length(allTaxa), ncol = length(sampleNames),
                dimnames = list(allTaxa, sampleNames))
  for (i in seq_along(perSampleVals)) {
    v <- perSampleVals[[i]]
    if (length(v) > 0)
      mat[names(v), sampleNames[i]] <- v
  }

  ## Scale relab (0-100 %) to integer pseudo-counts so the count-shaped
  ## filters (sample_sums, > 1 thresholds, CPM / TSS normalisation) keep
  ## working unchanged. relab * 1e4 -> sums on the order of 1e6 per sample.
  if (!needsCounts) {
    mat <- round(mat * 1e4)
  } else {
    mat <- round(mat)
  }
  storage.mode(mat) <- "integer"

  ## 7-rank taxonomy
  taxList <- lapply(allTaxa, .parse_metaphlan_lineage)
  taxMat  <- do.call(rbind, taxList)
  rownames(taxMat) <- allTaxa
  taxDf <- as.data.frame(taxMat, stringsAsFactors = FALSE)

  ## Build the BIOM and write to HDF5 (same shape as kraken-biom output).
  biomObj <- biomformat::make_biom(
    data                 = mat,
    observation_metadata = taxDf,
    sample_metadata      = if (file.exists(metadataFile)) {
                              md <- read.delim(metadataFile, sep = "\t",
                                               header = TRUE,
                                               stringsAsFactors = FALSE,
                                               check.names = FALSE)
                              rownames(md) <- md[[1]]
                              md <- md[, -1, drop = FALSE]
                              md[sampleNames, , drop = FALSE]
                            } else NULL
  )
  biomformat::write_biom(biomObj, outBiomFile)
  invisible(biomObj)
}

ezMethodDiffShot <- function(input = NA, output = NA, param = NA) {
  require(withr)
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  defer(setwd(cwd))

  ## The conda env (`gi_qiime2-amplicon-2026.4`) is already active because
  ## the SUSHI Ruby `commands` block sources it before launching R. That
  ## env provides R itself plus kraken-biom / biom / h5py / python3, so
  ## downstream ezSystem() calls find them on PATH.

  ## ---------------------------------------------------------------
  ## 0. Decide input profile type (Bracken vs MetaPhlAn) from the
  ##    dataset columns. The four DiffShot SUSHI apps advertise an
  ##    XOR-of-AND on @required_columns, so exactly one of the two
  ##    will be present.
  ## ---------------------------------------------------------------
  sampleNames <- input$getNames()
  hasBracken    <- "BrackenReport"   %in% input$colNames
  hasMetaPhlAn  <- "MetaPhlAnProfile" %in% input$colNames
  if (!hasBracken && !hasMetaPhlAn)
    stop("Input dataset has neither a BrackenReport nor a MetaPhlAnProfile column.")
  if (hasBracken && hasMetaPhlAn)
    stop("Input dataset carries both BrackenReport and MetaPhlAnProfile columns; ",
         "ambiguous - keep only one.")
  sourceType <- if (hasBracken) "bracken" else "metaphlan"
  param$sourceType <- sourceType   # exposed to the Rmd via param

  ## ---------------------------------------------------------------
  ## 1. Build sample-metadata.tsv from the Factor-tagged columns only.
  ##    File/Link/B-Fabric columns are intentionally dropped. Shared
  ##    by both input paths.
  ## ---------------------------------------------------------------
  factorMask <- input$columnHasTag("Factor")
  if (!any(factorMask)) {
    stop("Input dataset has no [Factor] columns - differential abundance ",
         "needs at least one (the grouping variable).")
  }
  factorMeta <- input$meta[, factorMask, drop = FALSE]
  colnames(factorMeta) <- input$colNames[factorMask]
  cleanMeta <- data.frame(
    SampleID = sampleNames,
    factorMeta,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadataFile <- "sample-metadata.tsv"
  write.table(cleanMeta, metadataFile, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")

  biomFile <- "bracken_table.biom"   # kept regardless of source so the Rmd
                                     # path doesn't have to switch on it.

  if (sourceType == "bracken") {

    ## -------------------------------------------------------------
    ## 2a. Stage Bracken reports under deterministic filenames so the
    ##     short SampleIDs that kraken-biom derives match input$getNames().
    ## -------------------------------------------------------------
    reportsDir <- "bracken_reports"
    dir.create(reportsDir, showWarnings = FALSE, recursive = TRUE)
    brackenReports <- input$getFullPaths("BrackenReport")
    for (i in seq_along(sampleNames)) {
      file.copy(
        brackenReports[i],
        file.path(reportsDir, paste0(sampleNames[i], ".bracken.report.txt")),
        overwrite = TRUE
      )
    }

    ## -------------------------------------------------------------
    ## 3a. Build the BIOM table (HDF5) from the staged Bracken reports.
    ## -------------------------------------------------------------
    maxRank  <- if (!is.null(param$brackenMaxRank)) param$brackenMaxRank else "D"
    minRank  <- if (!is.null(param$brackenMinRank)) param$brackenMinRank else "S"
    reportFiles <- list.files(reportsDir, full.names = TRUE)
    ezSystem(paste(
      "kraken-biom",
      paste(shQuote(reportFiles), collapse = " "),
      "--max", shQuote(maxRank),
      "--min", shQuote(minRank),
      "--fmt hdf5",
      "-o", shQuote(biomFile),
      "--verbose"
    ))

    ## -------------------------------------------------------------
    ## 4a. Normalise BIOM sample IDs (strip suffixes after first '.') and
    ##     embed the clean metadata via the biom Python API.
    ## -------------------------------------------------------------
    pyHelper <- "biom_embed_metadata.py"
    writeLines(c(
      "import sys, os, csv",
      "from biom import load_table",
      "from biom.util import biom_open",
      "biom_fp, meta_fp = sys.argv[1], sys.argv[2]",
      "t = load_table(biom_fp)",
      "cur = list(t.ids('sample'))",
      "short = {s: s.split('.')[0] for s in cur}",
      "if len(set(short.values())) != len(cur):",
      "    sys.exit('ERROR: shortening sample IDs produced collisions; check filenames.')",
      "t.update_ids(short, axis='sample', inplace=True)",
      "if os.path.isfile(meta_fp):",
      "    with open(meta_fp, newline='') as fh:",
      "        rows = list(csv.reader(fh, delimiter='\\t'))",
      "    header = rows[0]",
      "    md = {}",
      "    for r in rows[1:]:",
      "        if not r or not r[0].strip():",
      "            continue",
      "        md[r[0].strip()] = {header[i]: (r[i].strip() if i < len(r) else '')",
      "                            for i in range(1, len(header))}",
      "    present = set(t.ids('sample'))",
      "    missing = present - set(md)",
      "    if missing:",
      "        print('WARNING: no metadata row for: ' + ', '.join(sorted(missing)))",
      "    t.add_metadata(md, axis='sample')",
      "with biom_open(biom_fp, 'w') as fh:",
      "    t.to_hdf5(fh, 'EzAppDiffShot')",
      "print('Sample IDs now: ' + ', '.join(list(t.ids('sample'))[:5]) + ' ...')"
    ), pyHelper)
    ezSystem(paste("python3", shQuote(pyHelper),
                   shQuote(biomFile), shQuote(metadataFile)))

  } else {  # sourceType == "metaphlan"

    ## -------------------------------------------------------------
    ## 2b. Build the BIOM table directly from MetaPhlAn profiles.
    ##     For count-based DA methods (ALDEx2 / ANCOM-BC2) we extract
    ##     'estimated_number_of_reads_from_the_clade' and abort if it
    ##     is missing (= MetaPhlAn was run WITHOUT
    ##     'Estimation of read counts mapped to clade'). For
    ##     abundance-based methods (LEfSe / MaAsLin3) we use
    ##     'relative_abundance' (always present); to keep the
    ##     downstream count-shaped pipeline happy the relab values
    ##     are scaled to integer pseudo-counts (relab% * 1e4).
    ## -------------------------------------------------------------
    methodNorm  <- toupper(gsub("[^A-Za-z0-9]", "", as.character(param$daMethod %||% "")))
    needsCounts <- methodNorm %in% c("ALDEX2", "ANCOMBC", "ANCOMBC2")
    profilePaths <- input$getFullPaths("MetaPhlAnProfile")
    build_biom_from_metaphlan(
      profilePaths   = profilePaths,
      sampleNames    = sampleNames,
      needsCounts    = needsCounts,
      metadataFile   = metadataFile,
      outBiomFile    = biomFile,
      methodLabelStr = methodNorm
    )
  }

  ## ---------------------------------------------------------------
  ## 5. Stage the per-method child Rmd into the working dir. The main
  ##    DiffShot.Rmd uses knitr `child=` dispatch on param$daMethod and
  ##    makeRmdReport only copies the top-level template, so the child
  ##    has to be alongside it before render runs.
  ## ---------------------------------------------------------------
  methodNorm <- toupper(gsub("[^A-Za-z0-9]", "", as.character(param$daMethod %||% "")))
  methodChild <- switch(methodNorm,
                        ALDEX2   = "DiffShot_aldex.Rmd",
                        ANCOMBC  = "DiffShot_ancombc.Rmd",
                        ANCOMBC2 = "DiffShot_ancombc.Rmd",
                        MAASLIN3 = "DiffShot_maaslin3.Rmd",
                        MAASLIN  = "DiffShot_maaslin3.Rmd",
                        LEFSE    = "DiffShot_lefse.Rmd",
                        stop("Unknown daMethod '", param$daMethod, "'; expected ",
                             "one of ALDEx2, ANCOMBC, MaAsLin3, LEfSe."))
  childSrc <- system.file("templates", methodChild, package = "ezRun")
  if (!nzchar(childSrc) || !file.exists(childSrc))
    stop("Child template '", methodChild, "' not installed in ezRun. ",
         "Reinstall ezRun into the active R library.")
  file.copy(childSrc, methodChild, overwrite = TRUE)

  ## ---------------------------------------------------------------
  ## 6. Hand control over to the Rmd. Paths are deterministic relative
  ##    to the working dir, so the template loads them by name.
  ## ---------------------------------------------------------------
  methodLabel <- switch(methodNorm,
                        ALDEX2 = "ALDEx2", ANCOMBC = "ANCOM-BC2",
                        ANCOMBC2 = "ANCOM-BC2", MAASLIN3 = "MaAsLin3",
                        MAASLIN = "MaAsLin3", LEFSE = "LEfSe",
                        as.character(param$daMethod))
  comparisonLabel <- if (!is.null(param$comparison) && nzchar(param$comparison))
                       param$comparison
                     else if (!is.null(param$grouping) && nzchar(param$grouping))
                       param$grouping
                     else "Group"
  reportTitle <- paste0("Shotgun DA (", methodLabel, ") - ", comparisonLabel)

  makeRmdReport(
    output      = output,
    param       = param,
    rmdFile     = "DiffShot.Rmd",
    reportTitle = reportTitle
  )

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodDiffShot(input=NA, output=NA, param=NA)
##' @description Differential abundance analysis on shotgun metagenomics
##' taxonomic profiles. Accepts EITHER Bracken/Kraken2 (`BrackenReport`)
##' OR MetaPhlAn (`MetaPhlAnProfile`) input; the back-end builds the BIOM
##' table via `kraken-biom` for Bracken or directly via `biomformat::make_biom`
##' for MetaPhlAn, then runs a single DA method selected by `param$daMethod`
##' (ALDEx2, ANCOM-BC2, MaAsLin3, or LEfSe). Count-based methods (ALDEx2,
##' ANCOM-BC2) with MetaPhlAn input require profiles produced with
##' `-t rel_ab_w_read_stats`; otherwise the run aborts with a clear message.
##' Each of the four SUSHI apps (DiffShotALDEx, DiffShotANCOMBC,
##' DiffShotMaAsLin3, DiffShotLEfSe) hardcodes the method and dispatches to
##' a per-method child Rmd. Produces a single self-contained HTML report
##' (00index.html) with an Abundance section and a volcano + table tabset
##' for the chosen method.
EzAppDiffShot <-
  setRefClass(
    "EzAppDiffShot",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodDiffShot
        name <<- "EzAppDiffShot"
        appDefaults <<- rbind(
          daMethod = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "DA method to run (hardcoded by the parent SUSHI app). One of: ALDEx2, ANCOMBC, MaAsLin3, LEfSe."
          ),
          grouping = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Metadata [Factor] column to test on (categorical)."
          ),
          sampleGroup = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Value of 'grouping' treated as the contrast group; must differ from refGroup."
          ),
          refGroup = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Value of 'grouping' used as the reference level; coefficient signs read sample-vs-ref."
          ),
          daCovariates = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Comma-separated list of additional [Factor] metadata columns to include as covariates. Ignored by the LEfSe app (LEfSe does not accept covariates)."
          ),
          prevalenceMin = ezFrame(
            Type = "numeric",
            DefaultValue = 0.10,
            Description = "Prevalence filter: keep taxa seen (>1 read) in at least this fraction of samples."
          ),
          relabundMin = ezFrame(
            Type = "numeric",
            DefaultValue = 0.10,
            Description = "Relative-abundance floor: keep taxa whose pooled read count is >= this percentage of total."
          ),
          libCut = ezFrame(
            Type = "integer",
            DefaultValue = 1000,
            Description = "ANCOM-BC2 minimum library size (lib_cut). Used only by the ANCOMBC app."
          ),
          pValueThresh = ezFrame(
            Type = "numeric",
            DefaultValue = 0.05,
            Description = "p-value cut-off used on volcano plots."
          ),
          log2RatioThresh = ezFrame(
            Type = "numeric",
            DefaultValue = 1,
            Description = "log2 fold-change cut-off used on volcano plots."
          ),
          maaslinMaxSig = ezFrame(
            Type = "numeric",
            DefaultValue = 0.1,
            Description = "MaAsLin3 max_significance threshold (q-value cutoff for the headline plots). Used only by the MaAsLin3 app."
          ),
          samplesToDrop = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Comma-separated list of SampleID values to exclude entirely (e.g. failed libraries)."
          ),
          brackenMinRank = ezFrame(
            Type = "character",
            DefaultValue = "S",
            Description = "kraken-biom --min: lowest taxonomic rank collapsed to it (D/P/C/O/F/G/S/S1). Should match the rank Bracken was run at."
          ),
          brackenMaxRank = ezFrame(
            Type = "character",
            DefaultValue = "D",
            Description = "kraken-biom --max: highest taxonomic rank captured (D = Domain captures all)."
          ),
          cmdOptions = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Free-form options reserved for future use; currently unused."
          )
        )
      }
    )
  )
