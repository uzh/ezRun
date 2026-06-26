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
  ## 1. Stage Bracken reports under deterministic filenames so the
  ##    short SampleIDs that kraken-biom derives match input$getNames().
  ## ---------------------------------------------------------------
  reportsDir <- "bracken_reports"
  dir.create(reportsDir, showWarnings = FALSE, recursive = TRUE)
  sampleNames    <- input$getNames()
  brackenReports <- input$getFullPaths("BrackenReport")
  for (i in seq_along(sampleNames)) {
    file.copy(
      brackenReports[i],
      file.path(reportsDir, paste0(sampleNames[i], ".bracken.report.txt")),
      overwrite = TRUE
    )
  }

  ## ---------------------------------------------------------------
  ## 2. Build sample-metadata.tsv from the Factor-tagged columns only.
  ##    File/Link/B-Fabric columns are intentionally dropped.
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

  ## ---------------------------------------------------------------
  ## 3. Build the BIOM table (HDF5) from the staged Bracken reports.
  ## ---------------------------------------------------------------
  biomFile <- "bracken_table.biom"
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

  ## ---------------------------------------------------------------
  ## 4. Normalise BIOM sample IDs (strip suffixes after first '.') and
  ##    embed the clean metadata via the biom Python API.
  ## ---------------------------------------------------------------
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
##' @description Differential abundance analysis on Bracken/Kraken2 shotgun
##' metagenomics profiles. Builds a BIOM table from per-sample Bracken
##' reports, then runs a single DA method selected by `param$daMethod`
##' (ALDEx2, ANCOM-BC2, MaAsLin3, or LEfSe). Each of the four SUSHI apps
##' (DiffShotALDEx, DiffShotANCOMBC, DiffShotMaAsLin3, DiffShotLEfSe)
##' hardcodes this parameter and dispatches to a per-method child Rmd.
##' Produces a single self-contained HTML report (00index.html) with a
##' volcano + table tabset for the chosen method.
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
