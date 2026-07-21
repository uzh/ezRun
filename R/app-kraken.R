###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodKraken = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  # This worker always loops over input$getNames(): in the classic SAMPLE-mode
  # path `input` carries a single sample (one job per sample); when KrakenApp's
  # `exclusive` option is on it runs in DATASET mode, `input` carries the whole
  # dataset, and the node is reserved with SLURM --exclusive.
  #
  # In the exclusive/dataset case we keep the Kraken2 database resident across
  # samples via the k2 classify daemon: the first `--use-daemon` call loads the DB
  # into memory and starts a background classifier, every later call reuses it, and
  # we stop it with `k2 clean --stop-daemon` at the end — so the (large) DB loads
  # once for the run instead of once per sample. The daemon addresses fixed files
  # in /tmp shared between jobs on a node, so it is only safe on a reserved node:
  # we start it ONLY when `exclusive` is on (and there is more than one sample to
  # benefit). Otherwise each sample runs a plain `k2 classify` (classic behaviour).
  sampleNames <- input$getNames()
  exclusive <- isTRUE(param$exclusive) ||
               identical(tolower(as.character(param$exclusive)), "true")
  useDaemon <- exclusive && length(sampleNames) > 1L
  daemonFlag <- if (useDaemon) "--use-daemon" else ""

  # ---- DB selection: identical for every sample in the dataset ----
  dbVals <- unlist(strsplit(as.character(param$krakenDBOpt), ","))
  dbVals <- trimws(dbVals)
  dbVals <- dbVals[nzchar(dbVals)]
  if (length(dbVals) == 0L) {
    stop("krakenDBOpt is empty: no database selected.")
  }
  multiDB <- isTRUE(param$multiDB) ||
             identical(tolower(as.character(param$multiDB)), "true")
  if (!multiDB && length(dbVals) > 1L) {
    dbVals <- dbVals[1]
  }
  dbPaths <- file.path("/srv/GT/databases/kraken2", dbVals)
  dbArg <- paste(dbPaths, collapse = ",")

  conOpt <- param$krakenConfidenceOpt
  phredOpt <- param$krakenPhredOpt
  minHitGroups <- if (!is.null(param$minimum_hit_groups)) param$minimum_hit_groups else "2"
  wantMinimizer <- isTRUE(param$report_minimizer_data) ||
                   identical(tolower(as.character(param$report_minimizer_data)), "yes")
  if (wantMinimizer && length(dbVals) > 1L) {
    # kraken2 wiki: --report-minimizer-data is unsupported in multi-DB mode.
    message("Disabling --report-minimizer-data: not supported with multiple DBs.")
    wantMinimizer <- FALSE
  }
  reportMinimizerFlag <- if (wantMinimizer) "--report-minimizer-data" else ""

  # --unclassified-out is optional and off by default. In paired mode kraken2
  # requires a '#' placeholder in the filename, which it expands to _1 / _2.
  saveUnclassified <- isTRUE(param$save_unclassified) ||
                     identical(tolower(as.character(param$save_unclassified)), "yes")

  # When the shared daemon is used, always stop it, even if a sample errors out
  # mid-loop, so no orphaned daemon or /tmp/classify.* files are left on the node.
  if (useDaemon) {
    on.exit(try(ezSystem("k2 clean --stop-daemon"), silent = TRUE), add = TRUE)
  }

  for (sampleName in sampleNames) {
    # Trim one sample at a time: ezMethodFastpTrim is single-sample.
    singleInput <- input$subset(sampleName)
    trimmedInput <- ezMethodFastpTrim(input = singleInput, param = param)

    outTxt <- paste0(sampleName, ".txt")
    outReport <- paste0(sampleName, ".report.txt")
    outHtml <- paste0(sampleName, ".html")
    outLog <- paste0(sampleName, ".log")

    if (param$paired) {
      read1 <- trimmedInput$getColumn("Read1")
      read2 <- trimmedInput$getColumn("Read2")
      readOpt <- paste(read1, read2)
      pairedFlag <- "--paired"
    } else {
      read1 <- trimmedInput$getColumn("Read1")
      readOpt <- paste(read1)
      pairedFlag <- ""
    }

    if (saveUnclassified) {
      unclassifiedPattern <- if (param$paired) {
        paste0(sampleName, "_unclassified#.fastq")
      } else {
        paste0(sampleName, "_unclassified.fastq")
      }
      unclassifiedArg <- paste("--unclassified-out", unclassifiedPattern)
    } else {
      unclassifiedArg <- ""
    }

    # daemonFlag is "--use-daemon" only in the exclusive/dataset path: the first
    # sample starts the background classifier and loads the DB into memory, and
    # subsequent samples reuse the in-memory DB. The call blocks until this sample's
    # classification is done, so processing stays sequential. When empty, this is a
    # plain per-sample `k2 classify` (classic behaviour).
    cmd <- paste(
      "k2 classify",
      daemonFlag,
      "--db", dbArg,
      pairedFlag,
      "--use-names",
      "--minimum-hit-groups", minHitGroups,
      reportMinimizerFlag,
      "--confidence", conOpt,
      "--minimum-base-quality", phredOpt,
      "--output", outTxt,
      "--report", outReport,
      unclassifiedArg,
      "--threads", ezThreads(),
      param$cmdOptions,
      readOpt,
      "1>", outLog
    )
    ezSystem(cmd)

    # Per-sample cleanup: the trimmed reads are transient intermediates (already
    # gzipped by fastp) consumed by classify and not surfaced downstream, so drop
    # them now — in DATASET mode every sample shares one scratch dir, and deleting
    # here keeps peak scratch usage to roughly one sample's trimmed reads instead
    # of all samples' at once. The raw upstream FASTQs (linked in the output
    # table) live elsewhere and are untouched.
    trimmedReads <- if (param$paired) c(read1, read2) else read1
    ezSystem(paste("rm -f", paste(shQuote(trimmedReads), collapse = " ")))

    # Krona must run on the uncompressed per-read output
    cmd <- paste(
      "ktImportTaxonomy -q 2 -t 3",
      "-n", shQuote(sampleName),
      outTxt,
      "-o", outHtml
    )
    ezSystem(cmd)

    ezSystem(paste("pigz -f --best -p", ezThreads(), outTxt))
    if (saveUnclassified) {
      unclassifiedFiles <- if (param$paired) {
        c(paste0(sampleName, "_unclassified_1.fastq"),
          paste0(sampleName, "_unclassified_2.fastq"))
      } else {
        paste0(sampleName, "_unclassified.fastq")
      }
      ezSystem(paste("pigz -f --best -p", ezThreads(),
                     paste(unclassifiedFiles, collapse = " ")))
    }
  }

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppKraken <-
  setRefClass(
    "EzAppKraken",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodKraken
        name <<- "EzAppKraken"
        appDefaults <<- rbind(
          krakenDBOpt = ezFrame(
            Type = "character",
            DefaultValue = "bacteria",
            Description = "kraken database options: viruses bacteria. Default is bacteria"
          ),
          krakenConfidenceOpt = ezFrame(
            Type = "numeric",
            DefaultValue = "0.0",
            Description = "Confidence score threshold (default: 0.0); must be in [0, 1]."
          ),
          krakenPhredOpt = ezFrame(
            Type = "integer",
            DefaultValue = "0",
            Description = "minimum Phred quality, default 0;"
          )
        )
      }
    )
  )
