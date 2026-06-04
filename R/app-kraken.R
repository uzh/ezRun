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
  sampleName <- input$getNames()
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)

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
  unclassifiedFasta <- paste0(sampleName, "_unclassified.fasta")

  cmd <- paste(
    "k2 classify",
    "--db", dbArg,
    pairedFlag,
    "--use-names",
    "--minimum-hit-groups", minHitGroups,
    reportMinimizerFlag,
    "--confidence", conOpt,
    "--minimum-base-quality", phredOpt,
    "--output", outTxt,
    "--report", outReport,
    "--unclassified-out", unclassifiedFasta,
    "--threads", ezThreads(),
    param$cmdOptions,
    readOpt,
    "1>", outLog
  )
  ezSystem(cmd)

  # Krona must run on the uncompressed per-read output
  cmd <- paste(
    "ktImportTaxonomy -q 2 -t 3",
    "-n", shQuote(sampleName),
    outTxt,
    "-o", outHtml
  )
  ezSystem(cmd)

  ezSystem(paste("gzip -f", outTxt))
  ezSystem(paste("gzip -f", unclassifiedFasta))

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
