###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodMageckCount <- function(input, output, param) {
  require(Herper)
  local_CondaEnv("gi_mageck", pathToMiniConda = "/usr/local/ngseq/miniforge3")
  sampleName <- input$getNames()
  inputFile <- input$getFullPaths("Read1")

  ## Resolve the library dict / control-sgRNA files, materialising them from the
  ## basic library csv on first use (idempotent + race-safe, see below).
  param <- getMageckReference(param)
  if (length(param[['dictFile']]) != 1L) {
    stop(
      "expected exactly one MAGeCK dict file in library '",
      param[['libName']],
      "'; found ",
      length(param[['dictFile']])
    )
  }

  ## mageck count -- one of two flag variants (with/without control-sgRNA file) is
  ## always run, selected internally, not by a user param. Run through ezSystem so
  ## a non-zero exit stops the job (raw system2 does not) and the cmd is logged.
  hasCtrl <- length(param[['ctrlFile']]) == 1L && nzchar(param[['ctrlFile']])
  cmd <- paste(
    "mageck count",
    "-l",
    shQuote(param[['dictFile']]),
    if (hasCtrl) paste("--control-sgrna", shQuote(param[['ctrlFile']])) else "",
    "--fastq",
    paste(shQuote(inputFile), collapse = " "),
    "-n",
    shQuote(sampleName),
    if (ezIsSpecified(param[['cmdOptions']])) param[['cmdOptions']] else ""
  )
  ezSystem(cmd)

  ## Verify the expected outputs were actually produced before declaring success.
  countFile <- paste0(sampleName, ".count.txt")
  summaryFile <- paste0(sampleName, ".countsummary.txt")
  if (!file.exists(countFile)) {
    stop("mageck count did not produce the expected count file: ", countFile)
  }

  ## Per-sample count QC report from the countsummary mageck already writes
  ## (Gini index, zero-count sgRNAs, %-mapped, total reads). Wrapped so a report
  ## failure never discards an otherwise-good count.
  if (file.exists(summaryFile)) {
    param[['sampleName']] <- sampleName
    ## Per-sample report name: SAMPLE-mode jobs all copy into one shared result
    ## dir, so a fixed name (e.g. 00index.html) would collide across samples.
    ## Registered as a [File] column in the .rb so the framework rsyncs it back
    ## ([Link] columns are not copied).
    tryCatch(
      makeQuartoReport(
        param = param,
        htmlFile = paste0(sampleName, ".html"),
        qmdFile = "MageckCountQC.qmd",
        reportTitle = paste0("MAGeCK Count QC - ", sampleName)
      ),
      error = function(e) {
        ezLog(paste("MageckCountQC report failed:", conditionMessage(e)))
      }
    )
  } else {
    ezLog(paste("no countsummary file found, skipping QC report:", summaryFile))
  }

  return("Success")
}

##' @title Locate the MAGeCK library dict / control-sgRNA files
##' @description Returns \code{param} with \code{dictFile} and \code{ctrlFile}
##'   pointing at the per-library \code{*_MAGeCK.csv} / \code{*_MAGeCK_Ctrl.csv}
##'   files, materialising them from the basic library csv on first use.
getMageckReference <- function(param) {
  dictFile <- list.files(
    param[['libName']],
    pattern = '_MAGeCK\\.csv$',
    full.names = TRUE
  )
  if (length(dictFile) == 0L) {
    ## One-time bootstrap; idempotent and safe to call concurrently.
    prepareMageckLibrary(param[['libName']])
    dictFile <- list.files(
      param[['libName']],
      pattern = '_MAGeCK\\.csv$',
      full.names = TRUE
    )
  }
  param[['dictFile']] <- dictFile
  param[['ctrlFile']] <- list.files(
    param[['libName']],
    pattern = '_MAGeCK_Ctrl\\.csv$',
    full.names = TRUE
  )
  return(param)
}

##' @title Build the MAGeCK library files from the basic library csv
##' @description Converts the basic 4-column library csv
##'   (TranscriptName, Sequence, GeneSymbol, isControl) into the MAGeCK dict file
##'   (\code{*_MAGeCK.csv}) and, when control sgRNAs are flagged, the control file
##'   (\code{*_MAGeCK_Ctrl.csv}). Idempotent: does nothing if the dict already
##'   exists. Race-safe: writes are staged to a pid-suffixed temp file and moved
##'   into place with an atomic rename, so concurrent count jobs cannot observe a
##'   half-written file (no lock file, no stale-lock hazard).
prepareMageckLibrary <- function(libName) {
  dictFile <- list.files(libName, pattern = '_MAGeCK\\.csv$', full.names = TRUE)
  if (length(dictFile) >= 1L) {
    return(invisible(NULL)) # already prepared
  }

  ## The "basic" csv is any csv that is not one of the files we generate.
  allCsv <- list.files(libName, pattern = '\\.csv$', full.names = TRUE)
  generated <- list.files(
    libName,
    pattern = '_MAGeCK(_Ctrl)?\\.csv$',
    full.names = TRUE
  )
  basicFile <- setdiff(allCsv, generated)
  if (length(basicFile) != 1L) {
    stop(
      "no or multiple basic reference file(s) available in library: ",
      libName
    )
  }

  myRef <- ezRead.table(basicFile, row.names = NULL, sep = ',', header = FALSE)
  if (ncol(myRef) < 4L) {
    stop(
      "basic library csv must have 4 columns ",
      "(TranscriptName, Sequence, GeneSymbol, isControl); got ",
      ncol(myRef),
      ": ",
      basicFile
    )
  }
  myRef <- myRef[, 1:4]
  colnames(myRef) <- c('TranscriptName', 'Sequence', 'GeneSymbol', 'isControl')
  myRef[['ID']] <- paste(
    myRef[['TranscriptName']],
    myRef[['Sequence']],
    sep = '_'
  )

  ## Normalise the control flag: accept TRUE/T/1/yes/y (any case) as control.
  isControl <- tolower(trimws(as.character(myRef[['isControl']]))) %in%
    c('true', 't', '1', 'yes', 'y')

  refFile <- sub('\\.csv$', '_MAGeCK.csv', basicFile)
  .mageckAtomicWrite(
    myRef[, c('ID', 'Sequence', 'GeneSymbol')],
    refFile,
    col.names = FALSE,
    row.names = FALSE,
    sep = ','
  )

  if (any(isControl)) {
    ctrlFile <- sub('\\.csv$', '_MAGeCK_Ctrl.csv', basicFile)
    .mageckAtomicWrite(
      data.frame(ID = myRef[['ID']][isControl]),
      ctrlFile,
      col.names = FALSE,
      row.names = FALSE
    )
  }
  invisible(NULL)
}

## Write a table via a pid-suffixed temp file, then atomically rename into place.
.mageckAtomicWrite <- function(x, file, ...) {
  tmp <- paste0(file, ".tmp.", Sys.getpid())
  ezWrite.table(x, tmp, ...)
  file.rename(tmp, file)
}

##' @template app-template
##' @templateVar method ezMethodMageckCount(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppMageckCount <-
  setRefClass(
    "EzAppMageckCount",
    contains = "EzApp",
    methods = list(
      ## mageck count unconditional -- one of two flag variants (with/without
      ## control-sgRNA file) is always run, selected internally, not by a user param.
      citation = function() {
        c(
          "Li, W. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology 15, 554 (2014). https://doi.org/10.1186/s13059-014-0554-4"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodMageckCount
        name <<- "EzAppMageckCount"
        appDefaults <<- rbind(
          libName = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "sgRNA Library Name"
          ),
          cmdOptions = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "additional command line options passed to 'mageck count' (e.g. --sgrna-len, --count-n, --list-seq)"
          )
        )
      }
    )
  )
