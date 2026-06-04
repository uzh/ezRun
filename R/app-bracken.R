###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodBracken = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  sampleName <- input$getNames()
  krakenReport <- input$getFullPaths("KrakenReport")

  # brackenDBOpt is "<DBname>/<N>mers" (one entry per kmer_distrib file)
  dbOpt <- param$brackenDBOpt
  m <- regmatches(dbOpt, regexec("^(.+)/(\\d+)mers$", dbOpt))[[1]]
  if (length(m) != 3L) {
    stop(sprintf("Unexpected brackenDBOpt value %s; expected '<DBname>/<N>mers'.", dbOpt))
  }
  dbName <- m[2]
  readLen <- m[3]
  level <- if (!is.null(param$brackenLevel)) as.character(param$brackenLevel) else "S"
  threshold <- if (!is.null(param$brackenThreshold)) as.character(param$brackenThreshold) else "0"

  dbDir <- file.path("/srv/GT/databases/kraken2", dbName)

  outAbundance <- paste0(sampleName, ".bracken")
  outReport <- paste0(sampleName, ".bracken.report.txt")
  outLog <- paste0(sampleName, ".bracken.log")

  cmd <- paste(
    "bracken",
    "-d", dbDir,
    "-i", krakenReport,
    "-o", outAbundance,
    "-w", outReport,
    "-r", readLen,
    "-l", level,
    "-t", threshold,
    param$cmdOptions,
    "1>", outLog, "2>&1"
  )
  ezSystem(cmd)

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodBracken()
##' @templateVar htmlArg )
##' @description Use this reference class to run Bracken on a Kraken2 report
EzAppBracken <-
  setRefClass(
    "EzAppBracken",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBracken
        name <<- "EzAppBracken"
        appDefaults <<- rbind(
          brackenDBOpt = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Bracken DB + read length as '<DBname>/<N>mers'; resolves to /srv/GT/databases/kraken2/<DBname>/databaseNmers.kmer_distrib."
          ),
          brackenLevel = ezFrame(
            Type = "character",
            DefaultValue = "S",
            Description = "Taxonomic level for re-estimation: D, P, C, O, F, G, S, S1."
          ),
          brackenThreshold = ezFrame(
            Type = "integer",
            DefaultValue = "0",
            Description = "Minimum reads at the chosen level before re-estimation."
          )
        )
      }
    )
  )
