###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodMetaPhlAn = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  sampleName <- input$getNames()
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)

  indexName <- as.character(param$metaphlanIndex)
  if (!nzchar(indexName)) {
    stop("metaphlanIndex is empty: no database selected.")
  }
  dbDir <- "/srv/GT/databases/metaphlan_databases"
  if (!file.exists(file.path(dbDir, paste0(indexName, ".pkl")))) {
    stop(sprintf("MetaPhlAn .pkl not found for index '%s' under %s.", indexName, dbDir))
  }

  outProfile <- paste0(sampleName, "_metaphlan.txt")
  outMapout <- paste0(sampleName, ".bowtie2.bz2")
  outLog <- paste0(sampleName, ".metaphlan.log")

  if (param$paired) {
    read1 <- trimmedInput$getColumn("Read1")
    read2 <- trimmedInput$getColumn("Read2")
    readArg <- paste0(read1, ",", read2)
  } else {
    readArg <- trimmedInput$getColumn("Read1")
  }

  cmd <- paste(
    "metaphlan",
    readArg,
    "--db_dir", dbDir,
    "--index", indexName,
    "--input_type fastq",
    "--nproc", ezThreads(),
    "--tmp_dir .",
    "--mapout", outMapout,
    "-o", outProfile,
    param$cmdOptions,
    "1>", outLog, "2>&1"
  )
  ezSystem(cmd)

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMetaPhlAn()
##' @templateVar htmlArg )
##' @description Use this reference class to run MetaPhlAn taxonomic profiling
EzAppMetaPhlAn <-
  setRefClass(
    "EzAppMetaPhlAn",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodMetaPhlAn
        name <<- "EzAppMetaPhlAn"
        appDefaults <<- rbind(
          metaphlanIndex = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "MetaPhlAn bowtie2 index basename under /srv/GT/databases/metaphlan_databases/ (e.g. mpa_vJan25_CHOCOPhlAnSGB_202503)."
          )
        )
      }
    )
  )
