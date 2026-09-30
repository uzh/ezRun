###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Depletes host reads using bowtie2
##' @description Maps (optionally trimmed) reads to a host genome with bowtie2 and
##' keeps only the reads/read pairs for which nothing mapped, writing them out as
##' gzipped fastq files. Used to remove host contamination ahead of downstream
##' metagenomic/viral analyses (e.g. as a standalone version of the host-removal
##' step used inside EzAppVirDetect). Unless \code{param$refBuild} is itself
##' human, human contamination is always depleted first (same policy as
##' EzAppVirDetect) before depleting against the user-selected host build.
##' The number of unmapped (i.e. kept) read pairs from the final stage is
##' reported in a small stats file alongside the fastq output.
ezMethodDeHost <- function(input = NA, output = NA, param = NA) {
  param$fastpCompression <- 9
  sampleName <- input$getNames()
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)

  logFile <- paste0(sampleName, "_bowtie2.log")

  ## always deplete human contamination first, unless the user-selected host
  ## build already is human
  refChain <- c(Host = param$refBuild)
  if (!grepl("^Homo_sapiens", param$refBuild)) {
    refChain <- c(Human = DEFAULT_HUMAN_REFBUILD, refChain)
  }

  curR1 <- trimmedInput$getColumn("Read1")
  curR2 <- if (param$paired) trimmedInput$getColumn("Read2") else NULL
  nKept <- NA_integer_

  for (stageName in names(refChain)) {
    res <- depleteAgainstReference(
      read1 = curR1,
      read2 = curR2,
      param = param,
      refBuild = refChain[[stageName]],
      cmdOptions = param$cmdOptions,
      stageLabel = stageName,
      logFile = logFile
    )
    file.remove(curR1)
    if (param$paired) {
      file.remove(curR2)
    }
    curR1 <- res$read1
    curR2 <- res$read2
    nKept <- res$count
  }

  statsFile <- paste0(sampleName, "_dehost_stats.tsv")
  ezWrite.table(
    data.frame(
      Sample = sampleName,
      UnmappedReadPairs = nKept,
      check.names = FALSE
    ),
    file = statsFile,
    row.names = FALSE
  )
  ezSystem(paste0(
    "echo 'unmapped ",
    if (param$paired) "read pairs" else "reads",
    " kept: ",
    nKept,
    "' >> ",
    logFile
  ))

  ## move the final stage's kept fastq to the declared output file names
  file.rename(curR1, basename(output$getColumn("Read1")))
  if (param$paired) {
    file.rename(curR2, basename(output$getColumn("Read2")))
  }

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodDeHost(input=NA, output=NA, param=NA)
##' @description Use this reference class to deplete reads mapping to a host genome with bowtie2
##' @seealso \code{\link{getBowtie2Reference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppDeHost <-
  setRefClass(
    "EzAppDeHost",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodDeHost
        name <<- "EzAppDeHost"
      }
    )
  )
