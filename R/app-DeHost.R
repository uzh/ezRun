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
##' step used inside EzAppVirDetect). The number of unmapped (i.e. kept) read
##' pairs is reported in a small stats file alongside the fastq output.
ezMethodDeHost <- function(input = NA, output = NA, param = NA) {
  param$fastpCompression <- 9
  ref <- getBowtie2Reference(param)
  sampleName <- input$getNames()
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)

  defOpt <- paste("-p", param$cores)
  readGroupOpt <- paste0(
    "--rg-id ",
    sampleName,
    " --rg SM:",
    sampleName,
    " --rg LB:RGLB_",
    sampleName,
    " --rg PL:illumina",
    " --rg PU:RGPU_",
    sampleName
  )
  ## keep only the reads/pairs for which nothing mapped to the host genome;
  ## filter directly on bowtie2's output so the (typically much larger) fully
  ## aligned BAM is never written to disk
  if (param$paired) {
    flagFilter <- "-f 12 -F 256"
  } else {
    flagFilter <- "-f 4 -F 256"
  }
  cmd <- paste(
    "bowtie2",
    param$cmdOptions,
    defOpt,
    readGroupOpt,
    "-x",
    ref,
    if (param$paired) "-1",
    trimmedInput$getColumn("Read1"),
    if (param$paired) paste("-2", trimmedInput$getColumn("Read2")),
    "2>",
    paste0(sampleName, "_bowtie2.log"),
    "|",
    "samtools view -b",
    flagFilter,
    "- > host_unmapped.bam"
  )
  ezSystem(cmd)
  file.remove(trimmedInput$getColumn("Read1"))
  if (param$paired) {
    file.remove(trimmedInput$getColumn("Read2"))
  }

  ## count the kept reads/pairs
  nRecords <- as.integer(ezSystem(
    "samtools view -c host_unmapped.bam",
    intern = TRUE,
    stopOnFailure = TRUE
  ))
  nUnmapped <- if (param$paired) nRecords %/% 2L else nRecords
  statsFile <- paste0(sampleName, "_dehost_stats.tsv")
  ezWrite.table(
    data.frame(
      Sample = sampleName,
      UnmappedReadPairs = nUnmapped,
      check.names = FALSE
    ),
    file = statsFile,
    row.names = FALSE
  )
  ezSystem(paste0(
    "echo 'unmapped ",
    if (param$paired) "read pairs" else "reads",
    " kept: ",
    nUnmapped,
    "' >> ",
    paste0(sampleName, "_bowtie2.log")
  ))

  ## group mates back together (bowtie2 already emits them adjacently, but
  ## collate is cheap insurance) and write out gzipped fastq directly
  r1Fastq <- basename(output$getColumn("Read1"))
  if (param$paired) {
    r2Fastq <- basename(output$getColumn("Read2"))
    cmd <- paste(
      "samtools collate -Ou host_unmapped.bam",
      "|",
      "samtools fastq -@",
      param$cores,
      "-1",
      r1Fastq,
      "-2",
      r2Fastq,
      "-0 /dev/null -s /dev/null -"
    )
  } else {
    cmd <- paste("samtools fastq -@", param$cores, "-0", r1Fastq, "host_unmapped.bam")
  }
  ezSystem(cmd)
  file.remove("host_unmapped.bam")

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
