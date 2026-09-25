###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodFlash = function(input = NA, output = NA, param = NA) {
  opt = param$cmdOptions
  sampleName = input$getNames()
  param$fastpCompression = 9
  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  if (!param$skipFlash) {
    stopifnot((param$paired))
    cmd = paste(
      "flash",
      trimmedInput$getColumn("Read1"),
      trimmedInput$getColumn("Read2"),
      "-o",
      sampleName,
      '-t',
      ezThreads(),
      opt,
      "1>> ",
      paste0(sampleName, "_preprocessing.log")
    )
    ezSystem(cmd)
    cmd = paste0('pigz --best ', sampleName, '.extendedFrags.fastq')
    ezSystem(cmd)
    cmd = paste(
      'mv',
      paste0(sampleName, '.extendedFrags.fastq.gz'),
      paste0(sampleName, '.R1.fastq.gz')
    )
    ezSystem(cmd)
  } else {
    ezSystem(paste(
      'mv',
      paste0(sampleName, '-trimmed_R1.fastq.gz'),
      paste0(sampleName, '.R1.fastq.gz')
    ))
    if (param$paired) {
      ezSystem(paste(
        'mv',
        paste0(sampleName, '-trimmed_R2.fastq.gz'),
        paste0(sampleName, '.R2.fastq.gz')
      ))
    }
  }
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodFlash(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppFlash <-
  setRefClass(
    "EzAppFlash",
    contains = "EzApp",
    methods = list(
      ## FLASH defaults checked against FLASH 1.2.11 --help.
      ## fastp unconditional; FLASH gated on skipFlash = FALSE (paired only).
      citation = function() {
        c(
          "Magoč, T. & Salzberg, S.L. FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics 27, 2957-2963 (2011). https://doi.org/10.1093/bioinformatics/btr507",
          "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560"
        )
      },
      methods_facts = function(param = list()) {
        c(
          ## ezMethodFlash -> ezMethodFastpTrim (app-flash.R:12)
          methodsFastpFacts(param, "read merging"),
          ## app-flash.R:13-33
          if (isFALSE(as.logical(param$skipFlash))) "Overlapping mates were merged with FLASH using its defaults unless cmdOptions set them (minimum overlap 10 bp, maximum overlap 65 bp, maximum mismatch density 0.25, innie orientation only), and only the merged reads were delivered; pairs FLASH could not merge were discarded.",
          ## app-flash.R:34-46
          if (isTRUE(as.logical(param$skipFlash))) "FLASH was skipped and the fastp-trimmed reads were delivered unmerged."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodFlash
        name <<- "EzAppFlash"
        appDefaults <<- rbind(
          skipFlash = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "run or skip flash"
          )
        )
      }
    )
  )
