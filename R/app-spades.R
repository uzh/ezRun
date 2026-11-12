###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodSpades = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  opt = param$cmdOptions
  sampleName = input$getNames()
  ##stopifnot((param$paired))
  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  basicOpts <- param$spadesBasicOpt
  if (param$paired) {
    read1 = trimmedInput$getColumn("Read1")
    read2 = trimmedInput$getColumn("Read2")
    readOpt = paste("-1", read1, "-2", read2)
    cmd = paste(
      "spades.py",
      readOpt,
      basicOpts,
      param$spadesPipeOpt,
      "-m",
      param$ram,
      "-o",
      "spades",
      '-t',
      ezThreads(),
      opt,
      "1> ",
      paste0(sampleName, "_spades.log")
    )
    ezSystem(cmd)
  } else {
    read1 = trimmedInput$getColumn("Read1")
    readOpt = paste("-s", read1)
    cmd = paste(
      "spades.py",
      readOpt,
      basicOpts,
      param$spadesPipeOpt,
      "-m",
      param$ram,
      "-o",
      "spades",
      '-t',
      ezThreads(),
      opt,
      "1> ",
      paste0(sampleName, "_spades.log")
    )
    ezSystem(cmd)
  }
  wddir <- "."
  if (grepl("--rna", opt)) {
    sfile <- file.path(wddir, "spades/transcripts.fasta")
    cfile <- file.path(wddir, "spades/hard_filtered_transcripts.fasta")
    if (file.exists(sfile)) {
      ezSystem(paste(
        "cp",
        "spades/transcripts.fasta",
        basename(output$getColumn("Draft"))
      ))
    } else if (file.exists(cfile)) {
      ezSystem(paste(
        "cp",
        "spades/hard_filtered_transcripts.fasta",
        basename(output$getColumn("Draft"))
      ))
    }
    ezSystem(paste("mv", "spades", sampleName))
  } else {
    sfile <- file.path(wddir, "spades/scaffolds.fasta")
    cfile <- file.path(wddir, "spades/contigs.fasta")
    if (file.exists(sfile)) {
      ezSystem(paste(
        "cp",
        "spades/scaffolds.fasta",
        basename(output$getColumn("Draft"))
      ))
    } else if (file.exists(cfile)) {
      ezSystem(paste(
        "cp",
        "spades/contigs.fasta",
        basename(output$getColumn("Draft"))
      ))
    }
    ezSystem(paste("mv", "spades", sampleName))
  }
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppSpades <-
  setRefClass(
    "EzAppSpades",
    contains = "EzApp",
    methods = list(
      ## fastp: read trimming, always runs (ezMethodFastpTrim call in ezMethodSpades).
      ## SPAdes: always runs; Prjibelski 2020 is the paper the SPAdes README asks to cite.
      ## Bankevich 2012: the README asks for it with single-cell data (--sc).
      ## metaSPAdes / plasmidSPAdes: when --meta / --plasmid is in the SPAdes options
      ## (spadesBasicOpt, spadesPipeOpt or cmdOptions, all pasted into the spades.py call).
      ## rnaSPAdes: when --rna is in cmdOptions, the condition ezMethodSpades uses to deliver
      ## the transcripts.
      citation = function(param = list()) {
        spadesOpts <- paste(param$spadesBasicOpt, param$spadesPipeOpt, param$cmdOptions)
        hasOpt <- function(x, flag) {
          isTRUE(grepl(paste0("(^|\\s)", flag, "(\\s|$)"), x))
        }
        c(
          "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560",
          "Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A. & Korobeynikov, A. Using SPAdes de novo assembler. Current Protocols in Bioinformatics 70(1), e102 (2020). https://doi.org/10.1002/cpbi.102",
          if (hasOpt(spadesOpts, "--sc")) "Bankevich, A. et al. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of Computational Biology 19(5), 455-477 (2012). https://doi.org/10.1089/cmb.2012.0021",
          if (hasOpt(spadesOpts, "--meta")) "Nurk, S., Meleshko, D., Korobeynikov, A. & Pevzner, P.A. metaSPAdes: a new versatile metagenomic assembler. Genome Research 27(5), 824-834 (2017). https://doi.org/10.1101/gr.213959.116",
          if (hasOpt(spadesOpts, "--plasmid")) "Antipov, D., Hartwick, N., Shen, M., Raiko, M., Lapidus, A. & Pevzner, P.A. plasmidSPAdes: assembling plasmids from whole genome sequencing data. Bioinformatics 32(22), 3380-3387 (2016). https://doi.org/10.1093/bioinformatics/btw493",
          if (hasOpt(param$cmdOptions, "--rna")) "Bushmanova, E., Antipov, D., Lapidus, A. & Prjibelski, A.D. rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. GigaScience 8(9), giz100 (2019). https://doi.org/10.1093/gigascience/giz100"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSpades
        name <<- "EzAppSpades"
        appDefaults <<- rbind(
          spadesBasicOpt = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "spades basic options: --sc --meta --rna --plasmid. Default is empty for genome assembly without MDA"
          ),
          spadesPipeOpt = ezFrame(
            Type = "character",
            DefaultValue = "--careful",
            Description = "spades pipeline options: --only-assembler --careful"
          )
        )
      }
    )
  )
