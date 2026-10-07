###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodHifiasm = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  opt = param$cmdOptions
  sampleName = input$getNames()
  inputFileType = param$inputType
  if (inputFileType == "HiFi") {
    cmd = paste(
      "hifiasm",
      "-o",
      sampleName,
      paste("-t", ezThreads()),
      paste("--n-hap", param$ploidy),
      "--primary",
      opt,
      input$getFullPaths("Read1"),
      "2> ",
      paste0(sampleName, "_hifiasm.log")
    )
    ezSystem(cmd, wait = TRUE)
    cmd = paste(
      "awk",
      "'/^S/{print",
      paste0('"', '>', '"', '$2;print'),
      "$3}'",
      paste0(sampleName, ".p_ctg.gfa"),
      ">",
      paste0(sampleName, ".p_ctg.fa")
    )
    system(cmd, wait = TRUE)
    return("Success")
  } else if (inputFileType == "ONT") {
    cmd = paste(
      "hifiasm",
      "-o",
      sampleName,
      paste("-t", ezThreads()),
      paste("--n-hap", param$ploidy),
      "--primary",
      opt,
      "--ont",
      input$getFullPaths("Read1"),
      "2> ",
      paste0(sampleName, "_hifiasm.log")
    )
    ezSystem(cmd, wait = TRUE)
    cmd = paste(
      "awk",
      "'/^S/{print",
      paste0('"', '>', '"', '$2;print'),
      "$3}'",
      paste0(sampleName, ".p_ctg.gfa"),
      ">",
      paste0(sampleName, ".p_ctg.fa")
    )
    system(cmd, wait = TRUE)
    return("Success")
  }
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppHifiasm <-
  setRefClass(
    "EzAppHifiasm",
    contains = "EzApp",
    methods = list(
      ## hifiasm unconditional (Cheng 2021); ONT mode (inputType = ONT) is the Cheng 2026 paper.
      citation = function(param = list()) {
        c(
          "Cheng, H. et al. Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nature Methods 18, 170-175 (2021). https://doi.org/10.1038/s41592-020-01056-5",
          if (identical(param$inputType, "ONT")) "Cheng, H. et al. Efficient near-telomere-to-telomere assembly of nanopore simplex reads. Nature 655, 166-173 (2026). https://doi.org/10.1038/s41586-026-10105-6"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodHifiasm
        name <<- "EzAppHifiasm"
        appDefaults <<- rbind(
          inputType = ezFrame(
            Type = "character",
            DefaultValue = "HiFi",
            Description = "PacBio HiFi reads or Oxford Nanopore reads"
          ),
          ploidy = ezFrame(
            Type = "integer",
            DefaultValue = "2",
            Description = "number of haplotypes"
          )
        )
      }
    )
  )
