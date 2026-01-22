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
