###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodQuast = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  opt = param$cmdOptions
  sampleName = input$getNames()
  draft = input$getFullPaths("Draft")
  if (ezIsSpecified(param$refGenome)) {
    ref = param$refGenome
    if (ezIsSpecified(param$refGene)) {
      gene = param$refGene
      cmd = paste(
        "quast.py",
        "-R",
        ref,
        "-G",
        gene,
        "-o",
        sampleName,
        '-t',
        ezThreads(),
        opt,
        draft,
        "1> ",
        paste0(sampleName, "_quast.log")
      )
    } else {
      cmd = paste(
        "quast.py",
        "-R",
        ref,
        "-o",
        sampleName,
        '-t',
        ezThreads(),
        opt,
        draft,
        "1> ",
        paste0(sampleName, "_quast.log")
      )
    }
  } else {
    cmd = paste(
      "quast.py",
      "-o",
      sampleName,
      '-t',
      ezThreads(),
      opt,
      draft,
      "1> ",
      paste0(sampleName, "_quast.log")
    )
  }
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppQuast <-
  setRefClass(
    "EzAppQuast",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodQuast
        name <<- "EzAppQuast"
        appDefaults <<- rbind(
          refGenome = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "full path to a reference genome as a multi-fasta file"
          ),
          refGene = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "full path to a gene annotation file of the reference genome. Must be in gff or bed format"
          )
        )
      }
    )
  )
