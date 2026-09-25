###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodNanoPlot = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  opt = param$cmdOptions
  sampleName = input$getNames()
  cmd = paste(
    "NanoPlot",
    "-t",
    param$cores,
    "-p",
    paste0(sampleName, "."),
    "--title",
    sampleName,
    "-o",
    sampleName,
    opt,
    "--fastq",
    input$getFullPaths("Read1")
  )
  ezSystem(cmd)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodNanoPlot()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppNanoPlot <-
  setRefClass(
    "EzAppNanoPlot",
    contains = "EzApp",
    methods = list(
      ## NanoPlot unconditional (NanoPack2 is the current paper, NanoPack the original).
      citation = function() {
        c(
          "De Coster, W. & Rademakers, R. NanoPack2: population-scale evaluation of long-read sequencing data. Bioinformatics 39, btad311 (2023). https://doi.org/10.1093/bioinformatics/btad311",
          "De Coster, W. et al. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics 34, 2666-2669 (2018). https://doi.org/10.1093/bioinformatics/bty149"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodNanoPlot
        name <<- "EzAppNanoPlot"
      }
    )
  )
