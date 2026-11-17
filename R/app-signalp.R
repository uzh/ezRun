###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodSignalP = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  opt = param$cmdOptions
  org = param$org
  oformat = param$oformat
  pformat = param$pformat
  sampleName = input$getNames()
  proteins = input$getFullPaths("Proteins")
  cmd = paste(
    "signalp -fasta",
    proteins,
    "-format",
    oformat,
    "-gff3 -org",
    org,
    "-plot",
    pformat,
    "-prefix",
    sampleName,
    opt,
    "1>",
    paste0(sampleName, "_signalp.log")
  )
  ezSystem(cmd)
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodSignalP
##' @description Use this reference class to run
EzAppSignalP <-
  setRefClass(
    "EzAppSignalP",
    contains = "EzApp",
    methods = list(
      citation = function(param = list()) {
        c(
          "Almagro Armenteros, J.J. et al. SignalP 5.0 improves signal peptide predictions using deep neural networks. Nature Biotechnology 37(4), 420-423 (2019). https://doi.org/10.1038/s41587-019-0036-z"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSignalP
        name <<- "EzAppSignalP"
        appDefaults <<- rbind(
          org = ezFrame(
            Type = "character",
            DefaultValue = "gram-",
            Description = "type of organism: gram- or gram+ bacteria, arch (Archaea) or euk (Eukarya)"
          ),
          oformat = ezFrame(
            Type = "character",
            DefaultValue = "short",
            Description = "output format: long adds per-protein plots, short gives the predictions only"
          ),
          pformat = ezFrame(
            Type = "character",
            DefaultValue = "png",
            Description = "plot format for long output: png, eps or none"
          )
        )
      }
    )
  )
