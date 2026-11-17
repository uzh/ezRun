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
##' @templateVar method ezMethodPsortb
##' @description Use this reference class to run
EzAppPsortb <-
  setRefClass(
    "EzAppPsortb",
    contains = "EzApp",
    methods = list(
      ## PSORTb 3.0: always runs (psortb_app call in ezMethodPsortb, R/app-psortb.R). seqtk only
      ## subsets the protein FASTA and is not cited. This file is collated after app-psortb.R, so
      ## this EzAppPsortb definition is the one in effect; kept identical to the copy there.
      citation = function(param = list()) {
        c(
          "Yu, N.Y. et al. PSORTb 3.0: improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes. Bioinformatics 26(13), 1608-1615 (2010). https://doi.org/10.1093/bioinformatics/btq249"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodPsortb
        name <<- "EzAppPsortb"
        appDefaults <<- rbind(
          org = ezFrame(
            Type = "character",
            DefaultValue = "--negative",
            Description = "type of organism: gram negative/ gram positive bacteria or archea"
          )
        )
      }
    )
  )
