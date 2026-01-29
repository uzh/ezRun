###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodProdigal = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  opt = param$cmdOptions
  sampleName = input$getNames()
  draft = input$getFullPaths("Draft")
  genes = paste0(sampleName, ".genes")
  proteins = paste0(sampleName, ".proteins")
  cmd = paste(
    "prodigal -i",
    draft,
    "-o",
    genes,
    "-a",
    proteins,
    "-f",
    param$format,
    "-g",
    param$translation_table,
    opt,
    "1>",
    paste0(sampleName, "_prodigal.log")
  )
  ezSystem(cmd)
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodProdigal
##' @description Use this reference class to run
##' @seealso \code{\link{getPbmm2Reference}}
EzAppProdigal <-
  setRefClass(
    "EzAppProdigal",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodProdigal
        name <<- "EzAppProdigal"
        appDefaults <<- rbind(
          translation_table = ezFrame(
            Type = "integer",
            DefaultValue = "11",
            Description = "translation table to use (default 11)"
          ),
          format = ezFrame(
            Type = "character",
            DefaultValue = "gbk",
            Description = "Select output format (gbk, gff, or sco).  Default is gbk."
          )
        )
      }
    )
  )
