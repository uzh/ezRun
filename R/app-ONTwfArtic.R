# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodONTwfArtic <- function(input = NA, output = NA, param = NA) {
  
  dataset = input$getFullPaths("Read1")
  samplesheet = input$getFullPaths("SampleSheet")
  scheme_version_param <- param$schemeVersion
  normalise_param <- param$normalise
  
  #Run the new 
  setwdNew(basename(output$getColumn("ResultDir")))
  cmd = paste("nextflow run /srv/GT/software/epi2me-labs/wf-artic/", 
              "--fastq", dataset, 
              "--sample_sheet", samplesheet,
              "--normalise", normalise_param,
              "--scheme_version", scheme_version_param, 
              "--profile singularity")
  ezSystem(cmd)
  return("Success")
}
  
##' @template app-template
##' @templateVar method ezMethodONTwfArtic()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppONTwfArtic <-
  setRefClass("EzAppONTwfArtic",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodONTwfArtic
                  name <<- "EzAppONTwfArtic"
                  appDefaults <<- rbind(normalise = ezFrame(Type="integer",  DefaultValue="200",Description="depth ceiling for depth of coverage normalization"),
                                        schemeVersion = ezFrame(Type="character",  DefaultValue="ARTIC/V3",Description="Primer scheme version"))
                }
              )
  )
  