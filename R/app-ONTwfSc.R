# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodONTwfSc <- function(input = NA, output = NA, param = NA) {
  
  sampledataset = input$getFullPaths("Read1")
  samplename = input$getNames()
  refbuild = param$refbuild
  kitname = param$kitname
  kitversion = param$kitversion
  expCells = param$expCells
  
  cmd = paste("nextflow run /srv/GT/software/epi2me-labs/wf-single-cell/", 
              "-w", paste0(samplename,"/workspace"),
              "--fastq", sampledataset,
              "--kit_name", kitname,
              "--kit_version", kitversion,
              "--expected_cells", expCells,
              "--ref_genome_dir", paste0("/srv/GT/reference/ont-wf-single-cell/", refbuild),
              "--out_dir", samplename,
              "-profile singularity",
              "--plot_umaps"
              )
  ezSystem(cmd)
  return("Success")
}
  
##' @template app-template
##' @templateVar method ezMethodONTwfSc()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppONTwfSc <-
  setRefClass("EzAppONTwfSc",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodONTwfSc
                  name <<- "EzAppONTwfSc"
                  appDefaults <<- rbind(expected_cells = ezFrame(Type="integer",  DefaultValue="500",Description="number of expected cells"),
                                        kitname = ezFrame(Type="character",  DefaultValue="3prime",Description="10x kit name"))
                }
              )
  )
  