###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurDataAnalysisReport = function(input=NA, output=NA, param=NA, 
                                   htmlFile="00index.html"){
  
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  library(scales)
  require(gridExtra)
  dataset = input$meta

  
  ## Define input for rmd file
  RawDataSummary <- input$getFullPaths("RawDataSummary")
  lapply(RawDataSummary, function(x) ezSystem(paste("cp",basename(x),"./")))
  DeduppedSummary <- input$getFullPaths("DeduppedSummary")
  LenAndHomopSummary <- input$getFullPaths("LenAndHomopSummary")
  MapFiltSummary <- input$getFullPaths("MapFiltSummary")
  ChimeraPlot <- input$getFullPaths("ChimeraPlot")
  PreClusteredAndChimeraSummary <- input$getFullPaths("PreClusteredAndChimeraSummary")
  stepConvergence <- input$getFullPaths("stepConvergence")
  nRowsGrid <- param$rowsInPlotGrid
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "MothurDataAnalysisReport.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="MothurDataAnalysisReport.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
}

##' @template app-template
##' @templateVar method ezMethodMothurDataAnalysisReport()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurDataAnalysisReport <-
  setRefClass("EzAppMothurDataAnalysisReport",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurDataAnalysisReport
                  name <<- "EzAppMothurDataAnalysisReport"
                  appDefaults <<- rbind(rowsInPlotGrid = ezFrame(Type="integer",
                                                                DefaultValue="2",
                                                                Description="Rows to arrange plots")
                                        
                  )
                }
              )
  )