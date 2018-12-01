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
  require(grid)
  require(gridExtra)
  library(gtable)
  library(purrr)
  require(knitr)
  require(kableExtra)
  require(SummarizedExperiment)
  require(webshot)
  require(htmlwidgets)
  library(cowplot)
  
  ## Create list of summary files for tables
  
  dataset = input$meta
  relevantColumns <- gsub(" \\[File\\]","",grep("File",colnames(dataset), value = T))
  colnames(dataset) <-  gsub(" \\[File\\]","",colnames(dataset))
  allColumns <- dataset[,relevantColumns]
  plotLabels <- input$getNames()
  ## Copy all files locally
  copyLoopOverFiles <- function(x){ 
    lapply(x,function(x) ezSystem(paste("cp",file.path(DEMO_DATA_ROOT,x),"./")))
  }
  listOfListAllFiles <- as.list(allColumns)
  lapply(listOfListAllFiles,copyLoopOverFiles)


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
                }
              )
  )