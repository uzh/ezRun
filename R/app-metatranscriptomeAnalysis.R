###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMetatranscriptomeAnalysis = function(input=NA, output=NA, param=NA, 
                                          htmlFile="00index.html"){
  
  library(purrr)
  library(rtracklayer)
  library(ggplot2)
  library(RColorBrewer)
  library(GO.db)
  dataset = input$meta
  sampleNames = input$getNames() 
  numberOfTopNCategories = param$numberOfTopNCategories
  
  
  dataset = input$meta
  colnames(dataset) <-  gsub(" \\[File\\]","",colnames(dataset))
  annotationFiles <- input$getFullPaths("annotationFile")
  plotLabels <- input$getNames()

  ## Merge annotation files 
  listOfAnnotatedAbundTableOrg <- lapply(annotationFiles,convertDiamondAnnotationToAbund,
                                         feature="organism")
  listOfAnnotatedAbundTableFunc <- lapply(annotationFiles,convertDiamondAnnotationToAbund,
                                         feature="function")
  orgDFforHeatmap <- listOfAbundMerge(listOfAnnotatedAbundTableOrg,plotLabels)
  orgDFforHeatmap <- listOfAbundMerge(listOfAnnotatedAbundTableFunc,plotLabels)
  ##
  setwdNew(basename(output$getColumn("Report")))
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "metagenomeAnnotation.Rmd", 
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="metagenomeAnnotation.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
}
##' @template app-template
##' @templateVar method ezMethodMetatranscriptomeAnalysis()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMetatranscriptomeAnalysis<-
  setRefClass("EzAppMetatranscriptomeAnalysis",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMetatranscriptomeAnalysis
                  name <<- "EzAppMetatranscriptomeAnalysis"
                  appDefaults <<- rbind(numberOfTopNCategories = ezFrame(Type="integer",  DefaultValue="20",
                                                                         Description="How many top N GO and 
                                                                         prot families.")
                  )
                }
              )
  )







