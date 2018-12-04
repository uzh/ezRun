###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurStep1SampleReport = function(input=NA, output=NA, param=NA, 
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
  
  ### function to generate plost from mothur kabled summary tables
  plotFromMothurSumm <- function(x){
    tableDF <-x$mergedTable
    sampleInfo <- x$aboveHeader
    sampleIDs <- names(x$aboveHeader)
    xAxis <- rownames(tableDF)
    tableDF <- data.frame(tableDF, stringsAsFactors = F)
    sampleNameList <- list()
    for (k in 1:length(sampleInfo)) {
      cc <- sampleInfo[k]
      sampleNameList[[k]] <- rep(names(cc),as.numeric(cc))
    } 
    colWithSampleNames <- as.factor(unlist(sampleNameList))
    DFforPlot <- data.frame(percentile = xAxis, length=tableDF$nbases, 
                            ambigs=tableDF$ambigs,
                            homopol=tableDF$polymer, numSeqs=as.factor(tableDF$numSeqs),
                            sample=colWithSampleNames, stringsAsFactors = F)
    caz <- melt(DFforPlot)
    plot <- ggplot(caz, aes(x=reorder(numSeqs, value), y=value, color=variable))+
      geom_point() + facet_wrap(vars(sample),ncol = 2) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size =8), axis.title.x = element_blank())
    finalPlot <- plot + scale_x_discrete(breaks=quantile(as.numeric(as.character(caz$numSeqs)),
                                                         probs = seq(0, 1, 0.05), type=1))
    return(finalPlot)
  }
  
  
  ### prepare files for Rmd 
  summaryTablesToReport <- grep("Summary",names(listOfListAllFiles), value = T)
  summaryTablesToReport <- grep("stepConvergence", summaryTablesToReport, invert = T, value = T)
  summaryTablesToReport <- listOfListAllFiles[summaryTablesToReport]
  
  finalListOfSummaryTables <- imap(summaryTablesToReport,function(x,y)
    createSummaryTableForKableExtra(x))
  ## All files ready
  
  ## create final output dir
  setwdNew(basename(output$getColumn("Report")))
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "MothurStep1SampleReport.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="MothurStep1SampleReport.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
}

##' @template app-template
##' @templateVar method ezMethodMothurStep1SampleReport()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurStep1SampleReport <-
  setRefClass("EzAppMothurStep1SampleReport",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurStep1SampleReport
                  name <<- "EzAppMothurStep1SampleReport"
                }
              )
  )