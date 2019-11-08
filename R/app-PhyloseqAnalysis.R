###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodPhyloSeqAnalysis = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  require(DESeq2)
  library(limma)
  library(RColorBrewer)
  library(gplots)
  require(kableExtra)
  require(knitr)
  library(pheatmap)
  library(ggpubr)
  library(xfun)
  library(vegan)
  library(pryr)
  library(iNEXT)
  
  dataset = input$meta
  isGroupThere = param$group
  rank = param$taxonomicRank
  rawCount = param$rawCount
  sampleFraction = as.numeric(param$sampleFraction)
  numTopRanks=param$numTopRanks
### Analyzes results with phyloseq: preparing objects to be processed in the Rmd file

### load phyloseq object
  physeqObjectRData <- input$getFullPaths("RObjectPhyloseq")
  physeqObjectNoTreeUnfilt <- readRDS(physeqObjectRData)
  
### check how many cols sample data has; if only one, add a second otherwise the plot funciton is buggy
  N <- ncol(sample_data(physeqObjectNoTreeUnfilt))
  if (N==1) {
   areThereMultVar <- FALSE
    sample_data(physeqObjectNoTreeUnfilt)[,"dummy"] <- sample_data(physeqObjectNoTreeUnfilt)[,1]  
  }  else {
    areThereMultVar <- TRUE
   }

### Filtering step
  physeqFullObject <- phyloSeqPreprocess(physeqObjectNoTreeUnfilt,rawCount,sampleFraction)

### run report  
  setwdNew(basename(output$getColumn("Report")))
  
  if (isGroupThere){
    markdownFile <- "PhyloseqReport.Rmd"
  }else{
    markdownFile <- "PhyloseqReportNoGroup.Rmd"
  }
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", markdownFile, 
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input=markdownFile, envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  return("Success")
}
##' @template app-template
##' @templateVar method ezMethodPhyloSeqAnalysis()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppPhyloSeqAnalysis <-
  setRefClass("EzAppPhyloSeqAnalysis",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodPhyloSeqAnalysis
                  name <<- "EzAppPhyloSeqAnalysis"
                  appDefaults <<- rbind(representativeOTUs = ezFrame(Type="numeric",  DefaultValue="",Description="Number of core OTUs for  samples."),
                                        group = ezFrame(Type="logical",  DefaultValue="true",Description="Experiment with groups.")
                  )
                }
              )
  )
