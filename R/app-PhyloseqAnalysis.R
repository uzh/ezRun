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
  library(gridExtra)
  library(iNEXT)
  
  
  dataset = input$meta
  rank = param$taxonomicRank
  rawCount = param$rawCount
  sampleFraction = as.numeric(param$sampleFraction)
  numTopRanks=param$numTopRanks
  isGroupThere = param$grouping != ""
  if (isGroupThere) {
    group <- param$grouping
    if (param$sampleGroup == "" | param$refGroup == ""){
      stop("Both sample and reference groups must be specified")
    } else{
    sampleGroup <- param$sampleGroup
    refGroup <- param$refGroup
    }
  }
### Analyzes results with phyloseq: preparing objects to be processed in the Rmd file

### load  objects
  ## QC and chimera data
  QCChimeraObjectRData <- input$getFullPaths("RObjectQCChimera")
  QCChimeraObject <- readRDS(QCChimeraObjectRData)
  ## phyloseq object
  physeqObjectRData <- input$getFullPaths("RObjectPhyloseq")
  physeqObjectNoTreeUnfilt <- readRDS(physeqObjectRData)
  

  
### check how many cols sample data has; if only one, add a second otherwise the plot funciton is buggy
  if (isGroupThere){
  N <- ncol(sample_data(physeqObjectNoTreeUnfilt))
  if (N==1) {
   areThereMultVar <- FALSE
    sample_data(physeqObjectNoTreeUnfilt)[,"dummy"] <- sample_data(physeqObjectNoTreeUnfilt)[,1]  
  }  else {
    areThereMultVar <- TRUE
  } 
  ### also ensure the groups are factor 
  sample_data(physeqObjectNoTreeUnfilt)@.Data  <- lapply(sample_data(physeqObjectNoTreeUnfilt)@.Data, as.factor)
  } else {
    areThereMultVar <- FALSE
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
                  appDefaults <<- rbind(taxonomicRank = ezFrame(Type="character",  DefaultValue="Phylum",Description="Rank to be highlighted in the plots"),
                                        numTopRanks = ezFrame(Type="numeric",  DefaultValue="15",Description="Number of top OTUs for the plots"),
                                        rawCount = ezFrame(Type="numeric",  DefaultValue="15",Description="OTUs with fewer than these counts in less than the fraction of samples below will be removed."),
                                        sampleFraction = ezFrame(Type="numeric",  DefaultValue="15",Description="Minimum fraction of samples for which the raw count threshold above applies")
                  )
                }
              )
  )
