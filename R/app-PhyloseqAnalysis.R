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
  
  dataset = input$meta
  isGroupThere = param$group
### Analyzes results with phyloseq: preparing objects to be processed in the Rmd file

### OTUs
  ### create phyloseq OTU object
  otuObject <- phyloSeqOTUFromFile(input$getFullPaths("OTUsCountTable"))
  ### create phyloseq Taxa object
  taxaObject <- phyloSeqTaxaFromFile(input$getFullPaths("OTUsToTaxonomyFile"))
  
  ### Eventual groups
  if (isGroupThere){
  designMatrix <- ezRead.table(input$getFullPaths("sampleDescriptionFile"))
  sampleObject <- sample_data(designMatrix)
  physeqObjectNoTree = phyloseq(otuObject, taxaObject, sampleObject)
  }else{
    physeqObjectNoTree = phyloseq(otuObject, taxaObject)
  }
  ##prune OTUS
  pruneLevel <- param$representativeOTUs

### create, add trees, preprocess and prune phyloseq object 

  treeObject = rtree(ntaxa(physeqObjectNoTree), rooted=TRUE, tip.label=taxa_names(physeqObjectNoTree))
  physeqFullObject <- merge_phyloseq(physeqObjectNoTree,treeObject)
  physeqFullObject <- phyloSeqPreprocess(physeqFullObject)
  myTaxa = names(sort(taxa_sums(physeqFullObject), decreasing = TRUE)[1:pruneLevel])
  physeqFullObject <- prune_taxa(myTaxa,physeqFullObject)
  
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
