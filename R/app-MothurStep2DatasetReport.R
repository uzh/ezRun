###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurStep2DatasetReport = function(input=NA, output=NA, param=NA, 
                                           htmlFile="00index.html"){
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  require(DESeq2)
  library(Matrix)
  library(magic)
  library(ape)
  library(limma)
  library(RColorBrewer)
  library(gplots)
  require(knitr)
  require(kableExtra)
  require(SummarizedExperiment)
  require(webshot)
  require(htmlwidgets)
  
  dataset = input$meta
  fileNames <- as.vector(input$getNames())
  ### Further report on Mothur pipeline and analysis of the  results with phyloseq
  
  ### create phyloseq OTU object
  otuObject <- phyloSeqOTUFromFile(input$getFullPaths("OTUsCountTable"))
  ### create phyloseq Taxa object
  taxaObject <- phyloSeqTaxaFromFile(input$getFullPaths("OTUsToTaxonomyFile"))
  
  ### Add sample object (TODO, derive it from step1)
  if (param$group){
    designMatrix <- ezRead.table(input$getFullPaths("sampleDescriptionFile"))
    sampleObject <- sample_data(designMatrix)
    physeqPacBioNoTree = phyloseq(otuObject, taxaObject, sampleObject)
  }else{
    physeqPacBioNoTree = phyloseq(otuObject, taxaObject)
  }
  ##prune OTUS
  pruneLevel <- param$representativeOTUs
  
  ### create, add trees, preprocess and prune phyloseq object
  treeObject = rtree(ntaxa(physeqObjectNoTree), rooted=TRUE, tip.label=taxa_names(physeqObjectNoTree))
  physeqFullObject <- merge_phyloseq(physeqObjectNoTree,treeObject)
  physeqFullObject <- phyloSeqPreprocess(physeqFullObject)
  myTaxa = names(sort(taxa_sums(physeqFullObject), decreasing = TRUE)[1:pruneLevel])
  physeqFullObject <- prune_taxa(myTaxa,physeqFullObject)
  
  ## Copy the style files and templates
  RmarkdownFile <- "MothurStep2DatasetReport.Rmd"
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", RmarkdownFile, 
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input=RmarkdownFile, envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
}
##' @template app-template
##' @templateVar method ezMethodMothurStep2DatasetReport()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurStep2DatasetReport <-
  setRefClass("EzAppMothurStep2DatasetReport",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurStep2DatasetReport
                  name <<- "EzAppMothurStep2DatasetReport"
                  appDefaults <<- rbind(representativeOTUs = ezFrame(Type="numeric",  DefaultValue="",Description="Number of core OTUs for the samples."),
                                        group = ezFrame(Type="logical",  DefaultValue="true",Description="Experiment with groups.")
                  )
                }
              )
  )
