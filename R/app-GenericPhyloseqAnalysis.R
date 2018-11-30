###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodGenericPhyloSeqAnalysis = function(input=NA, output=NA, param=NA, 
                                          htmlFile="00index.html"){
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  require(DESeq2)
  dataset = input$meta
  
  ### Analyzes results with phyloseq: preparing objects to be processed in the Rmd file
  
  ### OTUs
  copyOTUTable <- ezSystem(paste("cp", input$getFullPaths("OTUsCountTable"),"./", sep = " "))
  otuObject <- phyloSeqOTU(basename(input$getFullPaths("OTUsCountTable")))
  
  ### taxonomy
  copytaxaTable <- ezSystem(paste("cp", input$getFullPaths("OTUsToTaxonomy"),"./", sep = " "))
  taxaObject <- phyloSeqTaxa(basename(input$getFullPaths("OTUsToTaxonomy")))

  ### pruning level
  pruneLevel <- param$representativeOTUs
  
  ### Samples
  designMatrix <- dataset[,c("Name","Group [Factor]")]
  sampleObject<- phyloSeqSample(designMatrix)

  ### create, add trees, preprocess and prune phyloseq object 
  physeqObjectNoTree = phyloseq(otuObject, taxaObject,sampleObject)
  treeObject = rtree(ntaxa(physeqObjectNoTree), rooted=TRUE, tip.label=taxa_names(physeqObjectNoTree))
  physeqFullObject <- merge_phyloseq(physeqObjectNoTree,treeObject)
  physeqFullObject <- phyloSeqPreprocess(physeqFullObject)
  physeqFullObject <- prune_taxa(taxa_names(physeqFullObject)[1:pruneLevel], physeqFullObject)

  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "MothurPhyloseqGenericAnalysis.Rmd", 
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="MothurPhyloseqGenericAnalysis.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
}
##' @template app-template
##' @templateVar method ezMethodGenericPhyloSeqAnalysis()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppGenericPhyloSeqAnalysis <-
  setRefClass("EzAppGenericPhyloSeqAnalysis",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodGenericPhyloSeqAnalysis
                  name <<- "EzAppGenericPhyloSeqAnalysis"
                  appDefaults <<- rbind(representativeOTUs = ezFrame(Type="numeric",  DefaultValue="",Description="Number of core OTUs for the samples.")
                  )
                }
              )
  )
