###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurPhyloSeqAnalysis = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  require(DESeq2)
  library(limma)
  dataset = input$meta
  isGroupThere = param$Group
### Analyzes results with phyloseq: preparing objects to be processed in the Rmd file

### OTUs
  ### create phyloseq OTU object
  otuObject <- phyloSeqOTUFromFile(input$getFullPaths("OTUsCountTable"))
  ### create phyloseq Taxa object
  taxaObject <- phyloSeqTaxaFromFile(input$getFullPaths("OTUsToTaxonomyFile"))
  
  ### Add sample object (TODO, derive it from step1)
  if (param$group){
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

  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "PhyloseqReport.Rmd", 
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="PhyloseqReport.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
}
##' @template app-template
##' @templateVar method ezMethodMothurPhyloSeqAnalysis()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurPhyloSeqAnalysis <-
  setRefClass("EzAppMothurPhyloSeqAnalysis",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurPhyloSeqAnalysis
                  name <<- "EzAppMothurPhyloSeqAnalysis"
                  appDefaults <<- rbind(RepresentativeOTUs = ezFrame(Type="numeric",  DefaultValue="",Description="Number of core OTUs for  samples.")
                  )
                }
              )
  )
