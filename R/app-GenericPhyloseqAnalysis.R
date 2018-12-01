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
  isGroupThere <- param$Group
  ### Analyzes results with phyloseq: preparing objects to be processed in the Rmd file
  
  relevantColumns <- c("OTUsCountTable","OTUsToTaxonomyFile")
  colnames(dataset) <-  gsub(" \\[File\\]","",colnames(dataset))
  colnames(dataset) <-  gsub(" \\[Factor\\]","",colnames(dataset))
  allColumns <- dataset[,relevantColumns]
  ## Copy all files locally
  copyLoopOverFiles <- function(x){ 
    lapply(x,function(x) ezSystem(paste("cp",file.path(DEMO_DATA_ROOT,x),"./")))
  }
  listOfListAllFiles <- as.list(allColumns)
  lapply(listOfListAllFiles,copyLoopOverFiles)
  
### create list of phyloseq OTU object

  k=0
  OTUsCount <- list()
  OTUsCountNoLabel <- list()
  otuObject <-list()
  Group <- vector()
    for (file in listOfListAllFiles["OTUsCountTable"]$OTUsCountTable){
    k=k+1
    rawFile <- basename(file)
    OTUsCount[[k]] <- read.delim(rawFile, header = T,stringsAsFactors = FALSE,check.names = FALSE)
    relCols <- grep("Group",colnames(OTUsCount[[k]]),invert = T)
    OTUsCountNoLabel[[k]] <- as.matrix(OTUsCount[[k]][,relCols])
    Group[k] <- fileNames[k]
    otuObject[[k]] <- phyloSeqOTU(data.frame(Group = Group[k],OTUsCountNoLabel[[k]],
                                             check.names = F))
}

  ### create list of phyloseq taxa object  
  k=0
  OTUsToTaxonomyDF <- list()
  taxaObject <- list()
  for (file in listOfListAllFiles["OTUsToTaxonomyFile"]$OTUsToTaxonomyFile){
    k=k+1
    rawFile <- basename(file)
    OTUsToTaxonomyDF[[k]] <- read.delim(rawFile, header = T,stringsAsFactors = FALSE)
    taxaObject[[k]] <- phyloSeqTaxa(OTUsToTaxonomyDF[[k]])
    }
  
  ###
  phyObjList <- mapply(phyloseq,otuObject, taxaObject)
  fullPhySeqObj <- Reduce(function(x, y) merge_phyloseq(x, y), phyObjList)
 
   ### Add sample object eventually with Group
 if (isGroupThere) {
   designMatrix <- data.frame(row.names  = rownames(dataset),
                             Group = as.factor(dataset[,c("Group")]),
                             stringsAsFactors = FALSE)
  sampleObject <- sample_data(designMatrix)
  physeqObjectNoTree = merge_phyloseq(fullPhySeqObj,sampleObject)
 }else{
   physeqObjectNoTree = fullPhySeqObj
 }

  pruneLevel <- param$representativeOTUs
  
  ### create, add trees, preprocess and prune phyloseq object 
  treeObject = rtree(ntaxa(physeqObjectNoTree), rooted=TRUE, tip.label=taxa_names(physeqObjectNoTree))
  physeqFullObject <- merge_phyloseq(physeqObjectNoTree,treeObject)
  physeqFullObject <- phyloSeqPreprocess(physeqFullObject)
  myTaxa = names(sort(taxa_sums(physeqFullObject), decreasing = TRUE)[1:pruneLevel])
  physeqFullObject <- prune_taxa(myTaxa,physeqFullObject)

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
