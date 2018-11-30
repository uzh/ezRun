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
  dataset = input$meta
  fileNames <- input$getNames()
  ### Analyzes results with phyloseq: preparing objects to be processed in the Rmd file
  
  relevantColumns <- c("OTUsCountTable","OTUsToTaxonomyFile")
  colnames(dataset) <-  gsub(" \\[File\\]","",colnames(dataset))
  allColumns <- dataset[,relevantColumns]
  ## Copy all files locally
  copyLoopOverFiles <- function(x){ 
    lapply(x,function(x) ezSystem(paste("cp",file.path(DEMO_DATA_ROOT,x),"./")))
  }
  listOfListAllFiles <- as.list(allColumns)
  lapply(listOfListAllFiles,copyLoopOverFiles)
  
  ### merge all count files into one file
  k=0
  OTUsCount <- list()
  OTUsCountNoLabel <- list()
  Group <- vector()
    for (file in listOfListAllFiles["OTUsCountTable"]$OTUsCountTable){
    k=k+1
    rawFile <- basename(file)
    OTUsCount[[k]] <- read.delim(rawFile, header = T,stringsAsFactors = FALSE)
    colnames(OTUsCount[[k]]) <- gsub("Otu[0-9]",paste0(fileNames[k],"_Otu"),colnames(OTUsCount[[k]]))
    relCols <- grep("Otu[0-9]",colnames(OTUsCount[[k]]),value = T)
    OTUsCountNoLabel[[k]] <- as.matrix(OTUsCount[[k]][,relCols])
    Group[k] <- OTUsCount[[k]]$Group
  }
  
  fullOTUCountTable <- cbind(Group,data.frame(do.call("adiag",OTUsCountNoLabel),
                                              stringsAsFactors = F))
  ### create phyloseq OTU object
  otuObject <- phyloSeqOTU(fullOTUCountTable)
  
  ### merge all taxa files into one file  
  k=0
  OTUsToTaxonomyDF <- list()
  for (file in listOfListAllFiles["OTUsToTaxonomyFile"]$OTUsToTaxonomyFile){
    k=k+1
    rawFile <- basename(file)
    OTUsToTaxonomyDF[[k]] <- read.delim(rawFile, header = T,stringsAsFactors = FALSE)
    OTUsToTaxonomyDF[[k]]$OTU <- paste(fileNames[k],OTUsToTaxonomyDF[[k]]$OTU,sep = "_")
  }
  fullTaxaTable <- do.call("rbind",OTUsToTaxonomyDF)

  ### create phyloseq taxa object
  taxaObject <- phyloSeqTaxa(fullTaxaTable)

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
