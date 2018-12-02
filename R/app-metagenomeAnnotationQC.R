###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMetagenomeAnnotationQC = function(input=NA, output=NA, param=NA, 
                                        htmlFile="00index.html"){
  
  library(Biostrings)
  library(rtracklayer)
  library(cowplot)
  library(kableExtra)
  library(plyr)
  library(dplyr)
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(ape)
  require(ggplot2)
  library(scales)
  library(RColorBrewer)
  dataset = input$meta
  sampleName = input$getNames() 
  numberOfTopNCategories = param$numberOfTopNCategories
  ## get input files
  prodigalGffFile <- input$getFullPaths("prodigalPredictionFile")
  IPSGffFile <- input$getFullPaths("interproscanFile")
  krakenFile <- input$getFullPaths("krakenLabelsFile")
  
  ezSystem(paste("cp", prodigalGffFile, basename(prodigalGffFile)))
  ezSystem(paste("cp", IPSGffFile, basename(IPSGffFile)))
  ezSystem(paste("cp", krakenFile, basename(krakenFile)))
  
  prodigalGffImport <- import.gff(prodigalGffFile)
  prodigalSummaryDF <- data.frame(mcols(prodigalGffImport), stringsAsFactors = F)
  prodigalSummaryDF$gc_cont <- as.numeric(prodigalSummaryDF$gc_cont)
  prodigalSummaryDF$conf <- as.numeric(prodigalSummaryDF$conf)
  
  IPSGffImport <- import.gff(IPSGffFile)
  description <- mcols(IPSGffImport)$signature_desc
  description[sapply(description,function(x) length(x)==0)] <- "NA"
  description <- sapply(description,function(x)x[1])
  ontology <- mcols(IPSGffImport)$Ontology_term
  ontology[sapply(ontology,function(x) length(x)==0)] <- "NA" 
  ontology <- sapply(ontology,function(x)x[1])
  IPSGffSummaryDF <- data.frame(score = as.numeric(mcols(IPSGffImport)$score),
                                description = description, 
                                GOterm = ontology,
                                   type = mcols(IPSGffImport)$type,
                                stringsAsFactors = F)
  IPSGffSummaryDF <- IPSGffSummaryDF[IPSGffSummaryDF$type == "protein_match",
                                     c("score","description","GOterm")]
  ## extract N entries with top frequency 
  extractTopN <- function(DF,column,N){
    col <- vector()
    tabNoNa <- DF[DF[[column]] != "NA",]
    tab <- table(tabNoNa[[column]])
  tab_s <- sort(tab)                                           
  col <- data.frame(tail(names(tab_s), N), stringsAsFactors = F)
  colnames(col) <- column
  topN <- data.frame(cbind(col, abundance = tail(as.data.frame(tab_s)$Freq, N)),
                     stringsAsFactors = F)
  topN <- topN[order(topN$abundance), ]
  topN[[column]] <- gsub("\"","",topN[[column]])
  return(topN)
  }
}
##' @template app-template
##' @templateVar method ezMethodMetagenomeAnnotationQC()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMetagenomeAnnotationQC<-
  setRefClass("EzAppMetagenomeAnnotationQC",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMetagenomeAnnotationQC
                  name <<- "EzAppMetagenomeAnnotationQC"
                  appDefaults <<- rbind(cutOffTaxonomy = ezFrame(Type="integer",  DefaultValue="80",Description="Cut-off for taxonomy assignment")
                  )
                  }
              )
  )







