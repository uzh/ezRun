###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMetagenomeAnnotationQC = function(input=NA, output=NA, param=NA, 
                                        htmlFile="00index.html"){
  
  library(purrr)
  library(rtracklayer)
  library(ggplot2)
  library(RColorBrewer)
  library(GO.db)
  library(iNEXT)
  library(pheatmap)
  library(cowplot)
  library(DESeq2)
  library(plyr)
  library(dplyr)
  library(tidyr)
  dataset = input$meta
  sampleNames = input$getNames() 
  numberOfTopNCategories = param$numberOfTopNCategories
  isGroupThere = param$grouping != ""
  group = param$grouping
  sampleGroup = param$sampleGroup
  refGroup = param$refGroup
  areThereRepl  <- nrow(dataset) > 3
    colnames(dataset) <-  gsub(" \\[File\\]","",colnames(dataset))
  allColumns <- dataset[,c("prodigalPredictionFile","interproscanFile","binSummaryFile")]
  plotLabels <- input$getNames()
  ## Copy all files locally
  copyLoopOverFiles <- function(x){ 
    lapply(x,function(x) ezSystem(paste("cp",x,".")))
  }
  listOfListAllFilesTemp <- list(input$getFullPaths("prodigalPredictionFile"),
                                    input$getFullPaths("interproscanFile"),
                                    input$getFullPaths("binSummaryFile"))
  lapply(listOfListAllFilesTemp,copyLoopOverFiles)
  listOfListAllFiles <- as.list(data.frame(apply(allColumns,2,basename),stringsAsFactors = F))
  namedList <- lapply(listOfListAllFiles,function(x){
    y=as.list(x);names(y)=sampleNames;return(y)
    })
  
  ## construct final lists
  ## prodigal
  k=0
  fullSumm <- list()
  subsetDataToPartial00 <-list()
  for (file in namedList$prodigalPredictionFile){
    k=k+1
    method <- sampleNames[k]
    fullSumm[[method]] <-  prodigalFileReport(file,method)$fullSumm
    subsetDataToPartial00[[method]] <- prodigalFileReport(file,method)$subsetDataToPartial00
  }

  fullSummPlot <-  do.call("rbind",fullSumm)
  subsetDataToPartial00Plot <- do.call("rbind",subsetDataToPartial00)

  ## IPS
  interproscanList <- mapply(interproscanFileReport,
                             namedList$interproscanFile,sampleNames,
                             N = numberOfTopNCategories, SIMPLIFY = FALSE)
  interproscanListForWrap <- lapply(names(interproscanList),function(x) {
    scores <- interproscanList[[x]][["summDF"]][["score"]]
    nProtId <- length(scores)
    data.frame(score = scores,
    sample = paste0(x," (N=",nProtId,")"))
  })
  interproscanListForPlot <- do.call("rbind",interproscanListForWrap)
  
  
  ### IPS for heatmaps 
  prepareIPSforHeatmap <- function(feature,colname){
  listForGO <- lapply(interproscanList,function(x) x[[feature]])
  listForGoForHeatmap <- lapply(names(listForGO),function(x) {
    z <- data.frame(listForGO[[x]][["abundance"]], row.names = listForGO[[x]][[colname]])
    names(z) <- x
    return(z)
  })
  finalDF <- listOfAbundMerge(listForGoForHeatmap, names=names(interproscanList))
  return(finalDF)
  }
  ### generate subset DF for heatmap
  GOTermMatrixForHeatmap <- prepareIPSforHeatmap("topN_GO","GOterm")
  descTermMatrixForHeatmap <- prepareIPSforHeatmap("topN_desc","description")
  ### generate full DF for diff. expr.
  GOTermMatrixForDiffExpr <- prepareIPSforHeatmap("full_GO","GOterm")
  descTermMatrixForDiffExpr <- prepareIPSforHeatmap("full_descrip","description")
  
  ## Summary bins
  mergedSummaryBinDFList <- lapply(namedList$binSummaryFile,mergeSummaryBinFiles)
  mergedSummaryBinDF <- do.call("rbind",mergedSummaryBinDFList)
  
  ##
  setwdNew(basename(output$getColumn("Report")))
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "metagenomeAnnotation.Rmd", 
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="metagenomeAnnotation.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
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
                  appDefaults <<- rbind(numberOfTopNCategories = ezFrame(Type="integer",  DefaultValue="20",
                                                                         Description="How many top N GO and 
                                                                         prot families.")
                  )
                  }
              )
  )







