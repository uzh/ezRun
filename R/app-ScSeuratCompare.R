###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratCompare <-
  setRefClass("EzAppScSeuratCompare",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScSeuratCompare
                  name <<- "EzAppScSeuratCompare"
                  appDefaults <<- rbind(DE.method=ezFrame(Type="charVector", DefaultValue="wilcox", 
                                                          Description="Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"),
                                        DE.regress=ezFrame(Type="charVector", DefaultValue="Batch", Description="Variables to regress out if the test LR is chosen"))
                }
              )
  )

ezMethodScSeuratCompare = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(HDF5Array)
  library(SingleCellExperiment)

  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  scDataURLs <- input$getColumn("Static Report")
  filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(scDataURLs)), "scData.rds")
  filePath_course <- file.path("/srv/GT/analysis/course_sushi/public/projects", dirname(scDataURLs), "scData.rds")
  
  if(!file.exists(filePath)) 
    filePath <- filePath_course
  
  scData <- readRDS(filePath)
  
  DefaultAssay(scData) = "SCT" 
  #subset the object to only contain the conditions we are interested in
  Idents(scData) <- scData[[param$grouping]]
  scData <- subset(scData, idents=c(param$sampleGroup, param$refGroup))
  
  pvalue_allMarkers <- 0.05
  
  #Before calculating the conserved markers and differentially expressed genes across conditions I will discard the clusters that were too small in at least one group
  Idents(scData) <- scData[[param$CellIdentity]]
  clusters_freq <- data.frame(table(scData@meta.data[[param$grouping]], Idents(scData)))
  small_clusters <- ""
  small_clusters <- unique(as.character(clusters_freq[clusters_freq[,"Freq"] < 10, 2]))
  
  diffGenes <- NULL
  consMarkers <- NULL
  
  #only subset object and calculate diff genes and conserved markers if there is at least one cluster shared among conditions
  if (!all(Idents(scData) %in% small_clusters)) {
    scData <- subset(scData, idents = small_clusters, invert = TRUE)
    
    ###PrepSCTFindMarkers crashes for empty models which occur if the data consists of more than 2 samples/conditions 
    if(length(slot(object = scData[['SCT']], name = "SCTModel.list")) > 2){
      toKeep <- which(sapply(SCTResults(object = scData[['SCT']], slot = "cell.attributes"), nrow) != 0)
      slot(object = scData[['SCT']], name = "SCTModel.list") = slot(object = scData[['SCT']], name = "SCTModel.list")[toKeep]
    }
    #conserved cluster markers
    scData <- PrepSCTFindMarkers(scData)
    consMarkers <- conservedMarkers(scData, grouping.var = param$grouping)
    #differentially expressed genes between clusters and conditions (in case of several conditions)
    diffGenes <- diffExpressedGenes(scData, param, grouping.var = param$grouping)
  }
  dataFiles = saveExternalFiles(list(differential_genes=diffGenes, conserved_markers=consMarkers))
  saveRDS(input, "input.rds")
  saveRDS(param, "param.rds")
  saveRDS(scData, file = "scData.rds")
  
  makeRmdReport(dataFiles=dataFiles, rmdFile = "ScSeuratCompare.Rmd", reportTitle = param$name) 
  return("Success")
}
