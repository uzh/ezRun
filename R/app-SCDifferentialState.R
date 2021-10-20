###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCDifferentialState <-
  setRefClass("EzAppSCDifferentialState",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCDifferentialState
                  name <<- "EzAppSCDifferentialState"
                  appDefaults <<- rbind(
  DE.method=ezFrame(Type="charVector", 
                                                          DefaultValue="wilcox", 
                                                          Description="Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"),
                                        DE.regress=ezFrame(Type="charVector", 
                                                           DefaultValue="Batch", 
                                                           Description="Variables to regress out if the test LR is chosen"))
                }
              )
  )

ezMethodSCDifferentialState = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(HDF5Array)
  library(SingleCellExperiment)

  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  sceURLs <- input$getColumn("Static Report")
  filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(sceURLs)), "sce_h5")
  filePath_course <- file.path(paste0("/srv/GT/analysis/course_sushi/public/projects/", input$getColumn("Report"), "/sce_h5"))
  
  #In case it is in hdf5 format the Seurat object is not stored in the metadata slot, so we have to build it and store it there, since the 
  #clustering functions below work with sce objects and take the seurat object from them.
  if(file.exists(filePath)) 
    filePath <- filePath
  else
    filePath <- filePath_course
  
  sce <- loadHDF5SummarizedExperiment(filePath)
  #create a seurat object from the sce object using the raw counts (I can't use the function sce_to_seurat() because it can't deal with delayed matrices)
  scData <- CreateSeuratObject(counts=counts(sce),meta.data=data.frame(colData(sce))) 
  #subset the object to only contain the conditions we are interested in
  Idents(scData) <- scData$Condition
  scData <- subset(scData, idents=c(param$sampleGroup, param$refGroup))
  
  pvalue_allMarkers <- 0.05
  
  #Before calculating the conserved markers and differentially expressed genes across conditions I will discard the clusters that were too small in at least one group
  clusters_freq <- data.frame(table(scData@meta.data[,c("Condition","ident")]))
  small_clusters <- ""
  small_clusters <- unique(as.character(clusters_freq[clusters_freq[,"Freq"] < 10, "ident"]))
  
  diffGenes <- NULL
  consMarkers <- NULL
  
  
  #only subset object and calculate diff genes and conserved markers if there is at least one cluster shared among conditions
  if (!all(scData$seurat_clusters %in% small_clusters)) {
     Idents(scData) <- scData$seurat_clusters
     scData <- subset(scData, idents = small_clusters, invert = TRUE)
     #conserved cluster markers
     consMarkers <- conservedMarkers(scData)
     #differentially expressed genes between clusters and conditions (in case of several conditions)
     diffGenes <- diffExpressedGenes(scData)
  }
  dataFiles = saveExternalFiles(list(differential_genes=diffGenes, conserved_markers=consMarkers, current_cells=data.frame(cell=colnames(scData))))
  saveRDS(input, "input.rds")
  saveRDS(param, "param.rds")
  
  makeRmdReport(dataFiles=dataFiles, rmdFile = "SCDifferentialState.Rmd", reportTitle = metadata(sce)$param$name) 
  return("Success")
  
}






