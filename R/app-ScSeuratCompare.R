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
  library(BiocParallel)
  
  
  if (param$cores > 1){
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    ## scDblFinder fails with many cells and MulticoreParam
    BPPARAM <- SerialParam() 
  }
  register(BPPARAM)
  require(future)
  plan("multicore", workers = param$cores)
  set.seed(38)
  future.seed = TRUE
  options(future.rng.onMisuse="ignore")
  options(future.globals.maxSize = param$ram*1024^3)
  
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  reportDir <- input$getFullPaths("Report")
  scData <- readRDS(file.path(reportDir, "scData.rds"))
  
  DefaultAssay(scData) = "SCT" 
  #subset the object to only contain the conditions we are interested in
  Idents(scData) <- scData[[param$grouping]]
  stopifnot(c(param$sampleGroup, param$refGroup) %in% Idents(scData))
  scData <- subset(scData, idents=c(param$sampleGroup, param$refGroup))
  
  pvalue_allMarkers <- 0.05
  
  #Before calculating the conserved markers and differentially expressed genes across conditions I will discard the clusters that were too small in at least one group
  Idents(scData) <- scData[[param$CellIdentity]]
  clusters_freq <- table(grouping=scData@meta.data[[param$grouping]], cellIdent=Idents(scData)) %>% data.frame()
  small_clusters <- clusters_freq[clusters_freq$Freq < 10, "cellIdent"] %>% as.character() %>% unique()
  big_clusters <- setdiff(Idents(scData), small_clusters)
  
  if(length(slot(scData[['SCT']], "SCTModel.list")) > 2){
    toKeep <- which(sapply(SCTResults(scData[['SCT']], slot = "cell.attributes"), nrow) != 0)
    slot(scData[['SCT']], "SCTModel.list") = slot(scData[['SCT']], "SCTModel.list")[toKeep]
  }
  scData <- PrepSCTFindMarkers(scData)
  consMarkers <- conservedMarkers(scData, grouping.var = param$grouping)
  #differentially expressed genes between clusters and conditions (in case of several conditions)
  diffGenes <- diffExpressedGenes(scData, param, grouping.var = param$grouping)
  
  writexl::write_xlsx(consMarkers, path="consMarkers.xlsx")
  writexl::write_xlsx(diffGenes, path="diffGenes.xlsx")
  
  makeRmdReport(param=param, output=output, scData=scData,
                rmdFile = "ScSeuratCompare.Rmd", reportTitle = paste0(param$name))
  return("Success")
}
