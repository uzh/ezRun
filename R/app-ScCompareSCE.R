###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScCompareSCE <-
  setRefClass("EzAppScCompareSCE",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScCompareSCE
                  name <<- "EzAppScCompareSCE"
                }
              )
  )

ezMethodScCompareSCE = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(HDF5Array)
  library(SingleCellExperiment)

  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  sceURLs <- input$getColumn("Static Report")
  filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(scDataURLs)), "sce_h5")
  sce <- loadHDF5SummarizedExperiment(filePath)
  #subset the object to only contain the conditions we are interested in
  sce <- sce[,sce$Condition %in% c(param$sampleGroup, param$refGroup)]
  
  #Diff expression 
 diffGenes <- scranDiffGenes(sce)

 dataFiles = saveExternalFiles(list(differential_genes=diffGenes))
 saveRDS(input, "input.rds")
 saveRDS(param, "param.rds")
  
  makeRmdReport(dataFiles=dataFiles, rmdFile = "ScCompareSCE.Rmd", reportTitle = param$name) 
  return("Success")
  
}






