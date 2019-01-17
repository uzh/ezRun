###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCReportMerging <-
  setRefClass("EzAppSCReportMerging",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCReportMerging
                  name <<- "EzAppSCReportMerging"
                  appDefaults <<- rbind(resolution=ezFrame(Type="numeric", 
                                                           DefaultValue="", 
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        markersToShow=ezFrame(Type="numeric", 
                                                              DefaultValue=10, 
                                                              Description="The markers to show in the heatmap of cluster marker genes"),
                                        markersToCheck=ezFrame(Type="charList", 
                                                               DefaultValue="", 
                                                               Description="The markers to check"),
                                        runPseudoTime=ezFrame(Type="logical", 
                                                              DefaultValue=FALSE, 
                                                              Description="Run PseudoTime for single cell data?"),
                                        all2allMarkers=ezFrame(Type="logical", 
                                                               DefaultValue=FALSE, 
                                                               Description="Run all against all cluster comparisons?"),
                                        batchCorrection=ezFrame(Type="character", 
                                                                DefaultValue="None",
                                                                Description="Which batch correction method to use?"),
                                        chosenClusters1=ezFrame(Type="charVector", 
                                                                DefaultValue="",
                                                                Description="Clusters to choose to merge in sample 1"),
                                        chosenClusters2=ezFrame(Type="charVector", 
                                                                DefaultValue="",
                                                                Description="Clusters to choose to merge in sample 2"))
                }
              )
  )

ezMethodSCReportMerging = function(input=NA, output=NA, param=NA, 
                            htmlFile="00index.html"){
  if(input$getLength() != 2L){
    stop("It only works for merging two samples at once!")
  }
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  sce1 <- readRDS(file.path(input$getFullPaths("Report"), "sce.rds")[1])
  sce2 <- readRDS(file.path(input$getFullPaths("Report"), "sce.rds")[2])
  
  if(ezIsSpecified(param$chosenClusters1)){
    chosenCells1 <- names(metadata(sce1)$scData@ident)[metadata(sce1)$scData@ident %in% 
                                                         param$chosenClusters1]
    sce1 <- sce1[, chosenCells1]
  }
  if(ezIsSpecified(param$chosenClusters2)){
    chosenCells2 <- names(metadata(sce2)$scData@ident)[metadata(sce2)$scData@ident %in% 
                                                         param$chosenClusters2]
    sce2 <- sce2[, chosenCells2]
  }
  ans <- list(sce1=sce1, sce2=sce2, param=param)
  saveRDS(ans, file="ans.rds")
  
}