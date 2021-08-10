###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppspatialSeurat <-
  setRefClass("EzAppspatialSeurat",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodspatialSeurat
                  name <<- "EzAppspatialSeurat"
                  appDefaults <<- rbind(scProtocol=ezFrame(Type="character", DefaultValue="10X", Description="Which single cell protocol?"),
                                        minReadsPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=5e4, 
                                                                Description="Minimal reads per cell of smart-Seq2 for Seurat filtering"),
                                        npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=20,
                                                    Description="The maximal dimensions to use for reduction"),
                                        pcGenes=ezFrame(Type="charVector", 
                                                        DefaultValue="", 
                                                        Description="The genes used in supvervised clustering"),
                                        SCT.regress=ezFrame(Type="character", 
                                                           DefaultValue="none", 
                                                           Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."),
                                        DE.method=ezFrame(Type="charVector", 
                                                          DefaultValue="wilcoxon", 
                                                          Description="Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.5,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        all2allMarkers=ezFrame(Type="logical", 
                                                               DefaultValue=FALSE, 
                                                               Description="Run all against all cluster comparisons?"),
                                        cellsFraction=ezFrame(Type="numeric", 
                                                                DefaultValue=0.05, 
                                                                Description="A gene will be kept if it is expressed in at least this percentage of cells"),
                                        nUMIs=ezFrame(Type="numeric", 
                                                      DefaultValue=1, 
                                                      Description='A gene will be kept if it has at least nUMIs in the fraction of cells specified before'),
                                        nmad=ezFrame(Type="numeric", 
                                                     DefaultValue=3, 
                                                     Description="Median absolute deviation (MAD) from the median value of each metric across all cells")
                                        )
                }
              )
  )

ezMethodspatialSeurat <- function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  scData <- Load10X_Spatial(input, param)
  scData <- seuratClusteringV3(scData, param, assay="SCT")
  
  pvalue_allMarkers <- 0.05
  
  #positive cluster markers
  clusterMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  
  #spatially variable genes
  spatialMarkers <- spatialMarkers(scData)
  
  #Save some results in external files 
  dataFiles = saveExternalFiles(list(cluster_markers=clusterMarkers, spatial_markers=spatialMarkers))
  saveRDS(scData, "scData.rds")
  
  
  makeRmdReport(dataFiles=dataFiles, rmdFile = "spatialSeurat.Rmd", reportTitle = metadata(sce)$param$name) 
  return("Success")
}


