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
                  appDefaults <<- rbind(minCellsPerGene=ezFrame(Type="numeric", 
                                                                DefaultValue=5, 
                                                                Description="Minimum number of cells per gene for creating Seurat object"),
                                        minGenesPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=500, 
                                                                Description="Minimal number of genes per cell for Seurat filtering"),
                                        maxGenesPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=3000, 
                                                                Description="Maximal number of genes per cell for Seurat filtering"),
                                        maxMitoFraction=ezFrame(Type="numeric", 
                                                                DefaultValue=0.25, 
                                                                Description="Maximal fraction of mitochondrial reads per cell for Seurat filtering"),
                                        minReadsPerCell=ezFrame(Type="numeric", 
                                                                DefaultValue=5e4, 
                                                                Description="Minimal reads per cell of smart-Seq2 for Seurat filtering"),
                                        pcs=ezFrame(Type="numeric", 
                                                    DefaultValue=10, 
                                                    Description="The maximal dimensions to use for reduction"),
                                        pcGenes=ezFrame(Type="charVector", 
                                                        DefaultValue="", 
                                                        Description="The genes used in supvervised clustering"),
                                        x.low.cutoff=ezFrame(Type="numeric", 
                                                             DefaultValue=0.1, 
                                                             Description="Bottom cutoff on x-axis for identifying variable genes"),
                                        x.high.cutoff=ezFrame(Type="numeric", 
                                                              DefaultValue=8, 
                                                              Description="Top cutoff on x-axis for identifying variable genes"),
                                        y.cutoff=ezFrame(Type="numeric", 
                                                         DefaultValue=1, 
                                                         Description="Bottom cutoff on y-axis for identifying variable genes"),
                                        vars.to.regress=ezFrame(Type="charVector", 
                                                                DefaultValue="nUMI,perc_mito", 
                                                                Description="Variables to regress out"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.8, 
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
  
  sce <- list(sce1=sce1, sce2=sce2)
  
  ## Merge samples without batch correction
  # scData1 <- metadata(sce1)$scData
  # scData2 <- metadata(sce2)$scData
  
}