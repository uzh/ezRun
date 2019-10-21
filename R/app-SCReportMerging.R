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
                  appDefaults <<- rbind(scProtocol=ezFrame(Type="character", DefaultValue="10x", Description="Which single cell protocol?"),
                                        x.low.cutoff=ezFrame(Type="numeric", 
                                                             DefaultValue=0.0125, 
                                                             Description="Bottom cutoff on x-axis for identifying variable genes"),
                                        x.high.cutoff=ezFrame(Type="numeric", 
                                                              DefaultValue=3,
                                                              Description="Top cutoff on x-axis for identifying variable genes"),
                                        y.cutoff=ezFrame(Type="numeric", 
                                                         DefaultValue=0.5,
                                                         Description="Bottom cutoff on y-axis for identifying variable genes"),
                                        vars.to.regress=ezFrame(Type="charVector", 
                                                                DefaultValue="nUMI,perc_mito", 
                                                                Description="Variables to regress out"),
                                        pcs=ezFrame(Type="numeric", 
                                                    DefaultValue=20, 
                                                    Description="The maximal dimensions to use for reduction"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        batchCorrection=ezFrame(Type="character", 
                                                                DefaultValue="CCA",
                                                                Description="Which batch correction method to use? None or CCA"),
                                        cc=ezFrame(Type="numeric",
                                                   DefaultValue=20,
                                                   Description="The number of CC for CCA analysis"),
                                        resolutionCCA=ezFrame(Type="numeric", 
                                                              DefaultValue=0.6,
                                                              Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        chosenClusters=ezFrame(Type="charList",
                                                               DefaultValue="",
                                                               Description="The clusters to choose from each sample.In the format of sample1=cluster1,cluster2;sample2=cluster1,cluster2."),
                                        all2allMarkers=ezFrame(Type="logical",
                                                               DefaultValue=FALSE, 
                                                               Description="Run all against all cluster comparisons?"),
                                        computeDifferentialExpression=ezFrame(Type="logical",
                                                               DefaultValue=TRUE, 
                                                               Description="Compute differential expression between samples"),
                                        markersToShow=ezFrame(Type="numeric", 
                                                              DefaultValue=10, 
                                                              Description="The markers to show in the heatmap of cluster marker genes"),
                                        maxSamplesSupported=ezFrame(Type="numeric", 
                                                              DefaultValue=5, 
                                                              Description="Maximum number of samples to compare"))
                }
              )
  )

ezMethodSCReportMerging = function(input=NA, output=NA, param=NA, 
                            htmlFile="00index.html"){
  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  if(input$getLength() > param$maxSamplesSupported){
    stop(paste("It only works for", param$maxSamplesSupported, "at most"))
  }
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  param$name <- paste(param$name, paste(input$getNames(), collapse=", "),
                      sep=": ")
  
  sceURLs <- input$getColumn("Static Report")
  
  saveRDS(param, file = "param.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCReportMerging.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCReportMerging.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  return("Success")
  # sce1 <- readRDS(file.path(input$getFullPaths("Report"), "sce.rds")[1])
  # sce2 <- readRDS(file.path(input$getFullPaths("Report"), "sce.rds")[2])
  # 
  # if(ezIsSpecified(param$chosenClusters1)){
  #   chosenCells1 <- names(metadata(sce1)$scData@ident)[metadata(sce1)$scData@ident %in% 
  #                                                        param$chosenClusters1]
  #   sce1 <- sce1[, chosenCells1]
  # }
  # if(ezIsSpecified(param$chosenClusters2)){
  #   chosenCells2 <- names(metadata(sce2)$scData@ident)[metadata(sce2)$scData@ident %in% 
  #                                                        param$chosenClusters2]
  #   sce2 <- sce2[, chosenCells2]
  # }
  # ans <- list(sce1=sce1, sce2=sce2, param=param)
  # saveRDS(ans, file="ans.rds")
}