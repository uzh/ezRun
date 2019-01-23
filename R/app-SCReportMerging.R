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
                  appDefaults <<- rbind(scProtocol=ezFrame(Type="character", DefaultValue="", Description="Which single cell protocol?"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6, 
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        batchCorrection=ezFrame(Type="character", 
                                                                DefaultValue="CCA",
                                                                Description="Which batch correction method to use? None or CCA"),
                                        cc=ezFrame(Type="numeric",
                                                   DefaultValue=20,
                                                   Description="The number of CC for CCA analysis"),
                                        chosenClusters1=ezFrame(Type="charVector",
                                                                DefaultValue="",
                                                                Description="Clusters to choose to merge in sample 1"),
                                        chosenClusters2=ezFrame(Type="charVector",
                                                                DefaultValue="",
                                                                Description="Clusters to choose to merge in sample 2"),
                                        all2allMarkers=ezFrame(Type="logical",
                                                               DefaultValue=FALSE, 
                                                               Description="Run all against all cluster comparisons?"),
                                        markersToShow=ezFrame(Type="numeric", 
                                                              DefaultValue=10, 
                                                              Description="The markers to show in the heatmap of cluster marker genes"))
                }
              )
  )

ezMethodSCReportMerging = function(input=NA, output=NA, param=NA, 
                            htmlFile="00index.html"){
  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  if(input$getLength() != 2L){
    stop("It only works for merging two samples at once!")
  }
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  sce1URL <- input$getColumn("Static Report")[1]
  sce2URL <- input$getColumn("Static Report")[2]
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