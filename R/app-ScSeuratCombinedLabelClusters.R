###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratCombinedLabelClusters <-
  setRefClass("EzAppScSeuratCombinedLabelClusters",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScSeuratCombinedLabelClusters
                  name <<- "EzAppScSeuratCombinedLabelClusters"
                  appDefaults <<- rbind(enrichrDatabase=ezFrame(Type = "charVector", 
                                                                DefaultValue = "", 
                                                                Description="enrichR databases to search"),
                                        computePathwayTFActivity=ezFrame(Type="logical", 
                                                                 DefaultValue="TRUE",
                                                                 Description="Whether we should compute pathway and TF activities."),
                                        DE.method=ezFrame(
                                          Type="charVector", 
                                          DefaultValue="wilcox", 
                                          Description="Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"),
                                        DE.regress=ezFrame(
                                          Type="charVector", 
                                          DefaultValue="Batch", 
                                          Description="Variables to regress out if the test LR is chosen"),
                                        min.pct = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = 0.1,
                                          Description = "Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations."
                                        ),
                                        logfc.threshold = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = 0.25,
                                          Description = "Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
                                        ))
                }
              )
  )

ezMethodScSeuratCombinedLabelClusters = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(rlist)
  library(HDF5Array)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(AUCell)
  library(enrichR)
  library(decoupleR)
  library(Azimuth)
  library(qs2)
  library(BiocParallel)
  
  BPPARAM <- MulticoreParam(workers = param$cores)
  register(BPPARAM)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  #the individual sce objects can be in hdf5 format (for new reports) or in rds format (for old reports)
  filePath <- input$getFullPaths("SeuratObject")
  filePath_course <- file.path("/srv/GT/analysis/course_sushi/public/projects", input$getColumn("SeuratObject"))
  
  if (!file.exists(filePath[1])) {
    filePath <- filePath_course
  }
  names(filePath) <- input$getNames()
  stopifnot("App only supports single integrated dataset!" = length(input$getNames()) == 1)
  
  # Load previous dataset
  scData <- ezLoadRobj(filePath, nthreads=param$cores)
  oldParams <- ezLoadRobj(file.path(input$getFullPaths("Report"), "param.rds"))
  param <- ezUpdateMissingParam(param, oldParams)
  param$refBuild <- oldParams$refBuild

  # load cluster annotation file
  clusterAnnoFn <- file.path(param$dataRoot, param$ClusterAnnotationFile)
  if (ezIsSpecified(param$ClusterAnnotationFile)) {
    clusterAnnoFn <- file.path(param$dataRoot, param$ClusterAnnotationFile)
    stopifnot("The cluster annotation file does not exist or is not an .xlsx file!" = 
                file.exists(clusterAnnoFn) && str_ends(clusterAnnoFn, ".xlsx$"))
  } else {
    stop("Must supply cluster annotation file path.", call.=FALSE)
  }
  clusterAnno <- readxl::read_xlsx(clusterAnnoFn) %>% 
    as_tibble() %>%
    dplyr::select(1:3) %>% # remove all other columns
    dplyr::rename(c("_"=1, "Cluster"=2, "ClusterLabel"=3)) # we don't use the first column
  labelMap <- as.character(clusterAnno$ClusterLabel)
  names(labelMap) <- as.character(clusterAnno$Cluster)
  
  # Do the renaming
  scData$cellTypeIntegrated <- unname(labelMap[as.character(Idents(scData))])
  Idents(scData) <- scData$cellTypeIntegrated
  scData$ident <- Idents(scData)
  
  # perform all of the analysis
  anno <- getSeuratMarkersAndAnnotate(scData, param, BPPARAM = BPPARAM)
  
  # save the markers
  writexl::write_xlsx(anno$markers, path="posMarkers.xlsx")
  qs2::qs_save(scData, "scData.qs2", nthreads = param$cores)
  
  # Save some results in external files
  reportTitle <- 'SCReport - MultipleSamples based on Seurat'
  makeRmdReport(param=param, output=output, scData=scData, 
                enrichRout=anno$enrichRout, TFActivity=anno$TFActivity, 
                pathwayActivity=anno$pathwayActivity, aziResults=anno$aziResults,
                cells.AUC=anno$cells.AUC, singler.results=anno$singler.results,
                rmdFile = "ScSeuratCombine.Rmd", reportTitle = reportTitle) 
  return("Success")
}
