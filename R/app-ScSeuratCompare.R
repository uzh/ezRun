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
                  appDefaults <<- rbind(
                    DE.method=ezFrame(Type="charVector", DefaultValue="wilcox", 
                                      Description="Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"),
                    DE.regress=ezFrame(Type="charVector", DefaultValue="Batch", Description="Variables to regress out if the test LR is chosen"),
                    sccomp.variability=ezFrame(Type="logical", DefaultValue="FALSE", Description="Whether to test for differential variability in sccomp"))
                }
              )
  )

ezMethodScSeuratCompare = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(HDF5Array)
  library(SingleCellExperiment)
  library(qs2)
  library(tidyverse)
  library(cmdstanr)
  
  ## Setup sccomp
  # Create scratch directory
  scratch_dir <- "/scratch/sccomp_output"
  dir.create(scratch_dir, recursive = TRUE, mode = "0777", showWarnings = FALSE)
  
  # Load sccomp and set up cmdstan
  library(sccomp, lib.loc = '/srv/GT/databases/writable_R_package')
  # cmdstanr::set_cmdstan_path("/misc/ngseq12/packages/Dev/R/4.4.2/lib/R/cmdstan-2.36.0")
  
  # check if in pseudobulk mode
  pseudoBulkMode <- ezIsSpecified(param$replicateGrouping) && param$pseudoBulkMode == "true"
  
  # Load model
  sccomp:::load_model(
    name = "glm_multi_beta_binomial",
    threads = 4,
    cache_dir = '/srv/GT/databases/writable_R_package/sccomp/'
  )
  
  
  ###
  set.seed(38)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  reportDir <- input$getFullPaths("Report")
  
  # First try qs2 format
  message("Attempting to load Seurat data...")
  qs_path <- file.path(reportDir, "scData.qs2")
  rds_path <- file.path(reportDir, "scData.rds")
  
  if(file.exists(qs_path)) {
    message("Found .qs2 file, attempting to load...")
    scData <- qs_read(qs_path, nthreads=param$cores)
  } else if(file.exists(rds_path)) {
    message("Found .rds file, attempting to load...")
    scData <- readRDS(rds_path)
  } else {
    stop("Neither scData.qs2 nor scData.rds found in ", reportDir)
  }
  
  DefaultAssay(scData) = "SCT" 
  #subset the object to only contain the conditions we are interested in
  Idents(scData) <- scData@meta.data[[param$grouping]]
  stopifnot(c(param$sampleGroup, param$refGroup) %in% Idents(scData))
  scData <- subset(scData, idents=c(param$sampleGroup, param$refGroup))
  
  # Only run sccomp if 'Sample' metadata exists
  if (pseudoBulkMode) {
    # Prepare data for sccomp using cellTypeIntegrated
    metadata <- scData@meta.data
    cell_counts <- table(metadata[[param$CellIdentity]], metadata[[param$replicateGrouping]]) %>%
      as.data.frame() %>%
      rename(cell_group = Var1, sample = Freq)
    
    # Run sccomp analysis with condition
    sccomp_res <- scData %>%
      sccomp_estimate(
        formula_composition = as.formula(paste("~", param$grouping)),  
        .sample = !!sym(param$replicateGrouping),   
        .cell_group = !!sym(param$CellIdentity),   
        cores = as.integer(param$cores), 
        output_directory = scratch_dir,
        verbose = TRUE
      )
    
    sccomp_res <- sccomp_res %>%
      sccomp_remove_outliers(cores = as.integer(param$cores)) %>%
      sccomp_test()
    
    # Save sccomp results
    saveRDS(sccomp_res, "sccomp_results.rds")
  } else {
    message("'Sample' metadata not found. Skipping sccomp analysis.")
  }
  
  
  pvalue_allMarkers <- 0.05
  
  #Before calculating the conserved markers and differentially expressed genes across conditions I will discard the clusters that were too small in at least one group
  Idents(scData) <- scData@meta.data[[param$CellIdentity]]
  clusters_freq <- table(grouping=scData@meta.data[[param$grouping]], cellIdent=Idents(scData)) %>% data.frame()
  small_clusters <- clusters_freq[clusters_freq$Freq < 10, "cellIdent"] %>% as.character() %>% unique()
  big_clusters <- setdiff(Idents(scData), small_clusters)
  
  if (length(slot(scData[['SCT']], "SCTModel.list")) > 2) {
    toKeep <- which(sapply(SCTResults(scData[['SCT']], slot = "cell.attributes"), nrow) != 0)
    slot(scData[['SCT']], "SCTModel.list") = slot(scData[['SCT']], "SCTModel.list")[toKeep]
  }
  if (pseudoBulkMode) {
    scData_agg <- AggregateExpression(scData, assays = "RNA", 
                                      return.seurat = TRUE, 
                                      group.by = c(param$grouping, param$replicateGrouping, param$CellIdentity))
    Idents(scData_agg) <- scData_agg@meta.data[[param$CellIdentity]]
    consMarkers <- conservedMarkers(scData_agg, grouping.var = param$grouping, 
                                    pseudoBulkMode = pseudoBulkMode)
    diffGenes <- diffExpressedGenes(scData_agg, param, grouping.var = param$grouping)
  } else {
    scData <- PrepSCTFindMarkers(scData)
    consMarkers <- conservedMarkers(scData, grouping.var = param$grouping, 
                                    pseudoBulkMode = pseudoBulkMode)
    diffGenes <- diffExpressedGenes(scData, param, grouping.var = param$grouping)
  }
  
  # Save the files for the report
  writexl::write_xlsx(consMarkers, path="consMarkers.xlsx")
  writexl::write_xlsx(diffGenes, path="diffGenes.xlsx")
  qs_save(scData, "scData.qs2", nthreads = as.integer(param$cores))
  makeRmdReport(param=param, output=output, scData=scData,
                rmdFile = "ScSeuratCompare.Rmd", reportTitle = paste0(param$name))
  return("Success")
}
