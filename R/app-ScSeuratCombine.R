###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeuratCombine <-
  setRefClass("EzAppScSeuratCombine",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodScSeuratCombine
                  name <<- "EzAppScSeuratCombine"
                  appDefaults <<- rbind(npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=30, 
                                                    Description="The maximal dimensions to use for reduction"),
                                        pcGenes = ezFrame(
                                          Type = "charVector",
                                          DefaultValue = "",
                                          Description = "The genes used in supvervised clustering"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        batchCorrection=ezFrame(Type="logical", 
                                                                DefaultValue="TRUE",
                                                                Description="Perform batch correction."),
                                        integrationMethod=ezFrame(Type="character", 
                                                                  DefaultValue="Classic", 
                                                                  Description="Choose integration method in Seurat (Classic or RPCA)"),
                                        SCT.regress=ezFrame(Type="character", 
                                                            DefaultValue="none", 
                                                            Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."),
                                        DE.method=ezFrame(Type="charVector", 
                                                          DefaultValue="wilcox", 
                                                          Description="Method to be used when calculating gene cluster markers and differentially expressed genes between conditions. Use LR to take into account the Batch and/or CellCycle"),
                                        DE.regress=ezFrame(Type="charVector", 
                                                           DefaultValue="Batch", 
                                                           Description="Variables to regress out if the test LR is chosen"))
                }
              )
  )

ezMethodScSeuratCombine = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(rlist)
  library(HDF5Array)
  library(scanalysis)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(AUCell)
  library(enrichR)
  library(decoupleR)
  library(BiocParallel)
  
  BPPARAM <- MulticoreParam(workers = param$cores)
  register(BPPARAM)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  #the individual sce objects can be in hdf5 format (for new reports) or in rds format (for old reports)
  filePath <- file.path("/srv/gstore/projects", input$getColumn("SC Cluster Report"), 'scData.rds')
  filePath_course <- file.path("/srv/GT/analysis/course_sushi/public/projects", input$getColumn("SC Cluster Report"), 'scData.rds')
  
  if(!file.exists(filePath[1])) 
    filePath <- filePath_course
  names(filePath) <- input$getNames()
  
  # Load the data and prepare metadata for integration
  scDataList <- lapply(names(filePath), function(sm) {
    scData <- readRDS(filePath[sm])
    scData$Sample <- sm
    # If we have new information in the Condition column, add it to the dataset
    if (all(scData$Condition == "NA" | scData$Condition == "")) {
      if (any(startsWith(colnames(input$meta), "Condition"))) {
        scData$Condition <- unname(input$getColumn("Condition")[sm])
      } else {
        scData$Condition <- scData$Sample
      }
    }
    # Also add the other factors in the input dataset to the objects
    if (!is.null(param$harmonyFactors)) {
      harmonyFactors <- str_split(param$harmonyFactors, ",", simplify=TRUE)[1,]
      metaFactorNames <- paste0("har_", harmonyFactors) %>% str_replace(., " ", ".")
      names(metaFactorNames) <- harmonyFactors
      for (hf in harmonyFactors) {
        scData[[metaFactorNames[hf]]] <- unname(input$getColumn(hf)[sm])
      }
    }
    # Rename the cells and add original sample-level clusters back in
    scData <- RenameCells(scData, new.names = paste0(scData$Sample, "-", colnames(scData)))
    scData$sample_seurat_clusters <- paste0(scData$Sample, "-", sprintf("%02d", scData$seurat_clusters))
    return(scData)
  })
  
  pvalue_allMarkers <- 0.05
  nrSamples <- length(scDataList)
  
  if(ezIsSpecified(param$chosenClusters)){
    for(eachSample in names(param$chosenClusters)){
      chosenCells <- names(Idents(scDataList[[eachSample]]))[Idents(scDataList[[eachSample]]) %in% param$chosenClusters[[eachSample]]]
      scDataList[[eachSample]] <- scDataList[[eachSample]][, chosenCells]
    }
  }
  
  scData_noCorrected <- cellClustNoCorrection(scDataList, param)
  scData = scData_noCorrected
  if (param$batchCorrection) {
    scData_corrected = cellClustWithCorrection(scDataList, param)
    #in order to compute the markers we switch again to the original assay
    DefaultAssay(scData_corrected) <- "SCT"
    scData <- scData_corrected
  }
  scData@reductions$tsne_noCorrected <- Reductions(scData_noCorrected, "tsne")
  scData@reductions$umap_noCorrected <- Reductions(scData_noCorrected, "umap")
  scData@meta.data$ident_noCorrected <- Idents(scData_noCorrected)
  scData <- PrepSCTFindMarkers(scData)
  
  # positive cluster markers
  markers <- FindAllMarkers(object=scData, test.use = param$DE.method, only.pos=TRUE)
  ## Significant markers
  markers <- markers[ ,c("gene","cluster","pct.1", "pct.2", "avg_log2FC","p_val_adj")]
  markers$cluster <- as.factor(markers$cluster)
  markers$diff_pct = abs(markers$pct.1-markers$pct.2)
  markers <- markers[order(markers$diff_pct, decreasing = TRUE),]
  writexl::write_xlsx(markers, path="posMarkers.xlsx")
  
  # we do cell type identification using AUCell and SingleR
  cells_AUC <- NULL
  singler.results <- NULL
  #cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if(species == "Human" | species == "Mouse") {
      genesPerCluster <- split(markers$gene, markers$cluster)
      enrichRout <- querySignificantClusterAnnotationEnrichR(genesPerCluster, unlist(strsplit(param$enrichrDatabase, ',')))
      cells.AUC <- cellsLabelsWithAUC(GetAssayData(scData, "counts"), species, param$tissue, BPPARAM = BPPARAM)
      singler.results <- cellsLabelsWithSingleR(GetAssayData(scData, "data"), Idents(scData), species, BPPARAM = BPPARAM)
      for (r in names(singler.results)) {
          scData[[paste0(r,"_single")]] <- singler.results[[r]]$single.fine$labels
          scData[[paste0(r,"_cluster")]] <- singler.results[[r]]$cluster.fine$labels[match(Idents(scData), rownames(singler.results[[r]]$cluster.fine))]
      }
    saveRDS(cells.AUC, file="cells.AUC.rds")
    saveRDS(singler.results, file="singler.results.rds")
  } else {
      cells.AUC <- NULL
      singler.results <- NULL
      enrichRout <- NULL
  }
  
  ## SCpubr advanced plots
  pathwayActivity <- computePathwayActivityAnalysis(cells = scData, species = species)
  TFActivity <- computeTFActivityAnalysis(cells = scData, species = species)
  
  geneMeans <- geneMeansCluster(scData)
  
  saveRDS(param, file="param.rds")
  saveRDS(output, file="output.rds")
  saveRDS(scData, file = "scData.rds")
  saveRDS(enrichRout, file = "enrichRout.rds")
  saveRDS(pathwayActivity, file = "pathwayActivity.rds")
  saveRDS(TFActivity, file = "TFActivity.rds")
  
  # Save some results in external files
  dataFiles <- saveExternalFiles(list(gene_means = as_tibble(as.data.frame(geneMeans), rownames = "gene_name")))
  reportTitle <- 'SCReport - MultipleSamples based on Seurat'
  makeRmdReport(dataFiles=dataFiles, rmdFile = "ScSeuratCombine.Rmd", reportTitle = reportTitle) 
  return("Success")
  
}