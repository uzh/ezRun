###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCFeatBarcoding <-
  setRefClass("EzAppSCFeatBarcoding",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCFeatBarcoding
                  name <<- "EzAppSCFeatBarcoding"
                  appDefaults <<- rbind(npcs=ezFrame(Type="numeric", 
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
                                        nreads = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = Inf,
                                          Description = "Low quality cells have less than \"nreads\" reads. Only when applying fixed thresholds."
                                        ),
                                        ngenes = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = Inf,
                                          Description = "Low quality cells have less than \"ngenes\" genes. Only when applying fixed thresholds."
                                        ),
                                        perc_mito = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = Inf,
                                          Description = "Low quality cells have more than \"perc_mito\" percent of mitochondrial genes. Only when applying fixed thresholds."
                                        ),
                                        perc_ribo = ezFrame(
                                          Type = "numeric",
                                          DefaultValue = Inf,
                                          Description = "Low quality cells have more than \"perc_ribo\" percent of ribosomal genes. Only when applying fixed thresholds."
                                        ),
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

ezMethodSCFeatBarcoding <- function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  library(scanalysis)
  library(HDF5Array)
  library(AUCell)
  library(GSEABase)
  library(SingleR)
  library(Seurat)
  library(SingleCellExperiment)
  library(tidyverse)
  library(scanalysis)
  require(scDblFinder)
  library(BiocParallel)

  if (param$cores > 1){
        BPPARAM <- MulticoreParam(workers = param$cores)
    } else {
        ## scDblFinder fails with many cells and MulticoreParam
        BPPARAM <- SerialParam() 
    }
  register(BPPARAM)
    

    
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)

  scData <- load10xSC_seurat(input, param)
  hto_assay <- GetAssay(scData, "HTO") #Save assay to add it back to the scData object after filtering genes. Filtering genes remove non active assays.
  
  # Cells and genes filtering
  scData_list <- filterCellsAndGenes(scData, param) # return sce objects filtered and unfiltered to show the QC metrics later in the rmd
  scData <- scData_list$scData
  #Subset hto assay too and add it back to the seurat object
  hto_assay <- subset(hto_assay, cells=colnames(scData))
  scData[["HTO"]] <- hto_assay
  scData.unfiltered <- scData_list$scData.unfiltered
  rm(scData_list)
  
  #Normalize antibody counts using the CLR method
  scData <- NormalizeData(scData, assay = "HTO", normalization.method = "CLR")
  #Demultiplex cells based on antibody enrichments
  scData <- HTODemux(scData, assay = "HTO", positive.quantile = 0.99)
  
  #Clustering on RNA expression using only the Singlet cells
  Idents(scData) <- "HTO_classification.global"
  scData.singlet <- subset(scData, idents = "Singlet")
  
  # calculate cellcycle for the singlets scData object
  scData.singlet <- addCellCycleToSeurat(scData.singlet, param$refBuild, BPPARAM)
  
  #Clustering on protein levels and on RNA levels
  scData.singlet <- seuratClusteringHTO(scData.singlet)
  scData.singlet <- seuratClusteringV3(scData.singlet, param)
  
  #positive cluster markers
  pvalue_allMarkers <- 0.05
  DefaultAssay(scData.singlet) <- "SCT"
  posMarkers <- posClusterMarkers(scData.singlet, pvalue_allMarkers, param)
 
  
  #cell types annotation is only supported for Human and Mouse at the moment
  species <- getSpecies(param$refBuild)
  if(species == "Human" | species == "Mouse") {
     cells_AUC <- cellsLabelsWithAUC(GetAssayData(scData.singlet, "counts"), species, param$tissue, BPPARAM=BPPARAM)
     singler.results <- cellsLabelsWithSingleR(GetAssayData(scData.singlet, "counts"), Idents(scData.singlet), species)
     for (r in names(singler.results)) {
         scData.singlet[[paste0(r,"_single")]] <- singler.results[[r]]$single.fine$labels
         scData.singlet[[paste0(r,"_cluster")]] <- singler.results[[r]]$cluster.fine$labels[match(Idents(scData.singlet), rownames(singler.results[[r]]$cluster.fine))]
     }
     saveRDS(cells_AUC, file="cells.AUC.rds")
     saveRDS(singler.results, file="singler.results.rds")
  } else {
      cells_AUC <- NULL
      singler.results <- NULL
  }
  
  #Convert scData to Single Cell experiment Object
  sce.unfiltered <- scData.unfiltered %>% seurat_to_sce()
  sce.singlets = scData.singlet %>% seurat_to_sce(default_assay = "SCT") #SCT as default assay for visualization
  # Doublets prediction (no removal)
  sce.singlets <- scDblFinder(sce.singlets, clusters=TRUE, BPPARAM = BPPARAM)
  metadata(sce.singlets)$PCA_stdev <- Reductions(scData.singlet, "pca")@stdev   
  metadata(sce.singlets)$cells_AUC <- cells_AUC
  metadata(sce.singlets)$singler.results <- singler.results
  metadata(sce.singlets)$output <- output
  metadata(sce.singlets)$param <- param
  metadata(sce.singlets)$param$name <- paste(metadata(sce.singlets)$param$name,paste(input$getNames(), collapse=", "),sep=": ")
  
  
  #Save some results in external files 
  saveRDS(scData, "scData.rds") #this sce object contains all Doublets, Negative and Singlet cells. Only used for demultiplexing
  saveHDF5SummarizedExperiment(sce.singlets, dir="sce_h5") #singlets used for downstream analyses
  saveHDF5SummarizedExperiment(sce.unfiltered, dir="sce.unfiltered_h5") #all cells with QC metrics to show on the rmd
  geneMeans <- geneMeansCluster(sce.singlets)
  dataFiles = saveExternalFiles(list(pos_markers=posMarkers, gene_means=as_tibble(as.data.frame(geneMeans), rownames="gene_name")))
  
  makeRmdReport(dataFiles=dataFiles, rmdFile = "SCFeatBarcoding.Rmd", reportTitle = metadata(sce.singlets)$param$name)
  #remove no longer used objects
  rm(scData, sce.singlets, sce.unfiltered)
  gc()
  return("Success")
}


