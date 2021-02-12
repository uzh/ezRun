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
                                        cellsFraction=ezFrame(Type="numeric", 
                                                                DefaultValue=0.05, 
                                                                Description="A gene will be kept if it is expressed in at least this percentage of cells"),
                                        nUMIs=ezFrame(Type="numeric", 
                                                      DefaultValue=1, 
                                                      Description='A gene will be kept if it has at least nUMIs in the fraction of cells specified before'),
                                        nmad=ezFrame(Type="numeric", 
                                                     DefaultValue=3, 
                                                     Description="Median absolute deviation (MAD) from the median value of each metric across all cells"),
                                        species=ezFrame(Type="character", DefaultValue="Human", Description="Organism")
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
  library(tibble)
  library(scanalysis)
  require(scDblFinder)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)

  data10X <- Read10X(input$getFullPaths("CountMatrix"))
  
  #Create a SingleCellExperiment using the Gene expression counts
  data10X_GeneExp <- data10X$`Gene Expression`
  cellCycleFn <- list.files(path=input$getFullPaths("CountMatrix"),
                            pattern="CellCyclePhase\\.txt$", recursive=TRUE,
                            full.names=TRUE)
  cellCycle <- ezRead.table(cellCycleFn)
  colData <- DataFrame(CellCycle=cellCycle$Phase,
                       CellCycleG1=cellCycle$G1,
                       CellCycleS=cellCycle$S,
                       CellCycleG2M=cellCycle$G2M,
                       Batch=input$getNames(),
                       Condition=try(input$getColumn("Condition"), silent = TRUE))
  
 
  sce <- SingleCellExperiment(assays=list(counts=data10X_GeneExp), 
                                      metadata=list(param=param),
                                      rowData = DataFrame(gene_name=rownames(data10X_GeneExp)),
                                      colData=colData)
  
  #Cells and genes filtering
  sce_list <- filterCellsAndGenes(sce, param)  #return sce objects filtered and unfiltered to show the QC metrics later in the rmd 
  sce <- sce_list$sce
  sce.unfiltered <- sce_list$sce.unfiltered
  
  # the Seurat object is built from the Gene expression filtered sce object
  scData <- CreateSeuratObject(counts=counts(sce), project=param$name, meta.data=data.frame(colData(sce)))
  
  #Create an assay for the hashtags and add it to the seurat object
  scData[["HTO"]] <- CreateAssayObject(data10X[["Antibody Capture"]][, colnames(x = scData)])
  #Normalize antibody counts using the CLR method
  scData <- NormalizeData(scData, assay = "HTO", normalization.method = "CLR")
  #Demultiplex cells based on antibody enrichments
  scData <- HTODemux(scData, assay = "HTO", positive.quantile = 0.99)
  
  #Clustering on RNA expression using only the Singlet cells
  Idents(scData) <- "HTO_classification.global"
  scData.singlet <- subset(scData, idents = "Singlet")
  #Clustering on protein levels and on RNA levels
  scData.singlet <- seuratClusteringHTO(scData.singlet)
  scData.singlet <- seuratClusteringV3(scData.singlet, param)
  
  #positive cluster markers
  pvalue_allMarkers <- 0.05
  DefaultAssay(scData.singlet) <- "SCT"
  posMarkers <- posClusterMarkers(scData.singlet, pvalue_allMarkers, param)
 
  cells_AUC <- NULL
  singler.results <- NULL
  #cell types annotation is only supported for Human and Mouse at the moment
  if(param$species == "Human" | param$species == "Mouse") {
     cells_AUC <- cellsLabelsWithAUC(scData.singlet, param)
     singler.results <- cellsLabelsWithSingleR(GetAssayData(scData.singlet, "counts"), Idents(scData.singlet), param)
  }
  
  #Convert scData to Single Cell experiment Object
  #sce = scData %>% seurat_to_sce() #this sce object contains all Doublets, Negative and Singlet cells. No downstream analyses were performed. Only used for demultiplexing
  sce.singlets = scData.singlet %>% seurat_to_sce(default_assay = "SCT") #SCT as default assay for visualization
  sce.singlets <- scDblFinder(sce.singlets)
  metadata(sce.singlets)$PCA_stdev <- Reductions(scData.singlet, "pca")@stdev   
  metadata(sce.singlets)$cells_AUC <- cells_AUC
  metadata(sce.singlets)$singler.results <- singler.results
  metadata(sce.singlets)$output <- output
  metadata(sce.singlets)$param <- param
  metadata(sce.singlets)$param$name <- paste(metadata(sce.singlets)$param$name,paste(input$getNames(), collapse=", "),sep=": ")
  
  
  #Save some results in external files 
  saveRDS(scData, "scData.rds")
  saveHDF5SummarizedExperiment(sce.singlets, dir="sce_h5")
  saveExternalFiles(sce.singlets, list(pos_markers=posMarkers))
 # rowData(sce) = rowData(sce)[, c("gene_id", "biotypes", "description")]
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCFeatBarcoding.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  while (dev.cur()>1) dev.off()
  rmarkdown::render(input="SCFeatBarcoding.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)

  return("Success")
}


