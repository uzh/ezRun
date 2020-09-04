###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCMergeLargeDatasets <-
  setRefClass("EzAppSCMergeLargeDatasets",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCMergeLargeDatasets
                  name <<- "EzAppSCMergeLargeDatasets"
                  appDefaults <<- rbind(k=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="A numeric value specifying the number of nearest neighbors to consider during graph construction"),
                                        block=ezFrame(Type="character", DefaultValue="Batch", Description="Blocking factor for each cell to consider during DE test"),
                                        species=ezFrame(Type="character", DefaultValue="Human", Description="Organism"))
                }
              )
  )

ezMethodSCMergeLargeDatasets = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(HDF5Array)
  library(scran)
  library(harmony)
  library(scater)
  library(readr)
  library(dplyr)
  
  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  sceURLs <- input$getColumn("Static Report")
  filePath <- file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(sceURLs)), "sce_h5")
  if(file.exists(filePath)) {
    sceList <- lapply(filePath,loadHDF5SummarizedExperiment)
    names(sceList) <- names(sceURLs)
  }
  
  #Only use common genes among all samples
  common_genes = Reduce(intersect, lapply(sceList, rownames)) 
  sceList = lapply(sceList, function(sce) {rowData(sce) = NULL; sce[common_genes,]})
  sce =  Reduce(SingleCellExperiment::cbind, sceList) 
  
  set.seed(1000)
  
  # Normalization
  clusters <- quickCluster(sce)
  sce = computeSumFactors(sce, cluster = clusters)
  sce = logNormCounts(sce)
  
  # Variance modelling
  dec <- modelGeneVar(sce, block=sce$Batch)
  
  # Extract variable genes for downtream analysis
  top.hvgs <- getTopHVGs(dec, n=2000)
  
  # PCA
  sce <- denoisePCA(sce, dec, subset.row=top.hvgs)
  
  # Integration with Harmony using the PCA embeddings
  harmony_PCA <- HarmonyMatrix(
    data_mat  = reducedDim(sce, "PCA"),
    meta_data = colData(sce),
    vars_use  = "Batch", 
    do_pca = FALSE, 
  )
  reducedDim(sce, "harmony_PCA") = harmony_PCA 
   
 # Run UMAP using uncorrected PCAs and the corrected PCAs obtained from Harmony
  sce <- runUMAP(sce, dimred = "PCA", name="uncorrected_UMAP")
  sce <- runUMAP(sce, dimred = "harmony_PCA", name="harmony_UMAP")              
  
 # Clustering 
  g <- buildSNNGraph(sce, use.dimred="harmony_PCA")
  cluster <- igraph::cluster_walktrap(g)$membership
  sce$cluster <- factor(cluster)
                 
  #positive cluster markers
  markers_any <- findMarkers(sce, sce$cluster, pval.type="any", direction="up", lfc=1, block=sce$Batch)
  markersList <- list()
  for(cluster in names(markers_any)){
    markersPerCluster <- data.frame(gene_name = rownames(markers_any[[cluster]]), markers_any[[cluster]], row.names = NULL)
    markersPerCluster <- markersPerCluster[markersPerCluster$FDR<0.05, ]
    markersList[[cluster]] <- markersPerCluster
  }
   markers_any <- bind_rows(markersList, .id="Cluster")
  
  
  #we do cell type identification only with SingleR since it implemments  block-processing
  singler.results <- NULL
  if(param$species == "Human" | param$species == "Mouse") 
    singler.results <- cellsLabelsWithSingleR(counts(sce), sce$cluster, param)
  
   #Convert scData to Single Cell experiment Object
   # Use the SCT logcounts for visualization instead of the RNA logcounts.
   # TODO: save all the assays (RNA, SCT and integrated) in the sce object using the package keshavmot2/scanalysis. The function from Seurat doesn't save everything.
   DefaultAssay(scData) <- "SCT" 
   sce <- as.SingleCellExperiment(scData)
   metadata(sce)$singler.results <- singler.results
   metadata(sce)$output <- output
   metadata(sce)$param <- param
   metadata(sce)$param$name <- paste(param$name, paste(input$getNames(), collapse=", "), sep=": ")
   
   #Save some results in external files 
  saveExternalFiles(sce, list(pos_markers=markers_any))
  saveHDF5SummarizedExperiment(sce, dir="sce_h5")
  
  # Copy the style files and templates
   styleFiles <- file.path(system.file("templates", package="ezRun"), c("fgcz.css", "SCMergeLargeDatasets.Rmd", "fgcz_header.html", "banner.png"))
 # styleFiles <- paste0("/home/daymegr/workspaceR/dayme-scripts/sushi_scripts_mod/", c("fgcz.css", "SCMergeLargeDatasets.Rmd", "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCMergeLargeDatasets.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, clean = TRUE, quiet=TRUE)
  rm(sceList)
  return("Success")
  
}

# createMarkersList <- function(markers) {
#   markersList <- list()
#  for(cluster in names(markers)){
#    markersPerCluster <- data.frame(gene_name = rownames(markers[[cluster]]), markers[[cluster]], row.names = NULL)
#    markersPerCluster <- markersPerCluster[markersPerCluster$FDR<0.05, ]
#    markersList[[cluster]] <- markersPerCluster
#  }
#  markersList <- bind_rows(markersList, .id="Cluster")
#  return(markersList)
# }
