###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCMultipleSamplesAndGroups <-
  setRefClass("EzAppSCMultipleSamplesAndGroups",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCMultipleSamplesAndGroups
                  name <<- "EzAppSCMultipleSamplesAndGroups"
                  appDefaults <<- rbind(npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=30, 
                                                    Description="The maximal dimensions to use for reduction"),
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        batchCorrection=ezFrame(Type="logical", 
                                                                DefaultValue="TRUE",
                                                                Description="Which batch correction method to use? None or CCA"),
                                        chosenClusters=ezFrame(Type="charList",
                                                               DefaultValue="",
                                                               Description="The clusters to choose from each sample.In the format of sample1=cluster1,cluster2;sample2=cluster1,cluster2."),
                                        all2allMarkers=ezFrame(Type="logical",
                                                               DefaultValue=FALSE, 
                                                               Description="Run all against all cluster comparisons?"),
                                        markersToShow=ezFrame(Type="numeric", 
                                                              DefaultValue=10, 
                                                              Description="The markers to show in the heatmap of cluster marker genes"),
                                        maxSamplesSupported=ezFrame(Type="numeric", 
                                                              DefaultValue=5, 
                                                              Description="Maximum number of samples to compare"))
                }
              )
  )

ezMethodSCMultipleSamplesAndGroups = function(input=NA, output=NA, param=NA, htmlFile="00index.html") {
  library(Seurat)
  library(rlist)
  library(tibble)
  library(dplyr)
  library(readr)
  
  ## subset the selected sample names
  samples <- param$samples
  input <- input$subset(samples)
  
  if(input$getLength() > param$maxSamplesSupported){
    stop(paste("It only works for", param$maxSamplesSupported, "at most"))
  }
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  param$name <- paste(param$name, paste(input$getNames(), collapse=", "), sep=": ")
  
  sceURLs <- input$getColumn("Static Report")
  
  saveRDS(param, file = "param.rds")
  
  sceList <- lapply(file.path("/srv/gstore/projects", sub("https://fgcz-(gstore|sushi).uzh.ch/projects", "",dirname(sceURLs)), "sce.rds"),readRDS)
  names(sceList) <- names(sceURLs)
  sceList = lapply(sceList, update_seuratObjectVersion)
  
  pvalue_allMarkers <- 0.05
  pvalue_all2allMarkers <- 0.01
  nrSamples <- length(sceList)
  
  if(ezIsSpecified(param$chosenClusters)){
    for(eachSample in names(param$chosenClusters)){
      chosenCells <- names(Idents(metadata(sceList[[eachSample]])$scData))[Idents(metadata(sceList[[eachSample]])$scData) %in% param$chosenClusters[[eachSample]]]
      sceList[[eachSample]] <- sceList[[eachSample]][, chosenCells]
      metadata(sceList[[eachSample]])$scData <-
        SubsetData(metadata(sceList[[eachSample]])$scData,
                   ident.use=param$chosenClusters[[eachSample]])
    }
  }
  
  scData_noCorrected = cellClustNoCorrection(sceList, param)
  scData = scData_noCorrected
  if (param$batchCorrection) {
     scData_corrected = cellClustWithCorrection(sceList, param)
     #in order to compute the markers we switch again to the original assay
     DefaultAssay(scData_corrected) <- "RNA"
     scData_corrected <- ScaleData(scData_corrected, verbose = FALSE)
     scData = scData_corrected
  }
     
  #positive cluster markers
  scData = posClusterMarkers(scData, pvalue_allMarkers)
  
  #perform all pairwise comparisons to obtain markers
  if(doEnrichr(param) && param$all2allMarkers) 
    scData = all2all(scData, pvalue_all2allMarkers, param)
  
  #conserved cluster markers
  scData_objects = conservedMarkers(scData)
  scData = scData_objects[[1]]
  scData_Clustersubset = scData_objects[[2]]
  
  #differentially expressed genes
  scData = diffExpressedGenes(scData, scData_Clustersubset)
  
  scData = saveExternalFiles(scData)
  
  scData_list = list()
  if(param$batchCorrection) {
    scData_corrected = scData
    scData_list = list.append(scData_list, scData_corrected = scData_corrected)
  } else {
      scData_noCorrected = scData
  }
  scData_list = list.append(scData_list, scData_noCorrected = scData_noCorrected) #this list will contain the non-corrected object and also the corrected object when calculated
  
  saveRDS(scData_list, "scData_ObjectList.rds")
  saveRDS(as.SingleCellExperiment(scData), "sce_iSEE.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCMultipleSamplesAndGroups.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCMultipleSamplesAndGroups.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  # rmarkdown::render(input="/home/daymegr/workspaceR/dayme-scripts/sushi_scripts_mod/SCMultipleSamplesAndGroups.Rmd", envir = new.env(),
  #                   output_dir=".", output_file=htmlFile, clean = TRUE, quiet=TRUE)
  rm(sceList)
  rm(scData)
  
  return("Success")
  
}

cellClustNoCorrection = function(sceList, param) {
  scData = Reduce(merge, lapply(sceList, function(x){metadata(x)$scData}))
  scData@project.name <- param$name
  scData <- NormalizeData(object = scData)
  scData <- FindVariableFeatures(object = scData)
  scData <- seuratStandardWorkflow(scData, param)
  return(scData)
}

cellClustWithCorrection = function (sceList, param) {
  seurat_objects = lapply(sceList, function(se) {metadata(se)$scData})
  scData = seuratIntegration(seurat_objects, param)
  # switch to integrated assay for downstream analyses
  DefaultAssay(scData) <- "integrated"
  # Run the standard workflow for visualization and clustering
  scData = seuratStandardWorkflow(scData, param)
  return(scData)
}

posClusterMarkers = function(scData, pvalue_allMarkers) {
  markers <- FindAllMarkers(object=scData, only.pos=TRUE, return.thresh = pvalue_allMarkers)
  ## Significant markers
  cm <- markers[ ,c("gene","cluster","avg_logFC","p_val_adj")]
  rownames(cm) <- NULL
  scData@misc$posMarkers = cm
  return(scData)
}

all2all = function(scData, pvalue_all2allMarkers, param) {
  clusterCombs <- combn(levels(Idents(scData)), m=2)
  all2allMarkers <- mcmapply(FindMarkers, as.integer(clusterCombs[1, ]), as.integer(clusterCombs[2, ]),
                             MoreArgs = list(object=scData,only.pos=FALSE),
                             mc.preschedule=FALSE,
                             mc.cores=min(4L, param$cores),
                             SIMPLIFY=FALSE)
  all2allMarkers <- lapply(all2allMarkers, function(x){
    x[x$p_val <= pvalue_all2allMarkers, ]
  })
  names(all2allMarkers) <- apply(clusterCombs, 2, paste, collapse="vs")
  scData@misc$all2allMarkers = all2allMarkers
  return(scData)
}


conservedMarkers = function(scData) {
#When finding the conserved markers and differentially expressed genes across conditions I will discard the clusters that were too small in at least one group
clusters_freq = data.frame(table(scData@meta.data[,c("orig.ident","seurat_clusters")]))
small_clusters = unique(as.character(clusters_freq[clusters_freq[,"Freq"] < 10, "seurat_clusters"]))
scData_Clustersubset = subset(scData, idents = small_clusters, invert = TRUE)

markers <- list()
for(eachCluster in levels(Idents(scData_Clustersubset))){
  markersEach <- try(FindConservedMarkers(scData_Clustersubset, ident.1=eachCluster, grouping.var="orig.ident", print.bar=FALSE, only.pos=TRUE,silent=TRUE))
  ## to skip some groups with few cells
  if(class(markersEach) != "try-error" && nrow(markersEach) > 0){
    markers[[eachCluster]] <- as_tibble(markersEach, rownames="gene")
  }
}
## some of the cluster have no significant conserved markers
markers <- markers[sapply(markers, nrow) != 0L] 
markers <- bind_rows(markers, .id="cluster")

scData@misc$conservedMarkers = markers
return(c(scData, scData_Clustersubset))
}

diffExpressedGenes = function(scData, scData_Clustersubset) {
  seurat_clusters = Idents(scData_Clustersubset)
  scData_Clustersubset@meta.data$cluster.condition <- paste0(seurat_clusters, "_", scData_Clustersubset@meta.data$orig.ident)
  Idents(scData_Clustersubset) = "cluster.condition"
  conditions <- unique(scData_Clustersubset@meta.data$orig.ident)

diffGenesFns <- c()
conditionsComb <- combn(conditions, m=2)
for(i in 1:ncol(conditionsComb)){
  diffGenes <- list()
  for(eachCluster in gtools::mixedsort(levels(seurat_clusters))){
    markersEach <- try(FindMarkers(scData_Clustersubset, ident.1=paste0(eachCluster, "_",
                                                          conditionsComb[2,i]),
                                   ident.2=paste0(eachCluster, "_", 
                                                  conditionsComb[1,i]),
                                   print.bar=FALSE), silent=TRUE)
    ## to skip some groups with few cells
    if(class(markersEach) != "try-error"){
      diffGenes[[eachCluster]] <- as_tibble(markersEach, rownames="gene")
    }
  }
  diffGenes <- bind_rows(diffGenes, .id="cluster")
}
diffGenes = diffGenes[diffGenes$p_val_adj < 0.05, ]
scData@misc$diffGenes = diffGenes
return(scData)
}


saveExternalFiles = function(scData) {
  library(readr)
  
  propCells_table = cellsProportion(scData)
  cellsPropPerClusterAndSampleFn = "cells_proportions.txt"
  scData@misc$cellsPropPerClusterAndSampleFn <- cellsPropPerClusterAndSampleFn
  ezWrite.table(propCells_table, cellsPropPerClusterAndSampleFn)

tr_cnts <- expm1(GetAssayData(scData))
geneMeans <- rowsum(t(as.matrix(tr_cnts)), group=Idents(scData))
geneMeans <- sweep(geneMeans, 1, STATS=table(Idents(scData))[rownames(geneMeans)], FUN="/")
geneMeans <- log1p(t(geneMeans))
colnames(geneMeans) <- paste("cluster", colnames(geneMeans), sep="_")
geneMeanPerClusterFn = "gene_means_per_cluster.txt"
scData@misc$geneMeanPerClusterFn <- geneMeanPerClusterFn
ezWrite.table(geneMeans, geneMeanPerClusterFn)

geneMeans <- Matrix::rowMeans(tr_cnts)
geneMeans <- log1p(geneMeans)
geneMeansFn = "gene_means.txt"
scData@misc$geneMeansFn <- geneMeansFn
ezWrite.table(geneMeans, geneMeansFn)

tSNE_data <- as_tibble(scData@reductions$tsne@cell.embeddings,
                       rownames="cells")
tSNE_data <- dplyr::rename(tSNE_data, X=`tSNE_1`, Y=`tSNE_2`)
tSNE_data$cluster <- Idents(scData)
tSNEFn = "tSNE_data.tsv"
scData@misc$tSNEFn <- tSNEFn
write_tsv(tSNE_data, path=tSNEFn)

posMarkersFn <- "pos_markers.tsv"
write_tsv(as_tibble(scData@misc$posMarkers), path=posMarkersFn)

return(scData)
}