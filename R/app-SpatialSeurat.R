###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSpatialSeurat <-
  setRefClass("EzAppSpatialSeurat",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSpatialSeurat
                  name <<- "EzAppSpatialSeurat"
                  appDefaults <<- rbind(npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=20,
                                                    Description="The maximal dimensions to use for reduction"),
                                        pcGenes=ezFrame(Type="charVector", 
                                                        DefaultValue="", 
                                                        Description="The genes used in supvervised clustering"),
                                        SCT.regress.CellCycle=ezFrame(
                                          Type = "logical", 
                                          DefaultValue = FALSE,
                                          Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
                                        ),
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
                                        spotClean = ezFrame(
                                            Type = "logical",
                                            DefaultValue = FALSE,
                                            Description = "Run spotClean method"
                                        )
                                        )
                }
              )
  )

ezMethodSpatialSeurat <- function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  library(Seurat)
  library(scater)
  
  cmDir <- input$getFullPaths("CountMatrix")
  featInfo <- ezRead.table(paste0(cmDir, "/features.tsv.gz"), header = FALSE, row.names = NULL)
  colnames(featInfo) <- c("gene_id", "gene_name", "type")
  featInfo$isMito = grepl( "(?i)^MT-", featInfo$gene_name)
  featInfo$isRiboprot = grepl(  "(?i)^RPS|^RPL", featInfo$gene_name)
  geneAnnoFile <- sub("byTranscript", "byGene", param$ezRef@refAnnotationFile)
  if (file.exists(geneAnnoFile)){
      geneAnno <- ezRead.table(geneAnnoFile)
      if (any(geneAnno$type == "rRNA")){
          featInfo$isRibosomal <- geneAnno[featInfo$gene_id, "type"] == "rRNA"
          if(any(is.na(featInfo[, "isRibosomal"]))){
              featInfo[, "isRibosomal"][which(is.na(featInfo[, "isRibosomal"]))] <- FALSE
          }
      }
  }
  rownames(featInfo) <- gsub("_", "-", uniquifyFeatureNames(ID=featInfo$gene_id, names=featInfo$gene_name)) 
  
  if(param$spotClean){
      scData <- load10xSpatialDataAndRunSpotClean(input, param)
      param$imageEnlargementFactor <- 1 
  } else {
      res <- load10xSpatialData(input, param)
      scData <- res[[1]]
      param <- res[[2]]
      remove(res)
  }
  scData$Condition <- input$getColumn("Condition")
  scData@meta.data$Sample <- input$getNames()
  scData[["Spatial"]] <- AddMetaData(object = scData[["Spatial"]], metadata = featInfo[rownames(scData), ])
  
  scData_list <- filterCellsAndGenes(scData, param) # return sce objects filtered and unfiltered to show the QC metrics later in the rmd
  scData <- scData_list$scData
  scData.unfiltered <- scData_list$scData.unfiltered
  rm(scData_list)
  
  scData <- addCellCycleToSeurat(scData, param$refBuild, BPPARAM, assay = 'Spatial')
  scData.unfiltered <- addCellCycleToSeurat(scData.unfiltered, param$refBuild, BPPARAM, assay = 'Spatial')
  scData <- seuratClusteringV3(scData, param, assay="Spatial")
  
  pvalue_allMarkers <- 0.05
  
  #positive cluster markers
  clusterMarkers <- posClusterMarkers(scData, pvalue_allMarkers, param)
  clusterMarkers[['isSpatialMarker']] = FALSE 
  #spatially variable genes
  spatialMarkersList <- list()
  res <- spatialMarkers(scData, selection.method = 'markvariogram')
  spatialMarkersList[['markvariogram']] <- data.frame(GeneSymbol = rownames(res), res, Method = 'Markvariogram')
  res <- spatialMarkers(scData, selection.method = 'moransi')
  spatialMarkersList[['moransi']] <- data.frame(GeneSymbol = rownames(res), res, Method = 'MoransI')
  spatialMarkers <- rbind(spatialMarkersList[['markvariogram']][,c('GeneSymbol', 'Rank','Method')], spatialMarkersList[['moransi']][,c('GeneSymbol', 'Rank','Method')])
  spatialMarkers <- spatialMarkers %>% spread(Method, Rank)
  spatialMarkers[['MeanRank']] <- apply(spatialMarkers[,c('Markvariogram','MoransI')],1,mean)
  spatialMarkers <- spatialMarkers[order(spatialMarkers$MeanRank),]
  
  spatialPosMarkers <- intersect(clusterMarkers$gene, spatialMarkers$GeneSymbol)
  clusterMarkers[which(clusterMarkers$gene %in% spatialPosMarkers), 'isSpatialMarker'] = TRUE
  
  #Save some results in external files
  library(scanalysis)
  scData_diet = DietSeurat(scData, dimreducs = c("pca", "tsne", "umap"))
  sce <- scData_diet %>% seurat_to_sce(default_assay = "SCT")
  geneMeansPerCluster <- geneMeansCluster(sce)
  geneMeans <- apply(logcounts(sce), 1, mean)
  geneMeans <- data.frame(logCount = geneMeans, row.names = names(geneMeans))
  dataFiles <- saveExternalFiles(list(cluster_markers=clusterMarkers, spatial_markers=data.frame(spatialMarkers), gene_means = geneMeans, gene_means_per_cluster = as_tibble(as.data.frame(geneMeansPerCluster), rownames = "gene_name")))
  
  saveRDS(scData, "scData.rds")
  allCellsMeta <- scData.unfiltered@meta.data
  allCellsMeta$Sample <- allCellsMeta$Batch
  allCellsMeta$useCell <- !allCellsMeta$discard
  saveRDS(allCellsMeta, 'allCellsMeta.rds')
  saveRDS(scData.unfiltered, "scData.unfiltered.rds")
  saveRDS(param, "param.rds")
  
  
  makeRmdReport(dataFiles=dataFiles, rmdFile = "SpatialSeurat.Rmd", reportTitle = param$name) 
  return("Success")
}


