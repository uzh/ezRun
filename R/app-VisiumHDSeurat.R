###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppVisiumHDSeurat <-
  setRefClass("EzAppVisiumHDSeurat",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodVisiumHDSeurat
                  name <<- "EzAppSeuratVisiumHD"
                  appDefaults <<- rbind(
                    nfeatures = ezFrame(
                      Type = "numeric",
                      DefaultValue = 3000,
                      Description = "number of variable genes for SCT"
                    ),
                    npcs=ezFrame(Type="numeric", 
                                 DefaultValue=50,
                                 Description="The maximal dimensions to use for reduction"),
                    pcGenes=ezFrame(Type="charVector", 
                                    DefaultValue="", 
                                    Description="The genes used in supvervised clustering"),
                    SCT.regress.CellCycle=ezFrame(
                      Type = "logical", 
                      DefaultValue = FALSE,
                      Description="Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
                    ),
                    enrichrDatabase = ezFrame(
                      Type = "charVector", DefaultValue = "", Description="enrichR databases to search"
                    ),
                    DE.method=ezFrame(Type="charVector", 
                                      DefaultValue="wilcoxon", 
                                      Description="Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."),
                    min.pct = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.1,
                      Description = "Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations."
                    ),
                    logfc.threshold = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.25,
                      Description = "Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
                    ),
                    clusterResolution=ezFrame(Type="numeric", 
                                              DefaultValue=0.6,
                                              Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                    cellsFraction=ezFrame(Type="numeric", 
                                          DefaultValue=0, 
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
                      Description = "Low quality cells have less than \"nUMI\" reads. Only when applying fixed thresholds."
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
                    ),
                    pvalue_allMarkers = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.01,
                      Description = "pValue for marker detection"
                    ),
                    featSelectionMethod = ezFrame(
                      Type = "character",
                      DefaultValue = "STACAS",
                      Description = "use default method or black list genes"
                    ),
                    nfeatures = ezFrame(
                      Type = "numeric",
                      DefaultValue = 3000,
                      Description = "number of variable genes for PCA etc"
                    ),
                    pt.size.factor = ezFrame(
                      Type = "numeric",
                      DefaultValue = NA,
                      Description = "pt.size.factor for spatial plots"
                    ),
                    binSize = ezFrame(
                      Type = "character",
                      DefaultValue = "",
                      Description = "binning for Visium HD data"
                    ),
                    lambda = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.8,
                      Description = "BANKSY lambda: spatial weighting parameter (0-1). Larger values (0.8) find spatial domains; smaller values (0.2) perform cell typing."
                    ),
                    Niche_resolution = ezFrame(
                      Type = "numeric",
                      DefaultValue = 0.5,
                      Description = "Value of the Niche resolution parameter for BANKSY clustering, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
                    )
                  )
                }
              )
  )

ezMethodVisiumHDSeurat <- function(input=NA, output=NA, param=NA,
                                   htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Visium HD Seurat")))
  on.exit(setwd(cwd), add=TRUE)
  library(Seurat)
  library(scater)
  library(enrichR)
  library(future)
  library(BiocParallel)
  library(sf)


  if (param$cores > 1){
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    BPPARAM <- SerialParam() 
  }
  register(BPPARAM)
  plan("multicore", workers = param$cores)
  set.seed(38)
  future.seed = TRUE
  options(future.rng.onMisuse="ignore")
  options(future.globals.maxSize = param$ram*1024^3)
  
  # Handle segmented vs binned outputs
  if(grepl("segmented", param$binSize, ignore.case = TRUE)) {
    # Segmented outputs: use parent directory and bin.size = "polygons"
    dataDir <- input$getFullPaths("SpaceRangerDir")
    scData <- Load10X_Spatial(data.dir = dataDir, image.name = "tissue_hires_image.png", bin.size = "polygons", use.names=FALSE)
    sf <- scData@images[[1]]@scale.factors
    scData@images[[1]]@scale.factors$lowres <- sf$hires
    matrixPath <- file.path(dataDir, param$binSize, "filtered_feature_cell_matrix")
  } else {
    # Binned outputs: use standard approach with specific bin directory
    dataDir <- file.path(input$getFullPaths("SpaceRangerDir"), param$binSize)
    scData <- Load10X_Spatial(data.dir = dataDir, image.name = "tissue_hires_image.png", use.names=FALSE)
    sf <- scData@images[[1]]@scale.factors
    scData@images[[1]]@scale.factors$lowres <- sf$hires
    matrixPath <- file.path(dataDir, "filtered_feature_bc_matrix")
  }
  
  
  featInfo <- ezRead.table(paste0(matrixPath, "/features.tsv.gz"),
                           header = FALSE, row.names = NULL)
  colnames(featInfo) <- c("gene_id", "gene_name", "type")
  stopifnot(length(rownames(scData)) == nrow(featInfo))
  featInfo$isMito = grepl( "(?i)^MT-", featInfo$gene_name)
  featInfo$isRiboprot = grepl(  "(?i)^RPS|^RPL", featInfo$gene_name)
  rownames(featInfo) <- gsub("_", "-", uniquifyFeatureNames(ID=featInfo$gene_id, names=featInfo$gene_name)) 
  rownames(scData) <- rownames(featInfo)
  featInfo$isRibosomal <- getRibosomalFlag(featInfo$gene_id, annoFile=param$ezRef@refAnnotationFile)
  myAssay <- DefaultAssay(scData)
  scData[[myAssay]] <- AddMetaData(object = scData[[myAssay]], metadata = featInfo[rownames(scData), ])
  scData@meta.data$Sample <- input$getNames()
  if(!('Batch' %in% colnames(scData@meta.data))) {
    scData$Batch <- scData$Sample
  }
  
  param$nreads <- param$numis ## needed by qc script
  scData <- addCellQcToSeurat(scData, param=param, BPPARAM = BPPARAM, ribosomalGenes = featInfo[rownames(scData), "isRibosomal"])
  ## make image name unique
  stopifnot(length(names(scData@images)) == 1)
  names(scData@images) <- paste0(input$getNames(),'_S1')
  scData_unfiltered <- scData
  
  scData <- subset(scData_unfiltered, cells=which(scData_unfiltered$useCell)) # %>% head(n=1000))
  
  scData <- addCellCycleToSeurat(scData, param$refBuild, BPPARAM, assay = DefaultAssay(scData))
  
  scData <- NormalizeData(scData)
  scData <- FindVariableFeatures(scData)
  scData <- ScaleData(scData)
  if(nrow(scData@meta.data) < 50000){
    scData <- RunPCA(scData, npcs = 80)
    stopifnot(param$npcs <= 80)
    scData <- FindNeighbors(scData, dims = 1:param$npcs)
    scData <- FindClusters(scData, resolution = param$clusterResolution)
    scData <- RunUMAP(scData, reduction = "pca", reduction.name = "umap", return.model = T, dims = 1:param$npcs)
  } else {
    # we select 50,0000 cells and create a new 'sketch' assay
    scData <- SketchData(object = scData, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch", 
                         features = VariableFeatures(scData))
    # switch analysis to sketched cells
    DefaultAssay(scData) <- "sketch"
    
    # perform clustering workflow
    scData <- FindVariableFeatures(scData)
    scData <- ScaleData(scData)
    scData <- RunPCA(scData, assay = "sketch", reduction.name = "pca.sketch", npcs = 80)
    scData <- FindNeighbors(scData, assay = "sketch", reduction = "pca.sketch", dims = 1:param$npcs)
    scData <- FindClusters(scData, resolution = param$clusterResolution, cluster.name="seurat_clusters.sketched")
    #scData$seurat_clusters.sketched <- scData$seurat_clusters
    scData <- RunUMAP(scData, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:param$npcs)
    scData <- ProjectData(
      object = scData,
      assay = myAssay,
      full.reduction = "full.pca.sketch",
      sketched.assay = "sketch",
      sketched.reduction = "pca.sketch",
      umap.model = "umap.sketch",
      dims = 1:param$npcs,
      refdata = list(seurat_clusters.projected = "seurat_clusters.sketched")
    )
    scData$seurat_clusters.projected <- factor(scData$seurat_clusters.projected, levels(scData$seurat_clusters.sketched))
    
    # switch to full dataset
    ## TODO: do we need seurat_clusters.projected at all
    Idents(scData) <- "seurat_clusters.projected"
    scData$seurat_clusters <- Idents(scData)
    DefaultAssay(scData) <- myAssay
  }
  

  # get markers and annotations
  vars.to.regress = NULL
  posMarkers <- FindAllMarkers(object=scData, test.use = param$DE.method, only.pos=TRUE, latent.vars = vars.to.regress,min.pct = param$min.pct, 
                            return.thresh = param$pvalue_allMarkers, logfc.threshold = param$logfc.threshold)
  ## Significant markers
  posMarkers <- posMarkers[ ,c("gene","cluster","pct.1", "pct.2", "avg_log2FC","p_val_adj")]
  posMarkers$cluster <- as.factor(posMarkers$cluster)
  diff_pct = abs(posMarkers$pct.1-posMarkers$pct.2)
  posMarkers$diff_pct <- diff_pct
  posMarkers <- posMarkers[order(posMarkers$diff_pct, decreasing = TRUE),] %>% mutate_if(is.numeric, round, digits=30)
  posMarkers <- posMarkers[posMarkers$p_val_adj < param$pvalue_allMarkers,]
  rownames(posMarkers) <- NULL
  writexl::write_xlsx(posMarkers, path="posMarkers.xlsx")


  tryCatch({
    lambda <- ifelse(is.null(param$lambda), 0.8, as.numeric(param$lambda))
    niche_res <- ifelse(is.null(param$Niche_resolution), 0.5, as.numeric(param$Niche_resolution))
    
    myDefAssay <- DefaultAssay(scData)
    myIdents <- Idents(scData)
    scData <- RunBanksy(scData, lambda = lambda, assay = myDefAssay,
                           layer = "data", features = "variable", k_geom = 30,
                           verbose = FALSE)
    DefaultAssay(scData) <- "BANKSY"
    scData <- RunPCA(scData, assay = "BANKSY", reduction.name = "pca.banksy",
                     features = rownames(scData), npcs = 30, verbose = FALSE)
    scData <- FindNeighbors(scData, reduction = "pca.banksy", dims = 1:12,
                            verbose = FALSE)
    scData <- FindClusters(scData, cluster.name = "banksy_cluster",
                           resolution = niche_res, verbose = FALSE)
    
    # Niche markers
    Idents(scData) <- "banksy_cluster"
    posMarkersBanksy <- FindAllMarkers(scData, only.pos = TRUE, min.pct = 0.25,
                                       logfc.threshold = 0.25, verbose = FALSE)
    if (nrow(posMarkersBanksy) > 0) {
      posMarkersBanksy$diff_pct <- abs(posMarkersBanksy$pct.1 - posMarkersBanksy$pct.2)
      posMarkersBanksy <- posMarkersBanksy[order(posMarkersBanksy$diff_pct,
                                                 decreasing = TRUE), ]
    }
    writexl::write_xlsx(posMarkersBanksy, "posMarkersBanksy.xlsx")
    
    # Reset default assay
    DefaultAssay(scData) <- myDefAssay
    Idents(scData) <- myIdents
  }, error = function(e) {
    #writexl::write_xlsx(data.frame(), "posMarkersBanksy.xlsx")
  })
  
  
    
  
  ## generate template for manual cluster annotation -----
  ## we only deal with one sample
  stopifnot(length(input$getNames()) == 1)
  clusterInfos <- ezFrame(Sample=input$getNames(), Cluster=levels(Idents(scData)), ClusterLabel="")
  nTopMarkers <- 10
  topMarkers <- posMarkers %>% group_by(cluster) %>%
    slice_max(n = nTopMarkers, order_by = avg_log2FC)
  topMarkerString <- sapply(split(topMarkers$gene, topMarkers$cluster), paste, collapse=", ")
  clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
  clusterInfoFile <- "clusterInfos.xlsx"
  writexl::write_xlsx(clusterInfos, path=clusterInfoFile)
  
  #Save some results in external files
  # bulkSignalPerCluster <- AggregateExpression(scData, group.by = 'ident', assays=myAssay)[[1]]
  # bulkSignalPerCluster <- data.frame(GeneSymbol = rownames(scData), Count = bulkSignalPerCluster)
  bulkSignalPerSample <- AggregateExpression(scData, group.by = 'Sample', assays=myAssay)[[1]]
  bulkSignalPerSample <-  data.frame(GeneSymbol = rownames(scData), Count = bulkSignalPerSample)
  writexl::write_xlsx(bulkSignalPerSample, path="bulkSignalPerSample.xlsx")
  
  
  
  makeRmdReport(param=param, output=output, input=input, scData=scData, scData_unfiltered=scData_unfiltered,
                rmdFile = "VisiumHDSeurat.Rmd", reportTitle = paste0(param$name, ": ",  input$getNames()), use.qs2 = TRUE)
  
  gc()
  return("Success")
}




