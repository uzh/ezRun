###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSpatialSeuratHD <-
  setRefClass("EzAppSpatialSeuratHD",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSpatialSeuratHD
                  name <<- "EzAppSpatialSeuratHD"
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
                                        resolution=ezFrame(Type="numeric", 
                                                           DefaultValue=0.6,
                                                           Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                        cellsFraction=ezFrame(Type="numeric", 
                                                                DefaultValue=0.0001, 
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
                                            Type = "numeric",
                                            DefaultValue = 16,
                                            Description = "bin size for spatial data"
                                        )
                                    )
                }
              )
  )

ezMethodSpatialSeuratHD <- function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  library(Seurat)
  library(scater)
  library(enrichR)
  library(future)
  library(BiocParallel)
  
  
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
  dataDir <- file.path(dirname(input$getFullPaths("CountMatrix")), 'binned_outputs')
  if(param$binSize == 8){
      dataDir <- file.path(dataDir, 'square_008um')
  } else if(param$binSize >= 10 & param$binSize < 100){
      dataDir <- file.path(dataDir, paste0('square_0', param$binSize, 'um'))
  } else {
      stop("Only bin sizes of 8 and 16 or even numbers between 10-100 are supported if the parameter --custom-bin-size was used for SpaceRanger before") 
  }
  scData <- Load10X_Spatial(data.dir = dataDir)
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
  
  scDataRes <- runBasicProcessingHD(scData,input, featInfo, param, BPPARAM)
  
  cellsPerGeneFraction <- scDataRes[['cellsPerGeneFraction']]
  scData <- scDataRes[['scData']]
  scData.unfiltered <- scDataRes[['scData.unfiltered']]
  remove(scDataRes)
  
  # get markers and annotations
  anno <- getSpatialSeuratMarkersAndAnnotateHD(scData, param, BPPARAM)
  
  # save markers
  clusterMarkers <- anno$clusterMarkers
  writexl::write_xlsx(clusterMarkers, path="posMarkers.xlsx")
  
  #spatialMarkers <- anno$spatialMarkers
  #writexl::write_xlsx(spatialMarkers, path="spatialMarkers.xlsx")
  
  
  ## generate template for manual cluster annotation -----
  ## we only deal with one sample
  stopifnot(length(input$getNames()) == 1)
  clusterInfos <- ezFrame(Sample=input$getNames(), Cluster=levels(Idents(scData)), ClusterLabel="")
  nTopMarkers <- 10
  topMarkers <- clusterMarkers %>% group_by(cluster) %>%
      slice_max(n = nTopMarkers, order_by = avg_log2FC)
  topMarkerString <- sapply(split(topMarkers$gene, topMarkers$cluster), paste, collapse=", ")
  clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
  clusterInfoFile <- "clusterInfos.xlsx"
  writexl::write_xlsx(clusterInfos, path=clusterInfoFile)
  
  #Save some results in external files
  bulkSignalPerCluster <- AggregateExpression(scData, group.by = 'ident')[[1]]
  bulkSignalPerCluster <- data.frame(GeneSymbol = rownames(scData), Count = bulkSignalPerCluster)
  bulkSignalPerSample <- AggregateExpression(scData, group.by = 'Sample')[[1]]
  bulkSignalPerSample <-  data.frame(GeneSymbol = rownames(scData), Count = bulkSignalPerSample)
  
  
 ###Add TPM, add HVG annotation to base tables, enrichR links for HVG and top 1000 genes per sample
  dataFiles <- saveExternalFiles(list(bulkSignalPerSample = bulkSignalPerSample, bulkSignalPerCluster = bulkSignalPerCluster))
  
  allCellsMeta <- scData.unfiltered@meta.data
  allCellsMeta$Sample <- allCellsMeta$Batch
  allCellsMeta$useCell <- !allCellsMeta$discard
  ## make image name unique
  names(scData@images)[1] <- paste0(input$getNames(),'_S1')
  
 
  makeRmdReport(dataFiles=dataFiles, param=param, output=output, input=input, scData.unfiltered = scData.unfiltered, 
                scData=scData, allCellsMeta=allCellsMeta, cellsPerGeneFraction = cellsPerGeneFraction, 
                rmdFile = "SpatialSeuratHD.Rmd", reportTitle = paste0(param$name, ": ",  input$getNames()), use.qs2 = TRUE)
  
  remove(scData.unfiltered, scData)
  gc()
  return("Success")
}

runBasicProcessingHD <- function(scData, input, featInfo, param, BPPARAM){
    #scData$Condition <- unname(input$getColumn("Condition"))
    scData@meta.data$Sample <- input$getNames()
    myAssay <- DefaultAssay(scData)
    scData[[myAssay]] <- AddMetaData(object = scData[[DefaultAssay(scData)]], metadata = featInfo[rownames(scData), ])
    
    scData_list <- filterCellsAndGenesHD(scData, param, myAssay) # return sce objects filtered and unfiltered to show the QC metrics later in the rmd
    scData <- scData_list$scData
    scData.unfiltered <- scData_list$scData.unfiltered
    cellsPerGeneFraction <- scData_list$cellsPerGeneFraction
    rm(scData_list)
    scData <- addCellCycleToSeurat(scData, param$refBuild, BPPARAM, assay = myAssay)
    scData.unfiltered <- addCellCycleToSeurat(scData.unfiltered, param$refBuild, BPPARAM, assay = myAssay)
    
    
    scData <- NormalizeData(scData)
    scData <- FindVariableFeatures(scData)
    scData <- ScaleData(scData)
    if(nrow(scData@meta.data) < 50000){
        scData <- RunPCA(scData, npcs = 80)
        scData <- FindNeighbors(scData, dims = 1:param$npcs)
        scData <- FindClusters(scData, cluster.name = "seurat_cluster", resolution = param$resolution)
        scData <- RunUMAP(scData, reduction = "pca", reduction.name = "umap", return.model = T, dims = 1:param$npcs)
    }
    else {
        # we select 50,0000 cells and create a new 'sketch' assay
        scData <- SketchData(object = scData, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch", features = VariableFeatures(scData))
        
        # switch analysis to sketched cells
        DefaultAssay(scData) <- "sketch"
        
        # perform clustering workflow
        scData <- FindVariableFeatures(scData)
        scData <- ScaleData(scData)
        scData <- RunPCA(scData, assay = "sketch", reduction.name = "pca.sketch", npcs = 80)
        scData <- FindNeighbors(scData, assay = "sketch", reduction = "pca.sketch", dims = 1:param$npcs)
        scData <- FindClusters(scData, cluster.name = "seurat_cluster.sketched", resolution = 3)
        scData <- RunUMAP(scData, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:param$npcs)
        
        scData <- ProjectData(
            object = scData,
            assay = myAssay,
            full.reduction = "full.pca.sketch",
            sketched.assay = "sketch",
            sketched.reduction = "pca.sketch",
            umap.model = "umap.sketch",
            dims = 1:param$npcs,
            refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
        )
        
        # switch to full dataset
        Idents(scData) <- "seurat_cluster.projected"
        DefaultAssay(scData) <- myAssay
    }
    return(list(scData = scData, scData.unfiltered = scData.unfiltered, cellsPerGeneFraction = cellsPerGeneFraction))
}

filterCellsAndGenesHD <- function(scData, param, myAssay) {
    library(scater)
    library(Seurat)
    # Cells filtering
    scData <- PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
    scData <- PercentageFeatureSet(scData, "(?i)^RPS|^RPL", col.name = "percent_riboprot")
    att_nCounts <- paste0("nCount_", myAssay) 
    att_nGenes <- paste0("nFeature_", myAssay)
    
    if (is.na(param$nreads)) {
        qc.lib <- isOutlier(scData@meta.data[,att_nCounts], log = TRUE, nmads = param$nmad, type = "lower")
    } else {
        qc.lib <- scData@meta.data[,att_nCounts] < param$nreads
    }
    if (is.na(param$ngenes)) {
        qc.nexprs <- isOutlier(scData@meta.data[,att_nGenes], nmads = param$nmad, log = TRUE, type = "lower")
    } else {
        qc.nexprs <- scData@meta.data[,att_nGenes] < param$ngenes
    }
    if (is.na(param$perc_mito)) {
        qc.mito <- isOutlier(scData@meta.data[,"percent_mito"], nmads = param$nmad, type = "higher")
    } else {
        qc.mito <- scData@meta.data[,"percent_mito"] > param$perc_mito
    }
    if (is.na(param$perc_ribo )) {
        qc.ribo <- isOutlier(scData@meta.data[,"percent_riboprot"], nmads = param$nmad, type = "higher")
    } else {
        qc.ribo <- scData@meta.data[,"percent_riboprot"] > param$perc_ribo
    }
    
    discard <- qc.lib | qc.nexprs | qc.mito | qc.ribo
    scData$discard <- discard
    scData$qc.lib <- qc.lib
    scData$qc.nexprs <- qc.nexprs
    scData$qc.mito <- qc.mito
    scData$qc.ribo <- qc.ribo
    scData.unfiltered <- scData
    if(any(discard))
        scData <- scData[, -which(discard)]
    
    # Genes filtering
    ## remove low expressed genes
    num.cells <- param$cellsFraction * ncol(scData) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
    cellsPerGene <- Matrix::rowSums(GetAssayData(scData, layer="counts") >= param$nUMIs)
    is.expressed <- cellsPerGene >= num.cells
    cellsPerGeneFraction <- data.frame(frac = cellsPerGene/ncol(scData), row.names = rownames(cellsPerGene))
    scData <- scData[is.expressed,]
    return(list(scData.unfiltered = scData.unfiltered, scData = scData, cellsPerGeneFraction = cellsPerGeneFraction))
}

getSpatialSeuratMarkersAndAnnotateHD <- function(scData, param, BPPARAM){
    clusterMarkers <- posClusterMarkersSpatialHD(scData, param)
    clusterMarkers[['isSpatialMarker']] = NA
    return(list(clusterMarkers=clusterMarkers, spatialMarkers=NULL))
}

posClusterMarkersSpatialHD <- function(scData, param) {
    vars.to.regress = NULL
    markers <- FindAllMarkers(object=scData, test.use = param$DE.method, only.pos=TRUE, latent.vars = vars.to.regress,min.pct = param$min.pct, 
                              return.thresh = param$pvalue_allMarkers, logfc.threshold = param$logfc.threshold)
    ## Significant markers
    cm <- markers[ ,c("gene","cluster","pct.1", "pct.2", "avg_log2FC","p_val_adj")]
    cm$cluster <- as.factor(cm$cluster)
    diff_pct = abs(cm$pct.1-cm$pct.2)
    cm$diff_pct <- diff_pct
    cm <- cm[order(cm$diff_pct, decreasing = TRUE),] %>% mutate_if(is.numeric, round, digits=30)
    cm <- cm[cm$p_val_adj < param$pvalue_allMarkers,]
    rownames(cm) <- NULL
    return(cm)
}