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
                  appDefaults <<- rbind(
                      nfeatures = ezFrame(
                          Type = "numeric",
                          DefaultValue = 3000,
                          Description = "number of variable genes for SCT"
                      ),
                      npcs=ezFrame(Type="numeric", 
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
                                        ),
                                        pvalue_allMarkers = ezFrame(
                                            Type = "numeric",
                                            DefaultValue = 0.05,
                                            Description = "pValue for marker detection"
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
  library(Azimuth)
  #library(loupeR)
  
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
  
  res <- load10xSpatialData(input, param)
  scData <- res[['scData']]
  scDataRaw <- res[['scDataRaw']]
  param <- res[['param']]
  remove(res)
  
  scDataRes <- runBasicProcessing(scData,input, featInfo, param)
  cellsPerGeneFraction <- scDataRes[['cellsPerGeneFraction']]
  scData <- scDataRes[['scData']]
  scData.unfiltered <- scDataRes[['scData.unfiltered']]
  remove(scDataRes)
  
  
  #positive cluster markers
  clusterMarkers <- posClusterMarkers(scData, param$pvalue_allMarkers, param)
  clusterMarkers[['isSpatialMarker']] = FALSE 
  #spatially variable genes
  spatialMarkersList <- list()
  message('Find spatial Markers using markvariogram')
  res <- spatialMarkers(scData, selection.method = 'markvariogram')
  spatialMarkersList[['markvariogram']] <- data.frame(GeneSymbol = rownames(res), res, Method = 'Markvariogram')
  message('Find spatial Markers using moransi method')
  res <- spatialMarkers(scData, selection.method = 'moransi')
  spatialMarkersList[['moransi']] <- data.frame(GeneSymbol = rownames(res), res, Method = 'MoransI')
  spatialMarkers <- rbind(spatialMarkersList[['markvariogram']][,c('GeneSymbol', 'Rank','Method')], spatialMarkersList[['moransi']][,c('GeneSymbol', 'Rank','Method')])
  spatialMarkers <- spatialMarkers %>% spread(Method, Rank)
  spatialMarkers[['MeanRank']] <- apply(spatialMarkers[,c('Markvariogram','MoransI')],1,mean)
  spatialMarkers <- spatialMarkers[order(spatialMarkers$MeanRank),]
  
  spatialPosMarkers <- intersect(clusterMarkers$gene, spatialMarkers$GeneSymbol)
  clusterMarkers[which(clusterMarkers$gene %in% spatialPosMarkers), 'isSpatialMarker'] = TRUE
  
  # run Azimuth
  if (ezIsSpecified(param$Azimuth) && param$Azimuth != "none"){
      environment(MyDietSeurat) <- asNamespace('Seurat')
      assignInNamespace("DietSeurat", MyDietSeurat, ns = "Seurat")
      rna_assay <- CreateAssay5Object(counts = GetAssayData(scData, assay = 'Spatial', layer = 'counts'))
      scData[["RNA"]] <- scData[["Spatial"]] #rna_assay
      
      scDataAzi <- RunAzimuth(scData, param$Azimuth, assay="RNA") ## TODO support ADT
      
      ##Rename annotation levels if neccessary:
      colnames(scDataAzi@meta.data) <- sub('level_', 'l', colnames(scDataAzi@meta.data))
      if(param$Azimuth=='mousecortexref'){
          aziNames <- setdiff(colnames(scDataAzi@meta.data), colnames(scData@meta.data))
          aziResults <- data.frame(
              Azimuth.celltype.l1=scDataAzi@meta.data[ , grep("predicted.class$", aziNames, value=TRUE)],
              Azimuth.celltype.l2=scDataAzi@meta.data[ , grep("predicted.cluster$", aziNames, value=TRUE)],
              Azimuth.celltype.l3=scDataAzi@meta.data[ , grep("predicted.subclass$", aziNames, value=TRUE)],
              Azimuth.celltype.l4=scDataAzi@meta.data[ , grep("predicted.cross_species_cluster$", aziNames, value=TRUE)],
              row.names=colnames(scDataAzi))
      } else {
            aziNames <- setdiff(colnames(scDataAzi@meta.data), colnames(scData@meta.data))
            aziResults <- data.frame(
              Azimuth.celltype.l1=scDataAzi@meta.data[ , grep("l1$", aziNames, value=TRUE)],
              Azimuth.celltype.l2=scDataAzi@meta.data[ , grep("l2$", aziNames, value=TRUE)],
              Azimuth.celltype.l3=scDataAzi@meta.data[ , grep("l3$", aziNames, value=TRUE)],
              Azimuth.celltype.l4=scDataAzi@meta.data[ , grep("l4$", aziNames, value=TRUE)],
          row.names=colnames(scDataAzi))
      }
      scData[["RNA"]] <- NULL
      saveRDS(aziResults, "aziResults.rds")
      scData <- AddMetaData(scData, aziResults)
      remove(scDataAzi)
  } else {
      aziResults <- NULL
  }
  
  #Save some results in external files
  geneMeansPerCluster <- data.frame(AverageExpression(scData, group.by = 'ident')$SCT)
  geneMeans <-  data.frame(AverageExpression(scData, group.by = 'Sample')$SCT)
  
  dataFiles <- saveExternalFiles(list(cluster_markers=clusterMarkers, spatial_markers=data.frame(spatialMarkers), 
                                      gene_means = geneMeans, gene_means_per_cluster = geneMeansPerCluster))
  
  saveRDS(scData, "scData.rds")
  allCellsMeta <- scData.unfiltered@meta.data
  allCellsMeta$Sample <- allCellsMeta$Batch
  allCellsMeta$useCell <- !allCellsMeta$discard
  saveRDS(allCellsMeta, 'allCellsMeta.rds')
  saveRDS(scData.unfiltered, "scData.unfiltered.rds")
  saveRDS(param, "param.rds")
  saveRDS(input, "input.rds")
  saveRDS(scDataRaw, "scData.raw.rds")
  saveRDS(cellsPerGeneFraction, "cellsPerGeneFraction.rds")
  remove(scDataRaw, scData.unfiltered, scData)
  gc()
  makeRmdReport(dataFiles=dataFiles, output=output, rmdFile = "SpatialSeurat.Rmd", reportTitle = param$name) 
  return("Success")
}

runBasicProcessing <- function(scData, input, featInfo, param){
    scData@meta.data$Condition <- input$getColumn("Condition")
    scData@meta.data$Sample <- input$getNames()
    scData[["Spatial"]] <- AddMetaData(object = scData[["Spatial"]], metadata = featInfo[rownames(scData), ])
    
    scData_list <- filterCellsAndGenes(scData, param) # return sce objects filtered and unfiltered to show the QC metrics later in the rmd
    scData <- scData_list$scData
    scData.unfiltered <- scData_list$scData.unfiltered
    cellsPerGeneFraction <- scData_list$cellsPerGeneFraction
    rm(scData_list)
    scData <- addCellCycleToSeurat(scData, param$refBuild, BPPARAM, assay = 'Spatial')
    scData.unfiltered <- addCellCycleToSeurat(scData.unfiltered, param$refBuild, BPPARAM, assay = 'Spatial')
    scData <- seuratClusteringV3(scData, param, assay="Spatial")
    return(list(scData = scData, scData.unfiltered = scData.unfiltered, cellsPerGeneFraction = cellsPerGeneFraction))
}

MyDietSeurat <- function (object, layers = NULL, features = NULL, assays = NULL, 
                          dimreducs = NULL, graphs = NULL, misc = TRUE, counts = deprecated(), 
                          data = deprecated(), scale.data = deprecated(), ...) 
{
    CheckDots(...)
    dep.args <- c(counts = counts, data = data, scale.data = scale.data)
    for (lyr in names(x = dep.args)) {
        if (is_present(arg = dep.args[[lyr]])) {
            if (is.null(x = layers)) {
                layers <- unique(x = unlist(x = lapply(X = Assays(object = object), 
                                                       FUN = function(x) {
                                                           return(Layers(object = object[[x]]))
                                                       })))
            }
            deprecate_soft(when = "5.0.0", what = paste0("DietSeurat(", 
                                                         lyr, " = )"), with = "DietSeurat(layers = )")
            layers <- if (isTRUE(x = dep.args[[lyr]])) {
                c(layers, lyr)
            }
            else {
                Filter(f = function(x) x != lyr, x = layers)
            }
        }
    }
    object <- UpdateSlots(object = object)
    assays <- assays %||% Assays(object = object)
    assays <- intersect(x = assays, y = Assays(object = object))
    if (!length(x = assays)) {
        abort(message = "No assays provided were found in the Seurat object")
    }
    if (!DefaultAssay(object = object) %in% assays) {
        abort(message = "The default assay is slated to be removed, please change the default assay")
    }
    layers <- layers %||% assays
    layers <- .PropagateList(x = layers, names = assays)
    for (assay in names(x = layers)) {
        layers[[assay]] <- tryCatch(expr = Layers(object = object[[assay]], 
                                                  search = layers[[assay]]), error = function(...) {
                                                      return(character(length = 0L))
                                                  })
    }
    layers <- Filter(f = length, x = layers)
    if (!length(x = layers)) {
        abort(message = "None of the requested layers found")
    }
    for (assay in Assays(object = object)) {
        if (!(assay %in% assays)) {
            object[[assay]] <- NULL
            next
        }
        layers.rm <- setdiff(x = Layers(object = object[[assay]]), 
                             y = layers[[assay]])
        if (length(x = layers.rm)) {
            if (inherits(x = object[[assay]], what = "Assay") && 
                all(c("counts", "data") %in% layers.rm)) {
                abort(message = "Cannot remove both 'counts' and 'data' from v3 Assays")
            }
            for (lyr in layers.rm) {
                suppressWarnings(object <- tryCatch(expr = {
                    object[[assay]][[lyr]] <- NULL
                    object
                }, error = function(e) {
                    if (lyr == "data") {
                        object[[assay]][[lyr]] <- sparseMatrix(i = 1, 
                                                               j = 1, x = 1, dims = dim(object[[assay]][[lyr]]), 
                                                               dimnames = dimnames(object[[assay]][[lyr]]))
                    }
                    else {
                        slot(object = object[[assay]], name = lyr) <- new(Class = "dgCMatrix")
                    }
                    message("Converting layer ", lyr, " in assay ", 
                            assay, " to empty dgCMatrix")
                    object
                }))
            }
        }
        if (!is.null(x = features)) {
            features.assay <- intersect(x = features, y = rownames(x = object[[assay]]))
            if (!length(x = features.assay)) {
                warn(message = paste0("No features found in assay ", 
                                      sQuote(x = assay), ", removing..."))
                object[[assay]] <- NULL
                next
            }
            
            fixModel = FALSE
            if('SCTModel.list' %in% slotNames(object[[assay]])){
                myModel <- object[[assay]]@SCTModel.list$refmodel
                fixModel = TRUE
            }
            suppressWarnings(object[[assay]] <- subset(x = object[[assay]], 
                                                       features = features.assay))
            if(fixModel){
                object[[assay]]@SCTModel.list$refmodel <- myModel
            }
        }
    }
    if (!isTRUE(x = misc)) {
        slot(object = object, name = "misc") <- list()
    }
    all.objects <- .FilterObjects(object = object, classes.keep = c("DimReduc", 
                                                                    "Graph"))
    objects.to.remove <- all.objects[!all.objects %in% c(dimreducs, 
                                                         graphs)]
    for (ob in objects.to.remove) {
        object[[ob]] <- NULL
    }
    cells.keep <- list()
    for (assay in Assays(object = object)) {
        cells.keep[[assay]] <- colnames(x = object[[assay]])
    }
    cells.keep <- intersect(colnames(x = object), unlist(cells.keep))
    if (length(cells.keep) <- ncol(x = object)) {
        object <- subset(object, cells = cells.keep)
    }
    return(object)
}