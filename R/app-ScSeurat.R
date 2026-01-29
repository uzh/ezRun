###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppScSeurat <-
  setRefClass(
    "EzAppScSeurat",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodScSeurat
        name <<- "EzAppScSeurat"
        appDefaults <<- rbind(
          nfeatures = ezFrame(
            Type = "numeric",
            DefaultValue = 3000,
            Description = "number of variable genes for SCT"
          ),
          npcs = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "The maximal dimensions to use for reduction"
          ),
          pcGenes = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "The genes used in unsupervised clustering"
          ),
          SCT.regress.CellCycle = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
          ),
          DE.method = ezFrame(
            Type = "charVector",
            DefaultValue = "wilcoxon",
            Description = "Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."
          ),
          min.pct = ezFrame(
            Type = "numeric",
            DefaultValue = 0.1,
            Description = "Used in calculating cluster markers: The minimum fraction of cells in either of the two tested populations."
          ),
          min.diff.pct = ezFrame(
            Type = "numeric",
            DefaultValue = 0,
            Description = "Used for filtering cluster markers: The minimum difference of cell fraction of the two tested populations."
          ),
          logfc.threshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.25,
            Description = "Used in calculating cluster markers: Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells."
          ),
          pvalue_allMarkers = ezFrame(
            Type = "numeric",
            DefaultValue = 0.01,
            Description = "Used for filtering cluster markers: adjusted pValue threshold for marker detection"
          ),
          resolution = ezFrame(
            Type = "numeric",
            DefaultValue = 0.5,
            Description = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
          ),
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
          perc_riboprot = ezFrame(
            Type = "numeric",
            DefaultValue = Inf,
            Description = "Low quality cells have more than \"perc_ribo\" percent of ribosomal genes. Only when applying fixed thresholds."
          ),
          keepDoublets = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Whether we should keep cells suspected of being doublets. Set to TRUE only for QC purposes."
          ),
          maxEmptyDropPValue = ezFrame(
            Type = "numeric",
            DefaultValue = 1,
            Description = "filter droplets based on DropletUtils::emptyDrops method"
          ),
          cellsFraction = ezFrame(
            Type = "numeric",
            DefaultValue = 0,
            Description = "A gene will be kept if it is expressed in at least this percentage of cells"
          ),
          nUMIs = ezFrame(
            Type = "numeric",
            DefaultValue = 1,
            Description = "A gene will be kept if it has at least nUMIs in the fraction of cells specified before"
          ),
          nmad = ezFrame(
            Type = "numeric",
            DefaultValue = 3,
            Description = "Median absolute deviation (MAD) from the median value of each metric across all cells"
          ),
          filterByExpression = ezFrame(
            Type = "character",
            DefaultValue = FALSE,
            Description = "Keep cells according to specific gene expression. i.e. Set > 1 | Pkn3 > 1"
          ),
          estimateAmbient = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "estimate contamination with ambient RNA"
          ),
          controlSeqs = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "control sequences to add"
          ),
          enrichrDatabase = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "enrichR databases to search"
          ),
          geneCountModel = ezFrame(
            Type = "character",
            DefaultValue = "GeneFull_ExonOverIntron",
            Description = "(STARsolo Input Only) The gene count model, i.e. Solo features, to use from the previous step"
          ),
          computePathwayTFActivity = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "Whether we should compute pathway and TF activities."
          ),
          excludeGenes = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "file path to txt file with gene symbols to exclude from the analysis"
          ),
          sctype.enabled = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Enable scType automatic cell type annotation (human and mouse supported)"
          ),
          sctype.tissue = ezFrame(
            Type = "character",
            DefaultValue = "auto",
            Description = "Tissue type for scType annotation. Select 'auto' for automatic detection"
          ),
          sctype.confidence.threshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.25,
            Description = "Confidence threshold for scType annotation"
          ),
          AzimuthPanHuman = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Enable Azimuth Pan-Human neural network-based cell type annotation (HUMAN DATASETS ONLY)"
          ),
          AzimuthPanHuman.confidence.threshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.5,
            Description = "Confidence threshold for Azimuth Pan-Human annotation (0.0-1.0)"
          )
        )
      }
    )
  )

ezMethodScSeurat <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  library(HDF5Array)
  library(AUCell)
  library(GSEABase)
  library(SingleR)
  library(Seurat)
  library(tidyverse)
  library(scDblFinder)
  library(BiocParallel)
  library(scuttle)
  library(DropletUtils)
  library(enrichR)
  library(decoupleR)
  library(Azimuth)

  if (param$cores > 1) {
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    ## scDblFinder fails with many cells and MulticoreParam
    BPPARAM <- SerialParam()
  }
  register(BPPARAM)
  require(future)
  plan("multicore", workers = param$cores)
  set.seed(38)
  future.seed = TRUE
  options(future.rng.onMisuse = "ignore")
  options(future.globals.maxSize = param$ram * 1024^3)

  cwd <- getwd()
  setwdNew(basename(output$getColumn("SC Cluster Report")))
  on.exit(setwd(cwd), add = TRUE)

  cmDir <- input$getFullPaths("CountMatrix")
  if (file.exists(file.path(cmDir, param$geneCountModel))) {
    cmDir <- file.path(cmDir, param$geneCountModel)
  }
  if (grepl('h5$', cmDir[1])) {
    param$cellbender = TRUE
  } else {
    param$cellbender = FALSE
  }
  if (!param$cellbender) {
    cts <- Read10X(cmDir, gene.column = 1)
    featInfo <- ezRead.table(
      paste0(cmDir, "/features.tsv.gz"),
      header = FALSE,
      row.names = NULL
    )
  } else if (param$cellbender) {
    # Read the cellbender H5 file with Ensembl IDs as rownames
    cts <- Read10X_h5(cmDir, use.names = FALSE)

    ## Handle multi-modality H5 files (e.g., Cell Ranger ARC/Multi with RNA + ATAC)
    ## Extract only the Gene Expression matrix
    if (is.list(cts)) {
      cts <- cts$`Gene Expression`
    }

    # Get path for raw matrix from current input
    if ("UnfilteredCountMatrix" %in% input$colNames) {
      countRawMatrix <- input$getFullPaths("UnfilteredCountMatrix")
    } else {
      countRawMatrix <- file.path(dirname(cmDir), 'cellbender_raw_seurat.h5')
    }

    # Read features from H5 file to get proper gene_id and gene_name mapping
    h5_features <- rhdf5::h5read(cmDir, "matrix/features")

    # Create featInfo from H5 features
    featInfo <- data.frame(
      gene_id = h5_features$id,
      gene_name = h5_features$name,
      type = h5_features$feature_type,
      stringsAsFactors = FALSE
    )

    # Set rownames to gene_id for proper alignment
    rownames(featInfo) <- featInfo$gene_id

    param[['cellrangerDir']] <- dirname(cmDir)
    param[['cellrangerCountFiltDir']] <- dirname(cmDir)
    param[['cellrangerCountRawDir']] <- dirname(countRawMatrix)

    # Pass features directly - no temp file needed
    param[['featInfo_h5']] <- featInfo

    # Align features with matrix rownames
    matchingIds <- intersect(rownames(cts), rownames(featInfo))
    cts <- cts[matchingIds, ]
    featInfo <- featInfo[matchingIds, ]
  }

  featInfo <- featInfo[, 1:3] # in cases where additional column exist, e.g. CellRangerARC output
  colnames(featInfo) <- c("gene_id", "gene_name", "type")
  featInfo$isMito = grepl("(?i)^MT-", featInfo$gene_name)
  featInfo$isRiboprot = grepl("(?i)^RPS|^RPL", featInfo$gene_name)
  featInfo$isRibosomal <- getRibosomalFlag(
    featInfo$gene_id,
    annoFile = param$ezRef@refAnnotationFile
  )

  ## if we have feature barcodes we keep only the expression matrix
  if (is.list(cts)) {
    cts <- cts$`Gene Expression`
    featInfo <- featInfo[featInfo$type == "Gene Expression", ]
  }
  if (param$cellbender) {
    rownames(featInfo) <- featInfo$gene_id
    matchingIds <- intersect(rownames(cts), rownames(featInfo))
    cts <- cts[matchingIds, ]
    featInfo <- featInfo[matchingIds, ]
  }

  ## underscores in genenames will become dashes
  rownames(cts) <- rownames(featInfo) <- gsub(
    "_",
    "-",
    uniquifyFeatureNames(ID = featInfo$gene_id, names = featInfo$gene_name)
  )
  scData <- CreateSeuratObject(counts = cts[rowSums2(cts > 0) > 0, ])
  scData$Condition <- unname(input$getColumn("Condition"))
  scData@meta.data$Sample <- input$getNames()
  scData[["RNA"]] <- AddMetaData(
    object = scData[["RNA"]],
    metadata = featInfo[rownames(scData), ]
  )
  scData$cellBarcode <- sub(".*_", "", colnames(scData))
  scData <- addCellQcToSeurat(
    scData,
    param = param,
    BPPARAM = BPPARAM,
    ribosomalGenes = featInfo[rownames(scData), "isRibosomal"]
  )

  ## use empty drops to test for ambient
  if ("UnfilteredCountMatrix" %in% input$colNames) {
    rawDir <- input$getFullPaths("UnfilteredCountMatrix")
    if (file.exists(file.path(rawDir, param$geneCountModel))) {
      rawDir <- file.path(rawDir, param$geneCountModel)
    }
  } else {
    # DEPRECATED; all input datasets should specify the path to the unfiltered count matrix
    rawDir <- sub("filtered_", "raw_", cmDir)
  }
  if (file.exists(rawDir) && rawDir != cmDir) {
    if (param$cellbender) {
      rawCts <- Read10X_h5(
        file.path(dirname(cmDir), 'cellbender_raw_seurat.h5'),
        use.names = FALSE
      )
    } else {
      rawCts <- Read10X(rawDir, gene.column = 1)
    }
    if (is.list(rawCts)) {
      rawCts <- rawCts$`Gene Expression`
      rawCts <- rawCts[featInfo$gene_id, ]
    }

    if (
      ("SCDataOrigin" %in% input$colNames) &&
        input$getColumn("SCDataOrigin") == 'BDRhapsody'
    ) {
      rawCts <- rawCts[featInfo$gene_id, ]
    }

    scData$qc.empty <- FALSE
    if (!param$cellbender) {
      if (length(setdiff(rownames(rawCts), featInfo$gene_id)) > 0) {
        rawCts <- rawCts[featInfo$gene_id, ]
      }
      stopifnot(rownames(rawCts) == featInfo$gene_id)
      emptyStats <- emptyDrops(
        rawCts[!featInfo$isMito & !featInfo$isRiboprot, ],
        BPPARAM = BPPARAM,
        niters = 1e5
      )
      scData$negLog10CellPValue <- -log10(emptyStats[
        colnames(scData),
        "PValue"
      ])
      emptyStats <- emptyDrops(rawCts, BPPARAM = BPPARAM, niters = 1e5)
      scData$negLog10CellPValue <- pmin(
        scData$negLog10CellPValue,
        -log10(emptyStats[colnames(scData), "PValue"])
      )
      scData@meta.data$negLog10CellPValue[is.na(scData$negLog10CellPValue)] <- 0

      if (param$maxEmptyDropPValue < 1) {
        scData$qc.empty[
          scData$negLog10CellPValue < -log10(param$maxEmptyDropPValue)
        ] <- TRUE
        scData$useCell[scData$qc.empty] <- FALSE
      }
    }
    remove(rawCts)
  }
  allCellsMeta <- scData@meta.data

  scData <- subset(scData, cells = rownames(allCellsMeta)[allCellsMeta$useCell]) # %>% head(n=1000))

  ## remove lowly expressed genes
  num.cells <- param$cellsFraction * ncol(scData) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  cellsPerGene <- Matrix::rowSums(
    GetAssayData(scData, layer = "counts") >= param$nUMIs
  )
  is.expressed <- cellsPerGene >= num.cells
  cellsPerGeneFraction <- data.frame(
    frac = cellsPerGene / ncol(scData),
    row.names = rownames(cellsPerGene)
  )
  scData <- scData[is.expressed, ]

  if (ezIsSpecified(param$excludeGenes) && param$excludeGenes != '') {
    genesToExclude <- ezRead.table(
      param$excludeGenes,
      header = FALSE,
      row.names = NULL
    )
    genesToExclude <- unique(genesToExclude$V1)
    genesToKeep <- setdiff(rownames(scData), genesToExclude)
    scData <- subset(scData, features = genesToKeep)
  }

  ## Add Cell Cycle information to Seurat object as metadata columns
  scData <- addCellCycleToSeurat(scData, param$refBuild, BPPARAM)

  ## Get information on which variables to regress out in scaling/SCT
  scData <- seuratStandardSCTPreprocessing(scData, param)
  ## defaultAssay is now SCT
  scData <- seuratStandardWorkflow(
    scData,
    param,
    ident.name = "seurat_clusters"
  )

  # estimate ambient first
  if (ezIsSpecified(param$estimateAmbient) && param$estimateAmbient) {
    scData <- addAmbientEstimateToSeurat(scData, rawDir = rawDir, param = param)
  }

  # get markers and annotations
  anno <- getSeuratMarkersAndAnnotate(scData, param, BPPARAM = BPPARAM)

  # save markers
  markers <- anno$markers
  writexl::write_xlsx(markers, path = "posMarkers.xlsx")

  ## generate template for manual cluster annotation -----
  ## we only deal with one sample
  stopifnot(length(input$getNames()) == 1)
  clusterInfos <- ezFrame(
    Sample = input$getNames(),
    Cluster = levels(Idents(scData)),
    ClusterLabel = ""
  )
  if (!is.null(anno$aziResults)) {
    for (nm in grep("celltype", colnames(anno$aziResults), value = TRUE)) {
      cellCounts <- table(
        cluster = scData$seurat_clusters,
        sample = anno$aziResults[[nm]]
      )
      cellPerc <- sweep(cellCounts, 1, rowSums(cellCounts), "/")
      percMat <- as.matrix(cellPerc)
      newLabels <- apply(percMat, 1, function(x) {
        colnames(percMat)[x > 0.5]
      }) %>%
        unlist()
      clusterInfos[[nm]] <- clusterInfos$Cluster %>%
        as.character() %>%
        recode(!!!newLabels)
    }
  }
  if (!is.null(anno$singler.results)) {
    clusterInfos$SinglerCellType <- anno$singler.results$singler.results.cluster[
      clusterInfos$Cluster,
      "pruned.labels"
    ]
  }
  nTopMarkers <- 10
  topMarkers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = nTopMarkers, order_by = avg_log2FC)
  topMarkerString <- sapply(
    split(topMarkers$gene, topMarkers$cluster),
    paste,
    collapse = ", "
  )
  clusterInfos[["TopMarkers"]] <- topMarkerString[clusterInfos$Cluster]
  clusterInfoFile <- "clusterInfos.xlsx"
  writexl::write_xlsx(clusterInfos, path = clusterInfoFile)

  # scType Integration using proven wrapper
  if (
    ezIsSpecified(param$sctype.enabled) &&
      (param$sctype.enabled == TRUE || param$sctype.enabled == "true")
  ) {
    tryCatch(
      {
        futile.logger::flog.info("Starting scType cell type annotation...")

        # Load HGNChelper explicitly or create fallback
        hgnc_available <- suppressPackageStartupMessages(
          suppressWarnings(require("HGNChelper", quietly = TRUE))
        )

        if (!hgnc_available) {
          futile.logger::flog.warn(
            "HGNChelper not available, installing fallback function"
          )
          # Install required package first if possible
          tryCatch(
            {
              if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages(
                  "BiocManager",
                  repos = "https://cloud.r-project.org"
                )
              }
              BiocManager::install("HGNChelper", quiet = TRUE)
              hgnc_available <- require("HGNChelper", quietly = TRUE)
            },
            error = function(e) {
              futile.logger::flog.warn(
                "Failed to install HGNChelper: %s",
                e$message
              )
            }
          )
        }

        # If still not available, create fallback
        if (!hgnc_available) {
          futile.logger::flog.warn("Using fallback checkGeneSymbols function")
          # Create fallback function that matches HGNChelper interface
          checkGeneSymbols <<- function(
            x,
            unmapped.as.na = TRUE,
            map = NULL,
            species = "human"
          ) {
            data.frame(
              x = x,
              Suggested.Symbol = x,
              Approved = TRUE,
              stringsAsFactors = FALSE
            )
          }
        } else {
          futile.logger::flog.info("HGNChelper loaded successfully")
        }

        # Source the proven scType wrapper
        source(
          "https://raw.githubusercontent.com/kris-nader/sc-type/master/R/sctype_wrapper.R"
        )

        # Determine tissue type
        tissue_type <- param$sctype.tissue
        if (!ezIsSpecified(tissue_type) || tissue_type == "auto") {
          futile.logger::flog.info("Using default tissue type: Immune system")
          tissue_type <- "Immune system"
        }

        # Run scType using the wrapper function
        scData <- run_sctype(
          scData,
          known_tissue_type = tissue_type,
          custom_marker_file = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",
          name = "sctype_classification",
          plot = FALSE
        )

        # Create simple results for saving (no complex tables)
        sctype_results <- list(
          scData = scData,
          tissue_type = tissue_type
        )

        saveRDS(sctype_results, "sctype_results.rds")
        futile.logger::flog.info(
          "scType annotation completed successfully using wrapper"
        )
      },
      error = function(e) {
        futile.logger::flog.error("scType annotation failed: %s", e$message)
      }
    )
  } else {
    futile.logger::flog.info("scType annotation disabled, skipping...")
  }

  # Azimuth Pan-Human Integration using CloudAzimuth
  if (
    ezIsSpecified(param$AzimuthPanHuman) &&
      (param$AzimuthPanHuman == TRUE || param$AzimuthPanHuman == "true")
  ) {
    tryCatch(
      {
        futile.logger::flog.info("Starting Azimuth Pan-Human annotation...")

        # Load AzimuthAPI package explicitly (CloudAzimuth is in AzimuthAPI, not Azimuth)
        if (!require("AzimuthAPI", quietly = TRUE)) {
          futile.logger::flog.error("AzimuthAPI package not available")
          stop("AzimuthAPI package required for CloudAzimuth function")
        }

        # Verify RNA normalization before Azimuth Pan-Human annotation
        if (
          !"data" %in% names(scData[["RNA"]]@layers) ||
            is.null(scData[["RNA"]]@layers[["data"]])
        ) {
          futile.logger::flog.info(
            "RNA normalization not found. Running NormalizeData for Azimuth Pan-Human annotation..."
          )
          scData <- NormalizeData(scData, assay = "RNA")
        }

        # Run CloudAzimuth - this handles everything automatically
        scData <- CloudAzimuth(scData)

        # Restore original seurat_clusters as default Idents (CloudAzimuth changes this)
        Idents(scData) <- scData$seurat_clusters
        futile.logger::flog.info(
          "Restored seurat_clusters as default Idents after CloudAzimuth"
        )

        # Create simple results for saving (no complex tables)
        azimuth_results <- list(
          scData = scData
        )

        saveRDS(azimuth_results, "azimuth_results.rds")
        futile.logger::flog.info(
          "Azimuth Pan-Human annotation completed successfully"
        )
      },
      error = function(e) {
        futile.logger::flog.error(
          "Azimuth Pan-Human annotation failed: %s",
          e$message
        )
      }
    )
  } else {
    futile.logger::flog.info(
      "Azimuth Pan-Human annotation disabled, skipping..."
    )
  }

  qs2::qs_save(scData, "scData.qs2", nthreads = param$cores)

  makeRmdReport(
    param = param,
    output = output,
    scData = scData,
    allCellsMeta = allCellsMeta,
    cellsPerGeneFraction = cellsPerGeneFraction,
    enrichRout = anno$enrichRout,
    cells.AUC = anno$cells.AUC,
    singler.results = anno$singler.results,
    aziResults = anno$aziResults,
    pathwayActivity = anno$pathwayActivity,
    TFActivity = anno$TFActivity,
    cellxgeneResults = anno$cellxgeneResults,
    rmdFile = "ScSeurat.Rmd",
    reportTitle = paste0(param$name, ": ", input$getNames())
  )
  #remove no longer used objects
  rm(scData)
  gc()
  return("Success")
}

addCellQcToSeurat <- function(
  scData,
  param = NULL,
  BPPARAM = NULL,
  ribosomalGenes = NULL
) {
  library(scater)

  # Cells filtering
  ## TODO: extract mito / riboprot and ribosomal genes from assay annotation
  scData <- PercentageFeatureSet(scData, "(?i)^MT-", col.name = "percent_mito")
  scData <- PercentageFeatureSet(
    scData,
    "(?i)^RPS|^RPL",
    col.name = "percent_riboprot"
  )
  if (!is.null(ribosomalGenes)) {
    scData <- PercentageFeatureSet(
      scData,
      features = ribosomalGenes,
      col.name = "percent_ribosomal"
    )
  }
  att_nCounts <- paste0("nCount_", DefaultAssay(scData))
  att_nGenes <- paste0("nFeature_", DefaultAssay(scData))

  if (!ezIsSpecified(param$nreads)) {
    scData$qc.lib <- isOutlier(
      scData@meta.data[, att_nCounts],
      log = TRUE,
      nmads = param$nmad,
      type = "lower"
    )
  } else {
    scData$qc.lib <- scData@meta.data[, att_nCounts] < param$nreads
  }
  if (!ezIsSpecified(param$ngenes)) {
    scData$qc.nexprs <- isOutlier(
      scData@meta.data[, att_nGenes],
      nmads = param$nmad,
      log = TRUE,
      type = "lower"
    )
  } else {
    scData$qc.nexprs <- scData@meta.data[, att_nGenes] < param$ngenes
  }
  if (!ezIsSpecified(param$perc_mito)) {
    scData$qc.mito <- isOutlier(
      scData@meta.data[, "percent_mito"],
      subset = !is.na(scData@meta.data[, "percent_mito"]),
      nmads = param$nmad,
      type = "higher"
    )
  } else {
    scData$qc.mito <- scData@meta.data[, "percent_mito"] > param$perc_mito
  }
  if (!ezIsSpecified(param$perc_riboprot)) {
    scData$qc.riboprot <- isOutlier(
      scData@meta.data[, "percent_riboprot"],
      subset = !is.na(scData@meta.data[, "percent_riboprot"]),
      nmads = param$nmad,
      type = "higher"
    )
  } else {
    scData$qc.riboprot <- scData@meta.data[, "percent_riboprot"] >
      as.numeric(param$perc_riboprot)
  }

  scData$useCell <- !(scData$qc.lib |
    scData$qc.nexprs |
    scData$qc.mito |
    scData$qc.riboprot)
  #TODO: consider also ribosomal????

  if (DefaultAssay(scData) == "RNA") {
    set.seed(38)
    doubletsInfo <- scDblFinder(
      GetAssayData(scData, layer = "counts")[, scData$useCell],
      returnType = "table",
      clusters = TRUE,
      BPPARAM = BPPARAM
    )
    scData$doubletScore <- doubletsInfo[colnames(scData), "score"]
    scData$doubletClass <- doubletsInfo[colnames(scData), "class"]
    scData$qc.doublet <- scData$doubletClass %in% "doublet"
    if (ezIsSpecified(param$keepDoublets) && param$keepDoublets) {
      futile.logger::flog.info("Keeping doublets...")
    } else {
      scData$useCell <- scData$useCell & scData$doubletClass %in% "singlet"
    }
  }
  return(scData)
}

querySignificantClusterAnnotationEnrichR <- function(
  genesPerCluster,
  dbs,
  overlapGeneCutOff = 3,
  adjPvalueCutOff = 0.001,
  reportTopN = 5,
  keepGenes = FALSE
) {
  enrichRout <- list()
  columnsToKeep <- c(
    "Term",
    "Cluster",
    "Overlap",
    "OverlapGenesN",
    "Adjusted.P.value",
    "Odds.Ratio",
    "Combined.Score"
  )
  if (keepGenes) {
    columnsToKeep <- c(columnsToKeep, "Genes")
  }
  for (cluster in unique(names(genesPerCluster))) {
    # Check if gene list is empty or contains only empty/NA values
    genes <- as.character(genesPerCluster[[cluster]])
    genes <- genes[!is.na(genes) & genes != ""]

    if (length(genes) == 0) {
      futile.logger::flog.warn(
        "Cluster %s has no genes for enrichment analysis, skipping",
        cluster
      )
      next
    }

    enriched <- enrichr(genes, dbs)

    for (db in names(enriched)) {
      enriched_db <- enriched[[db]]
      if (nrow(enriched_db) > 0 && colnames(enriched_db)[1] == "Term") {
        enriched_db$OverlapGenesN <- sub("/.*", "", enriched_db$Overlap) %>%
          as.numeric()
        enriched_db$Cluster <- cluster
        enriched_db <- enriched_db %>%
          filter(., Adjusted.P.value < adjPvalueCutOff) %>%
          filter(., OverlapGenesN > overlapGeneCutOff) %>%
          head(reportTopN)
        enrichRout[[cluster]][[db]] <- enriched_db[, columnsToKeep]
      }
    }
  }
  return(enrichRout)
}


computeTFActivityAnalysis <- function(cells, species) {
  species <- tolower(species)
  # Retrieve prior knowledge network.
  network <- decoupleR::get_dorothea(
    organism = species,
    levels = c("A", "B", "C")
  )

  # Run weighted means algorithm.
  activities <- decoupleR::run_wmean(
    mat = as.matrix(GetAssayData(cells)),
    network = network,
    .source = "source",
    .targe = "target",
    .mor = "mor",
    times = 100,
    minsize = 5
  )

  return(activities)
}


computePathwayActivityAnalysis <- function(cells, species) {
  species <- tolower(species)
  # Retrieve prior knowledge network.
  network <- decoupleR::get_progeny(organism = species)

  # Run weighted means algorithm.
  activities <- decoupleR::run_wmean(
    mat = as.matrix(GetAssayData(cells)),
    network = network,
    .source = "source",
    .targe = "target",
    .mor = "weight",
    times = 100,
    minsize = 5
  )

  return(activities)
}


getRibosomalFlag <- function(gene_id, annoFile) {
  geneAnnoFile <- sub("byTranscript", "byGene", annoFile)
  if (file.exists(geneAnnoFile)) {
    geneAnno <- ezRead.table(geneAnnoFile)
    if (any(geneAnno$type == "rRNA")) {
      isRibosomal <- geneAnno[gene_id, "type"] %in% "rRNA"
      if (any(isRibosomal)) {
        return(isRibosomal)
      } else {
        return(NULL)
      }
    }
  }
}
