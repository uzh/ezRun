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
          nUMI = ezFrame(
            Type = "numeric",
            DefaultValue = Inf,
            Description = "Low quality cells have less than \"nUMI\" UMIs. Only when applying fixed thresholds."
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
          geneMinUMI = ezFrame(
            Type = "numeric",
            DefaultValue = 1,
            Description = "A gene will be kept if it has at least this many UMIs in the fraction of cells specified before"
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
          ),
          CyteTypeR = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Enable CyteTypeR AI-powered cell type annotation"
          ),
          CyteTypeR.apiKey = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Nygen Analytics CyteType API key"
          ),
          CyteTypeR.studyContext = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Biological context for annotation"
          ),
          mLLMCelltype = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Enable mLLMCelltype cluster annotation on the FGCZ-internal vLLM server (no data leaves FGCZ)"
          ),
          mLLMCelltype.tissue = ezFrame(
            Type = "character",
            DefaultValue = "auto",
            Description = "Tissue context for mLLMCelltype. 'auto' falls back to sctype.tissue, then to the reference species alone"
          )
          ## NOTE: the vLLM endpoint is deliberately NOT a parameter. Every
          ## appDefault is written to the delivered parameters.tsv, so declaring
          ## it here would publish the internal host to every customer. It lives
          ## as a constant in fgczVllmEndpoint(), same convention as
          ## app-fastQC.R's AI_ENDPOINT.
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
  library(decoupleR)

  if (param$cores > 1) {
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    ## scDblFinder fails with many cells and MulticoreParam
    BPPARAM <- SerialParam()
  }
  ## Pin BLAS/OpenMP to one thread before forking (MulticoreParam/future) to
  ## avoid the fork-in-multithreaded-process deadlock (e.g. AUCell labeling).
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
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
    # fill = TRUE: multiome (CR ARC) features.tsv.gz mixes 6-field GEX rows
    # (id, name, "Gene Expression", chr, start, end) with 3-field Peaks rows
    # (id, name, "Peaks"). Without fill, read.table errors at the first short
    # row. Load featInfo first so we can detect multiome before calling Read10X.
    featInfo <- ezRead.table(
      paste0(cmDir, "/features.tsv.gz"),
      header = FALSE,
      row.names = NULL,
      fill = TRUE
    )
    # Detect multiome (Peaks present alongside Gene Expression). Seurat's
    # Read10X uses data.table::fread on features.tsv.gz without fill, which
    # truncates at the first short row and crashes during dimname assignment
    # ("length of Dimnames[[1]] (1) is not equal to Dim[1] (...)"). Side-step
    # by reading from the sibling .h5 — Read10X_h5 handles multi-feature_type
    # input natively (returns a named list keyed by feature_type). The
    # subsequent `is.list(cts)` filter at line ~297 keeps only Gene Expression.
    matrix_h5 <- paste0(cmDir, ".h5")
    has_peaks <- "Peaks" %in% as.character(featInfo[[3]])
    if (has_peaks && file.exists(matrix_h5)) {
      cts <- Read10X_h5(matrix_h5, use.names = FALSE)
    } else {
      cts <- Read10X(cmDir, gene.column = 1)
    }
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
    BPPARAM = BPPARAM
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
      # Same multiome guard as the filtered-matrix branch above: when the
      # unfiltered features.tsv.gz contains Peaks rows AND a sibling .h5
      # exists, prefer Read10X_h5 — Seurat::Read10X chokes on multiome
      # features.tsv.gz via internal data.table::fread without fill.
      raw_h5 <- paste0(rawDir, ".h5")
      raw_features <- file.path(rawDir, "features.tsv.gz")
      raw_has_peaks <- FALSE
      if (file.exists(raw_features)) {
        raw_feat <- tryCatch(
          ezRead.table(raw_features, header = FALSE, row.names = NULL,
                       fill = TRUE),
          error = function(e) NULL
        )
        if (!is.null(raw_feat) && ncol(raw_feat) >= 3) {
          raw_has_peaks <- "Peaks" %in% as.character(raw_feat[[3]])
        }
      }
      if (raw_has_peaks && file.exists(raw_h5)) {
        rawCts <- Read10X_h5(raw_h5, use.names = FALSE)
      } else {
        rawCts <- Read10X(rawDir, gene.column = 1)
      }
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
      # Skip emptyDrops if raw matrix has no additional (empty) barcodes
      # This happens with CellRanger Multi per_sample_outs where raw == filtered
      if (ncol(rawCts) > ncol(scData)) {
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
        scData@meta.data$negLog10CellPValue[is.na(
          scData$negLog10CellPValue
        )] <- 0

        if (param$maxEmptyDropPValue < 1) {
          scData$qc.empty[
            scData$negLog10CellPValue < -log10(param$maxEmptyDropPValue)
          ] <- TRUE
          scData$useCell[scData$qc.empty] <- FALSE
        }
      } else {
        futile.logger::flog.warn(
          "emptyDrops skipped: raw matrix has no additional empty barcodes (e.g. CellRanger Multi per_sample output)"
        )
        scData$negLog10CellPValue <- 0
      }
    }
    remove(rawCts)
  }
  allCellsMeta <- scData@meta.data

  scData <- subset(scData, cells = rownames(allCellsMeta)[allCellsMeta$useCell]) # %>% head(n=1000))

  ## remove lowly expressed genes
  num.cells <- param$cellsFraction * ncol(scData) # if we expect at least one rare subpopulation of cells, we should decrease the percentage of cells
  cellsPerGene <- Matrix::rowSums(
    GetAssayData(scData, layer = "counts") >= param$geneMinUMI
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

        # Save only the tissue type (scData already contains the annotation)
        sctype_results <- list(
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

  # mLLMCelltype annotation on the FGCZ-internal vLLM server.
  # Same marker-gene-to-LLM idea as CyteTypeR, but the endpoint is on-site, so
  # no expression data leaves FGCZ and no API key is needed.
  if (
    ezIsSpecified(param$mLLMCelltype) &&
      (param$mLLMCelltype == TRUE || param$mLLMCelltype == "true")
  ) {
    tryCatch(
      {
        futile.logger::flog.info("Starting mLLMCelltype cell type annotation...")

        if (!requireNamespace("mLLMCelltype", quietly = TRUE)) {
          stop(
            "mLLMCelltype package required. Install with: ",
            "pak::pkg_install('cafferychen777/mLLMCelltype/R')"
          )
        }

        ## The tissue string is the single biggest lever on label quality, far
        ## bigger than any model or prompting choice. On a FACS-sorted B-cell
        ## dataset the fallback context ("Human Immune system") produced
        ## confident monocyte / Treg / ILC2 / stromal calls for clusters that
        ## were all B cells - 4 of 14 right. Passing "human sorted B cells"
        ## against the same markers gave 13 of 14. So record whether the value
        ## was actually supplied, and warn in the report when it was not.
        tissueIsAuto <- !(
          ezIsSpecified(param$mLLMCelltype.tissue) &&
            param$mLLMCelltype.tissue != "auto"
        )
        tissue <- if (!tissueIsAuto) {
          param$mLLMCelltype.tissue
        } else if (
          ezIsSpecified(param$sctype.tissue) && param$sctype.tissue != "auto"
        ) {
          param$sctype.tissue
        } else {
          ""
        }
        if (tissueIsAuto) {
          futile.logger::flog.warn(
            paste0(
              "mLLMCelltype.tissue is 'auto', falling back to '%s'. Labels are ",
              "materially worse without a specific tissue/sort description."
            ),
            tissue
          )
        }

        mllm <- annotateClustersWithMLLMCelltype(
          topMarkers = topMarkers,
          tissueName = trimws(paste(getSpecies(param$refBuild), tissue)),
          vllmUrl = fgczVllmEndpoint()
        )

        ## unname() is required for Seurat v5 metadata assignment
        scData$mLLMCelltype <- unname(
          mllm$clusterLabels[as.character(scData$seurat_clusters)]
        )

        clusterInfos[["mLLMCelltype"]] <- unname(
          mllm$clusterLabels[as.character(clusterInfos$Cluster)]
        )
        clusterInfos[["mLLMCelltypeMarkers"]] <- unname(
          mllm$clusterSupport[as.character(clusterInfos$Cluster)]
        )
        writexl::write_xlsx(clusterInfos, path = clusterInfoFile)

        saveRDS(
          list(
            ## model name only - the endpoint is never written to a
            ## delivered file (this .rds ships to gstore with the report)
            model = mllm$model,
            tissue_name = mllm$tissueName,
            tissue_is_auto = tissueIsAuto,
            clusterLabels = mllm$clusterLabels,
            clusterSupport = mllm$clusterSupport
          ),
          "mllmcelltype_results.rds"
        )
        futile.logger::flog.info(
          "mLLMCelltype annotation completed successfully"
        )
      },
      error = function(e) {
        futile.logger::flog.error(
          "mLLMCelltype annotation failed: %s",
          redactVllmEndpoint(e$message)
        )
      }
    )
  } else {
    futile.logger::flog.info("mLLMCelltype annotation disabled, skipping...")
  }

  # Azimuth Pan-Human Integration using CloudAzimuth
  if (
    ezIsSpecified(param$AzimuthPanHuman) &&
      (param$AzimuthPanHuman == TRUE || param$AzimuthPanHuman == "true")
  ) {
    tryCatch(
      {
        futile.logger::flog.info("Starting Azimuth Pan-Human annotation...")

        # CloudAzimuth is in AzimuthAPI, not Azimuth.
        if (!requireNamespace("AzimuthAPI", quietly = TRUE)) {
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
        scData <- AzimuthAPI::CloudAzimuth(scData)

        # Restore original seurat_clusters as default Idents (CloudAzimuth changes this)
        Idents(scData) <- scData$seurat_clusters
        futile.logger::flog.info(
          "Restored seurat_clusters as default Idents after CloudAzimuth"
        )

        # Save a marker that Azimuth completed (scData already contains the annotation)
        azimuth_results <- list(completed = TRUE)

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

  # CyteTypeR AI-powered Annotation
  if (
    ezIsSpecified(param$CyteTypeR) &&
      (param$CyteTypeR == TRUE || param$CyteTypeR == "true")
  ) {
    tryCatch(
      {
        futile.logger::flog.info("Starting CyteTypeR cell type annotation...")

        if (!require("CyteTypeR", quietly = TRUE)) {
          futile.logger::flog.error("CyteTypeR package not available")
          stop(
            "CyteTypeR package required. ",
            "Install from: https://github.com/NygenAnalytics/CyteTypeR"
          )
        }

        # Resolve API key: parameter > environment variable
        api_key <- NULL
        if (ezIsSpecified(param$CyteTypeR.apiKey)) {
          api_key <- param$CyteTypeR.apiKey
        } else if (nzchar(Sys.getenv("CYTETYPE_API_KEY", ""))) {
          api_key <- Sys.getenv("CYTETYPE_API_KEY")
        }

        # Prepare markers (reuse existing markers from getSeuratMarkersAndAnnotate)
        cytetype_markers <- markers |>
          dplyr::mutate(cluster = as.character(cluster)) |>
          dplyr::filter(!is.na(cluster), p_val_adj < 0.05, avg_log2FC > 0.5) |>
          dplyr::group_by(cluster) |>
          dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE)

        # Work on a copy to avoid mutating scData before CyteTypeR succeeds
        scData_ct <- scData

        # Ensure seurat_clusters is character (local copy only)
        scData_ct$seurat_clusters <- as.character(scData_ct$seurat_clusters)

        # Remove cells with NA cluster
        cells_with_cluster <- colnames(scData_ct)[
          !is.na(scData_ct$seurat_clusters)
        ]
        if (length(cells_with_cluster) < ncol(scData_ct)) {
          futile.logger::flog.warn(
            "Removing %d cells with NA cluster assignment for CyteTypeR",
            ncol(scData_ct) - length(cells_with_cluster)
          )
          scData_ct <- subset(scData_ct, cells = cells_with_cluster)
        }

        # Handle case-insensitive duplicate columns (DuckDB requirement)
        cols <- colnames(scData_ct@meta.data)
        dupes <- cols[duplicated(tolower(cols))]
        if (length(dupes) > 0) {
          futile.logger::flog.warn(
            "Removing duplicate metadata columns for CyteTypeR: %s",
            paste(dupes, collapse = ", ")
          )
          scData_ct@meta.data[dupes] <- NULL
        }

        # Ensure clusters have enough markers (PrepareCyteTypeR needs >= 5 per cluster)
        marker_counts <- cytetype_markers |>
          dplyr::count(cluster) |>
          dplyr::filter(n >= 5)
        valid_clusters <- marker_counts$cluster
        seurat_clusters_all <- unique(scData_ct$seurat_clusters)
        excluded_clusters <- setdiff(seurat_clusters_all, valid_clusters)
        if (length(excluded_clusters) > 0) {
          futile.logger::flog.warn(
            "Clusters excluded from CyteTypeR (no/insufficient markers): %s",
            paste(excluded_clusters, collapse = ", ")
          )
          scData_ct <- subset(
            scData_ct,
            seurat_clusters %in% valid_clusters
          )
          cytetype_markers <- cytetype_markers |>
            dplyr::filter(cluster %in% valid_clusters)
        }

        # Prepare data for CyteTypeR
        prepped_data <- CyteTypeR::PrepareCyteTypeR(
          scData_ct,
          cytetype_markers,
          n_top_genes = 15,
          group_key = "seurat_clusters",
          aggregate_metadata = TRUE,
          coordinates_key = "umap"
        )

        # Fix NA entry in clusterMetadata (known PrepareCyteTypeR bug)
        if (any(is.na(names(prepped_data$clusterMetadata)))) {
          prepped_data$clusterMetadata <- prepped_data$clusterMetadata[
            !is.na(names(prepped_data$clusterMetadata))
          ]
        }

        # Auto-generate study_context if not provided
        species <- getSpecies(param$refBuild)
        study_context <- if (ezIsSpecified(param$CyteTypeR.studyContext)) {
          param$CyteTypeR.studyContext
        } else {
          tissue_hint <- if (
            ezIsSpecified(param$sctype.tissue) &&
              param$sctype.tissue != "auto"
          ) {
            param$sctype.tissue
          } else {
            "unspecified tissue"
          }
          paste0(
            species, " ", tissue_hint,
            " single-cell RNA-seq data with ",
            length(unique(scData$seurat_clusters)), " clusters and ",
            ncol(scData), " cells."
          )
        }

        # Submit annotation job (use copy to avoid mutating scData on failure)
        scData_ct <- CyteTypeR::CyteTypeR(
          obj = scData_ct,
          prepped_data = prepped_data,
          study_context = study_context,
          auth_token = api_key,
          require_artifacts = FALSE,
          metadata = list(
            title = paste0("ScSeurat CyteTypeR - ", input$getNames()),
            run_label = paste0("sushi_", format(Sys.Date(), "%Y%m%d")),
            experiment_name = param$name
          )
        )

        # Extract results from the copy
        cytetype_results <- scData_ct@misc$cytetype_results

        if (!is.null(cytetype_results) && nrow(cytetype_results) > 0) {
          # Map annotations to original scData (CRITICAL: use unname() for Seurat v5)
          ct_map <- setNames(
            cytetype_results$annotation,
            cytetype_results$clusterId
          )
          ct_vals <- unname(ct_map[as.character(scData$seurat_clusters)])
          ct_vals[is.na(ct_vals)] <- "Unassigned"
          scData$CyteTypeR_annotation <- ct_vals

          ct_granular_map <- setNames(
            cytetype_results$granularAnnotation,
            cytetype_results$clusterId
          )
          ct_gran_vals <- unname(ct_granular_map[as.character(scData$seurat_clusters)])
          ct_gran_vals[is.na(ct_gran_vals)] <- "Unassigned"
          scData$CyteTypeR_granular <- ct_gran_vals

          ct_onto_map <- setNames(
            cytetype_results$ontologyTerm,
            cytetype_results$clusterId
          )
          ct_onto_vals <- unname(ct_onto_map[as.character(scData$seurat_clusters)])
          ct_onto_vals[is.na(ct_onto_vals)] <- "Unassigned"
          scData$CyteTypeR_ontology <- ct_onto_vals

          # Transfer job details for report link
          scData@misc$cytetype_results <- cytetype_results
          scData@misc$cytetype_jobDetails <- scData_ct@misc$cytetype_jobDetails

          # Save results for Rmd template
          saveRDS(cytetype_results, "cytetype_results.rds")

          # Add CyteTypeR to clusterInfos.xlsx
          if (file.exists(clusterInfoFile)) {
            clusterInfos_update <- readxl::read_xlsx(clusterInfoFile)
            clusterInfos_update$CyteTypeR <- unname(
              ct_map[as.character(clusterInfos_update$Cluster)]
            )
            writexl::write_xlsx(clusterInfos_update, path = clusterInfoFile)
          }

          futile.logger::flog.info(
            "CyteTypeR annotation completed: %d cell types identified",
            length(unique(cytetype_results$annotation))
          )
        } else {
          futile.logger::flog.warn("CyteTypeR returned no results")
        }
        rm(scData_ct)
        gc()
      },
      error = function(e) {
        futile.logger::flog.error(
          "CyteTypeR annotation failed: %s",
          e$message
        )
      }
    )
  } else {
    futile.logger::flog.info("CyteTypeR annotation disabled, skipping...")
  }

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
  BPPARAM = NULL
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

  scData <- PercentageFeatureSet(
    scData,
    "(?i)^HB[^P]",
    col.name = "percent_hb"
  )

  att_nCounts <- paste0("nCount_", DefaultAssay(scData))
  att_nGenes <- paste0("nFeature_", DefaultAssay(scData))

  # Handle QC filtering - use fixed thresholds if specified,

  # otherwise use MAD-based filtering if nmad is specified
  use_mad <- ezIsSpecified(param$nmad)

  if (!ezIsSpecified(param$nUMI)) {
    if (use_mad) {
      scData$qc.lib <- isOutlier(
        scData@meta.data[, att_nCounts],
        log = TRUE,
        nmads = param$nmad,
        type = "lower"
      )
    } else {
      scData$qc.lib <- FALSE # No filtering if neither threshold nor nmad set
    }
  } else {
    scData$qc.lib <- scData@meta.data[, att_nCounts] < param$nUMI
  }
  if (!ezIsSpecified(param$ngenes)) {
    if (use_mad) {
      scData$qc.nexprs <- isOutlier(
        scData@meta.data[, att_nGenes],
        nmads = param$nmad,
        log = TRUE,
        type = "lower"
      )
    } else {
      scData$qc.nexprs <- FALSE
    }
  } else {
    scData$qc.nexprs <- scData@meta.data[, att_nGenes] < param$ngenes
  }
  if (!ezIsSpecified(param$perc_mito)) {
    if (use_mad) {
      scData$qc.mito <- isOutlier(
        scData@meta.data[, "percent_mito"],
        subset = !is.na(scData@meta.data[, "percent_mito"]),
        nmads = param$nmad,
        type = "higher"
      )
    } else {
      scData$qc.mito <- FALSE
    }
  } else {
    scData$qc.mito <- scData@meta.data[, "percent_mito"] > param$perc_mito
  }
  if (!ezIsSpecified(param$perc_riboprot)) {
    if (use_mad) {
      scData$qc.riboprot <- isOutlier(
        scData@meta.data[, "percent_riboprot"],
        subset = !is.na(scData@meta.data[, "percent_riboprot"]),
        nmads = param$nmad,
        type = "higher"
      )
    } else {
      scData$qc.riboprot <- FALSE
    }
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
    doubletsInfo <- tryCatch(
      {
        scDblFinder(
          GetAssayData(scData, layer = "counts")[, scData$useCell],
          returnType = "table",
          clusters = TRUE,
          BPPARAM = BPPARAM
        )
      },
      error = function(e) {
        futile.logger::flog.warn(paste(
          "scDblFinder failed (likely too few cells), skipping doublet detection:",
          conditionMessage(e)
        ))
        NULL
      }
    )
    if (!is.null(doubletsInfo)) {
      scData$doubletScore <- doubletsInfo[colnames(scData), "score"]
      scData$doubletClass <- doubletsInfo[colnames(scData), "class"]
      scData$qc.doublet <- scData$doubletClass %in% "doublet"
      if (ezIsSpecified(param$keepDoublets) && param$keepDoublets) {
        futile.logger::flog.info("Keeping doublets...")
      } else {
        scData$useCell <- scData$useCell & scData$doubletClass %in% "singlet"
      }
    } else {
      scData$doubletScore <- NA_real_
      scData$doubletClass <- "undetermined"
      scData$qc.doublet <- FALSE
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
  # Return NULL if no databases are specified
  if (!ezIsSpecified(dbs)) {
    futile.logger::flog.info(
      "No enrichR databases specified, skipping enrichment analysis"
    )
    return(NULL)
  }
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    stop("enrichR database requested, but package 'enrichR' is not available.")
  } else {
      require(enrichR)
  }

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


##' @title The FGCZ-internal vLLM endpoint
##' @description The endpoint is treated as sensitive: it is a constant here
##'   rather than an app parameter, is kept out of SLURM logs and out of every
##'   delivered file, and is never echoed by a message. Only the model name is
##'   exposed publicly. Same convention as `AI_ENDPOINT` in app-fastQC.R.
##'   Declaring it as an appDefault would publish it, because appDefaults are
##'   written to the parameters.tsv that ships with every job.
##' @return character, the chat-completions URL.
fgczVllmEndpoint <- function() {
  Sys.getenv(
    "FGCZ_VLLM_URL",
    unset = "http://fgcz-c-056:8000/v1/chat/completions"
  )
}

##' @title Redact the vLLM host from text that may reach a user
##' @param x character vector, possibly containing the endpoint or its host.
##' @return the same text with endpoint and host replaced by placeholders.
redactVllmEndpoint <- function(x) {
  if (is.null(x) || !length(x)) {
    return(x)
  }
  endpoint <- fgczVllmEndpoint()
  hostPort <- sub("^https?://([^/]+).*", "\\1", endpoint)
  schemeHost <- sub("^(https?://[^/]+).*", "\\1", endpoint)
  for (p in list(
    c(endpoint, "[REDACTED-LLM-ENDPOINT]"),
    c(schemeHost, "[REDACTED-LLM-HOST]"),
    c(hostPort, "[REDACTED-LLM-HOST]")
  )) {
    x <- gsub(p[1], p[2], x, fixed = TRUE)
  }
  x
}

##' @title Register the FGCZ-internal vLLM as an mLLMCelltype provider
##' @description mLLMCelltype's built-in providers cannot be used against the
##'   FGCZ vLLM as-is: the processor is chosen from the model-name prefix (so a
##'   served id that looks like nothing it knows fails outright), and every
##'   built-in processor hardcodes `max_tokens = 4096` with no way to disable a
##'   reasoning trace. DeepSeek-V4 spends that whole budget thinking and returns
##'   `content = null` for roughly one request in three at ~15 clusters, which
##'   surfaces as "Unexpected response format". Registering our own provider
##'   fixes both: it turns thinking off (vLLM ignores chat-template kwargs a
##'   model does not understand, so this is safe across models), raises the token
##'   ceiling, and retries.
##' @param modelId character, the model id served by the endpoint.
##' @return Invisibly, the normalised model id.
registerFgczVllmProvider <- function(modelId) {
  processFn <- function(prompt, model, api_key, base_url = NULL) {
    stopifnot(!is.null(base_url))
    body <- list(
      model = model,
      messages = list(list(role = "user", content = prompt)),
      temperature = 0, # cluster labels should not change between reruns
      seed = 42L, # temperature 0 alone is not reproducible under vLLM batching
      max_tokens = 8192L,
      stream = FALSE,
      chat_template_kwargs = list(thinking = FALSE)
    )
    for (attempt in seq_len(3L)) {
      content <- tryCatch(
        {
          resp <- httr::POST(
            base_url,
            httr::add_headers(Authorization = paste("Bearer", api_key)),
            body = body,
            encode = "json",
            httr::timeout(600)
          )
          httr::stop_for_status(resp)
          httr::content(resp, "parsed")$choices[[1]]$message$content
        },
        error = function(e) {
          ## httr errors embed the request URL, so redact before logging
          futile.logger::flog.warn(
            "FGCZ vLLM attempt %d/3 failed: %s",
            attempt,
            redactVllmEndpoint(e$message)
          )
          NULL
        }
      )
      if (!is.null(content) && nzchar(trimws(content))) {
        return(content)
      }
      Sys.sleep(2 * attempt)
    }
    stop("FGCZ vLLM returned no usable content after 3 attempts")
  }

  ## mLLMCelltype logs to "<cwd>/logs", and the app's cwd is the report
  ## directory that ships to gstore. get_logger() takes no arguments and its
  ## constructor CREATES the directory, so simply reassigning $log_dir
  ## afterwards is too late - the empty logs/ already exists in the delivery.
  ## Force that first construction to happen in the tempdir instead: nothing is
  ## ever written where it should not be, rather than written and cleaned up.
  ## (`f()$field <- v` is also not an assignable target in R, hence the local.)
  oldwd <- setwd(tempdir())
  mllmLogger <- tryCatch(mLLMCelltype::get_logger(), finally = setwd(oldwd))
  mllmLogger$log_dir <- file.path(tempdir(), "mLLMCelltype_logs")

  if (!"fgczvllm" %in% mLLMCelltype::list_custom_providers()) {
    mLLMCelltype::register_custom_provider(
      provider_name = "fgczvllm",
      process_fn = processFn,
      description = "FGCZ-internal vLLM (OpenAI-compatible chat completions)"
    )
  }
  if (!tolower(modelId) %in% mLLMCelltype::list_custom_models()) {
    mLLMCelltype::register_custom_model(
      model_name = modelId,
      provider_name = "fgczvllm"
    )
  }
  invisible(tolower(modelId))
}

##' @title Annotate clusters with mLLMCelltype on the FGCZ-internal vLLM
##' @description Sends the top marker gene names of each cluster to an
##'   OpenAI-compatible vLLM endpoint running inside FGCZ and returns one
##'   suggested cell type per cluster. Only gene names are transmitted, and the
##'   endpoint is on-site, so no expression data leaves FGCZ.
##' @param topMarkers data.frame with `cluster` and `gene` columns, already
##'   reduced to the top N markers per cluster.
##' @param tissueName character, biological context passed to the model
##'   (e.g. "human PBMC").
##' @param vllmUrl character, full OpenAI-compatible chat-completions URL.
##' @return list with `clusterLabels` (named character vector, cluster id ->
##'   label), `model` (the served model id) and `tissueName`.
annotateClustersWithMLLMCelltype <- function(topMarkers, tissueName, vllmUrl) {
  if (!requireNamespace("mLLMCelltype", quietly = TRUE)) {
    stop(
      "mLLMCelltype package required. Install with: ",
      "pak::pkg_install('cafferychen777/mLLMCelltype/R')"
    )
  }

  ## The served model id drifts (Qwen -> DeepSeek -> ...), so read it off the
  ## endpoint rather than hardcoding one.
  modelId <- jsonlite::fromJSON(
    sub("/chat/completions/*$", "/models", vllmUrl)
  )$data$id[1]
  ## model name only - the endpoint must not reach the SLURM log
  futile.logger::flog.info("mLLMCelltype using model %s", modelId)
  registerFgczVllmProvider(modelId)

  ## Pass the markers as a NAMED LIST rather than the raw data.frame: given a
  ## data.frame, mLLMCelltype re-sorts purely numeric cluster ids internally,
  ## which would silently desync labels from clusters. Non-numeric names are
  ## kept verbatim, so the returned vector is in exactly this order.
  geneLists <- lapply(
    split(as.character(topMarkers$gene), topMarkers$cluster),
    function(genes) list(genes = genes)
  )
  geneLists <- geneLists[lengths(lapply(geneLists, `[[`, "genes")) > 0]
  stopifnot(length(geneLists) > 0)
  names(geneLists) <- paste0("cluster", names(geneLists))

  ## Ask for the reasoning form, which returns a list KEYED BY CLUSTER ID rather
  ## than one label per line. The plain form is matched POSITIONALLY, so a reply
  ## with one line too few or too many silently shifts every label by a cluster;
  ## on a real 14-cluster run it returned 13, 15 and 13 lines on three
  ## consecutive attempts, so that is the common case, not an edge case. The
  ## keyed form also carries the marker genes the model cited per cluster, which
  ## is what makes a label checkable in the report.
  labels <- NULL
  support <- NULL
  for (attempt in seq_len(3L)) {
    res <- mLLMCelltype::annotate_cell_types(
      input = geneLists,
      tissue_name = tissueName,
      model = modelId,
      ## The internal endpoint takes no auth and vLLM ignores the value, but
      ## read it from the environment anyway so pointing at an authenticated
      ## endpoint later needs no code change.
      api_key = Sys.getenv("FGCZ_VLLM_API_KEY", unset = "no-auth-required"),
      base_urls = vllmUrl,
      return_reasoning = TRUE
    )
    cand <- vapply(res, function(x) as.character(x$cell_type)[1], character(1))
    ## on a parse failure the package returns "Unknown" for every cluster
    if (all(names(geneLists) %in% names(cand)) && !all(cand == "Unknown")) {
      labels <- cand[names(geneLists)]
      support <- vapply(
        res[names(geneLists)],
        function(x) paste(as.character(x$marker_genes), collapse = ", "),
        character(1)
      )
      break
    }
    futile.logger::flog.warn(
      "mLLMCelltype returned %d usable labels for %d clusters (attempt %d/3)",
      sum(cand != "Unknown"),
      length(geneLists),
      attempt
    )
  }
  if (is.null(labels)) {
    stop("mLLMCelltype never returned a usable label for every cluster")
  }

  clusterIds <- sub("^cluster", "", names(geneLists))
  list(
    clusterLabels = setNames(as.character(labels), clusterIds),
    ## the marker genes the model cited, so a reader can check the call
    clusterSupport = setNames(as.character(support), clusterIds),
    model = modelId,
    tissueName = tissueName
  )
}

computeTFActivityAnalysis <- function(cells, species) {
    species <- tolower(species)
    # Retrieve prior knowledge network.
    if (species == 'mouse') {
        data("dorothea_mm", package = "dorothea")
        network <- dorothea_mm |>
            dplyr::filter(confidence %in% c("A", "B", "C")) |>
            dplyr::rename(source = tf) |>
            dplyr::select(source, target, mor)
    } else if (species == 'human') {
        data("dorothea_hs", package = "dorothea")
        network <- dorothea_hs |>
            dplyr::filter(confidence %in% c("A", "B", "C")) |>
            dplyr::rename(source = tf) |>
            dplyr::select(source, target, mor)
    }
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
  if(species == 'human'){
  network <- decoupleR::get_progeny(organism = species)
  } else if(species == 'mouse'){
      network <- progeny::getModel(organism = "Mouse") |>
          as.data.frame() |>
          tibble::rownames_to_column("target") |>
          tidyr::pivot_longer(-target, names_to = "source", values_to = "weight") |>
          dplyr::filter(weight != 0) |>
          dplyr::select(source, target, weight)
    }
  
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
