###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppVisiumHDSeurat <-
  setRefClass(
    "EzAppVisiumHDSeurat",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodVisiumHDSeurat
        name <<- "EzAppSeuratVisiumHD"

        # Populate RCTD References
        rctd_refs <- tryCatch(
          {
            if (dir.exists("/srv/GT/databases/RCTD_References")) {
              list.files(
                "/srv/GT/databases/RCTD_References",
                pattern = ".rds$",
                full.names = FALSE
              )
            } else {
              c("Reference_Not_Found_Locally")
            }
          },
          error = function(e) {
            c("Error_Listing_References")
          }
        )

        if (length(rctd_refs) == 0) {
          rctd_refs <- c("None")
        }

        appDefaults <<- rbind(
          nfeatures = ezFrame(
            Type = "numeric",
            DefaultValue = 3000,
            Description = "number of variable genes for SCT"
          ),
          npcs = ezFrame(
            Type = "numeric",
            DefaultValue = 50,
            Description = "The maximal dimensions to use for reduction"
          ),
          pcGenes = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "The genes used in supvervised clustering"
          ),
          SCT.regress.CellCycle = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Choose CellCycle to be regressed out when using the SCTransform method if it is a bias."
          ),
          enrichrDatabase = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "enrichR databases to search"
          ),
          DE.method = ezFrame(
            Type = "charVector",
            DefaultValue = "wilcox",
            Description = "Method to be used when calculating gene cluster markers. Use LR if you want to include cell cycle in the regression model."
          ),
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
          clusterResolution = ezFrame(
            Type = "numeric",
            DefaultValue = 0.6,
            Description = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
          ),
          cellsFraction = ezFrame(
            Type = "numeric",
            DefaultValue = 0,
            Description = "A gene will be kept if it is expressed in at least this percentage of cells"
          ),
          nUMIs = ezFrame(
            Type = "numeric",
            DefaultValue = 1,
            Description = 'A gene will be kept if it has at least nUMIs in the fraction of cells specified before'
          ),
          nmad = ezFrame(
            Type = "numeric",
            DefaultValue = 3,
            Description = "Median absolute deviation (MAD) from the median value of each metric across all cells"
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
            DefaultValue = "vst",
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
          nicheResolution = ezFrame(
            Type = "numeric",
            DefaultValue = 0.5,
            Description = "Value of the Niche resolution parameter for BANKSY clustering, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."
          ),
          rctdReference = ezFrame(
            ## a single dropdown choice: "charVector" would split the value at
            ## the commas of its "(description)" suffix
            Type = "character",
            DefaultValue = rctd_refs[1],
            Description = "RCTD Reference to use"
          ),
          rctdFile = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Manual override: Full path to custom RCTD reference .rds file"
          ),
          rctdUMImin = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "Minimum UMI count for RCTD annotation"
          ),
          rctdEngine = ezFrame(
            Type = "character",
            DefaultValue = "rctd-py",
            Description = "RCTD implementation: rctd-py (GPU if the job has one, else CPU) or spacexr (R)"
          ),
          computeCellChat = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Distance-constrained CellChat v2 between BANKSY niches and between RCTD cell types"
          ),
          cellchatInteractionLength = ezFrame(
            Type = "numeric",
            DefaultValue = 250,
            Description = "Maximal ligand-receptor interaction distance in microns"
          ),
          cellchatNboot = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "Permutations for the CellChat and SpatialCellChat p-values"
          ),
          cellchatMaxPerGroup = ezFrame(
            Type = "numeric",
            DefaultValue = 1500,
            Description = "CellChat and liana use at most this many bins per group (stratified subsample)"
          ),
          computeMisty = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "MISTy on pseudo-spots: which cell types predict each other nearby (needs RCTD)"
          ),
          nicheGridSize = ezFrame(
            Type = "numeric",
            DefaultValue = 48,
            Description = "Side length in microns of the pseudo-spots MISTy runs on"
          ),
          mistyParaviewRange = ezFrame(
            Type = "numeric",
            DefaultValue = 200,
            Description = "Length scale in microns of the MISTy neighbourhood (paraview)"
          ),
          mistyNeighbours = ezFrame(
            Type = "numeric",
            DefaultValue = 100,
            Description = "Pseudo-spots per paraview neighbourhood (nearest neighbours; avoids an N x N distance matrix)"
          ),
          computeSpatialCellChat = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "SpatialCellChat on pseudo-spots, averaged to RCTD cell types (needs RCTD; ~18 min on 8 cores for a 16 um section with 100 um pseudo-spots)"
          ),
          spatialCellChatGridSize = ezFrame(
            Type = "numeric",
            DefaultValue = 100,
            Description = "Side length in microns of the SpatialCellChat pseudo-spots"
          ),
          computeLiana = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "liana consensus ligand-receptor ranking between RCTD cell types (cross-check of CellChat)"
          )
        )
      }
    )
  )

ezMethodVisiumHDSeurat <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Visium HD Seurat")))
  on.exit(setwd(cwd), add = TRUE)
  library(Banksy)
  library(Seurat)
  library(SeuratWrappers)
  library(scater)
  library(enrichR)
  library(spacexr)
  library(future)
  library(BiocParallel)
  library(sf)

  if (param$cores > 1) {
    BPPARAM <- MulticoreParam(workers = param$cores)
  } else {
    BPPARAM <- SerialParam()
  }
  ## Pin BLAS/OpenMP to one thread before forking (MulticoreParam/future) to
  ## avoid the fork-in-multithreaded-process deadlock (e.g. AUCell labeling).
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  register(BPPARAM)
  plan("multicore", workers = param$cores)
  set.seed(38)
  future.seed = TRUE
  options(future.rng.onMisuse = "ignore")
  options(future.globals.maxSize = param$ram * 1024^3)

  ## Phase timings: one line per step so a slow run shows where the time went.
  logStep <- function(step) {
    futile.logger::flog.info("VisiumHD: %s (%d bins)", step, ncol(scData))
  }

  # Handle segmented vs binned outputs
  if (grepl("segmented", param$binSize, ignore.case = TRUE)) {
    # Segmented outputs: use parent directory and bin.size = "polygons"
    dataDir <- input$getFullPaths("SpaceRangerDir")
    scData <- Load10X_Spatial(
      data.dir = dataDir,
      image.name = "tissue_hires_image.png",
      bin.size = "polygons",
      use.names = FALSE
    )
    sf <- scData@images[[1]]@scale.factors
    scData@images[[1]]@scale.factors$lowres <- sf$hires
    matrixPath <- file.path(
      dataDir,
      param$binSize,
      "filtered_feature_cell_matrix"
    )
  } else {
    # Binned outputs: use standard approach with specific bin directory
    dataDir <- file.path(input$getFullPaths("SpaceRangerDir"), param$binSize)
    scData <- Load10X_Spatial(
      data.dir = dataDir,
      image.name = "tissue_hires_image.png",
      use.names = FALSE
    )
    sf <- scData@images[[1]]@scale.factors
    scData@images[[1]]@scale.factors$lowres <- sf$hires ## todo: needed??
    matrixPath <- file.path(dataDir, "filtered_feature_bc_matrix")
  }

  featInfo <- ezRead.table(
    paste0(matrixPath, "/features.tsv.gz"),
    header = FALSE,
    row.names = NULL
  )
  colnames(featInfo) <- c("gene_id", "gene_name", "type")
  stopifnot(length(rownames(scData)) == nrow(featInfo))
  featInfo$isMito = grepl("(?i)^MT-", featInfo$gene_name)
  featInfo$isRiboprot = grepl("(?i)^RPS|^RPL", featInfo$gene_name)
  rownames(featInfo) <- gsub(
    "_",
    "-",
    uniquifyFeatureNames(ID = featInfo$gene_id, names = featInfo$gene_name)
  )
  rownames(scData) <- rownames(featInfo)
  geneAnnoFile <- sub("byTranscript", "byGene", param$ezRef@refAnnotationFile)
  if (file.exists(geneAnnoFile)) {
    geneAnno <- ezRead.table(geneAnnoFile)
    if (any(geneAnno$type == "rRNA")) {
      featInfo$isRibosomal <- geneAnno[featInfo$gene_id, "type"] == "rRNA"
      featInfo$isRibosomal[is.na(featInfo$isRibosomal)] <- FALSE
    } else {
      featInfo$isRibosomal <- FALSE
    }
  } else {
    featInfo$isRibosomal <- FALSE
  }
  myAssay <- DefaultAssay(scData)
  scData[[myAssay]] <- AddMetaData(
    object = scData[[myAssay]],
    metadata = featInfo[rownames(scData), ]
  )
  scData@meta.data$Sample <- input$getNames()

  logStep("loaded")
  param$nUMI <- as.numeric(param$numis) ## needed by addCellQcToSeurat
  scData <- addCellQcToSeurat(
    scData,
    param = param,
    BPPARAM = BPPARAM
  )
  ## make image name unique
  stopifnot(length(names(scData@images)) == 1)
  names(scData@images) <- make.names(input$getNames())
  scData_unfiltered <- scData
  scData <- subset(scData_unfiltered, cells = which(scData_unfiltered$useCell)) # %>% head(n=1000))

  logStep("QC done")
  scData <- NormalizeData(scData)
  scData <- FindVariableFeatures(scData)
  scData <- ScaleData(scData)
  logStep("normalised and scaled")

  scData <- addCellCycleToSeurat(
    scData,
    param$refBuild,
    BPPARAM,
    assay = DefaultAssay(scData),
    method = "seurat"
  )
  logStep("cell cycle done")

  if (nrow(scData@meta.data) < 50000) {
    scData <- RunPCA(scData, npcs = 80)
    stopifnot(param$npcs <= 80)
    scData <- FindNeighbors(scData, dims = 1:param$npcs)
    logStep("neighbours done")
    scData <- findClustersFast(scData, resolution = param$clusterResolution)
    logStep("clustering done")
    scData <- RunUMAP(
      scData,
      reduction = "pca",
      reduction.name = "umap",
      return.model = T,
      dims = 1:param$npcs
    )
  } else {
    # we select 50,0000 cells and create a new 'sketch' assay
    scData <- SketchData(
      object = scData,
      ncells = 50000,
      method = "LeverageScore",
      sketched.assay = "sketch",
      features = VariableFeatures(scData)
    )
    # switch analysis to sketched cells
    DefaultAssay(scData) <- "sketch"

    # perform clustering workflow
    scData <- FindVariableFeatures(scData)
    scData <- ScaleData(scData)
    scData <- RunPCA(
      scData,
      assay = "sketch",
      reduction.name = "pca.sketch",
      npcs = 80
    )
    scData <- FindNeighbors(
      scData,
      assay = "sketch",
      reduction = "pca.sketch",
      dims = 1:param$npcs
    )
    logStep("sketch neighbours done")
    scData <- findClustersFast(
      scData,
      resolution = param$clusterResolution,
      graph.name = "sketch_snn",
      cluster.name = "seurat_clusters.sketched"
    )
    logStep("sketch clustering done")
    #scData$seurat_clusters.sketched <- scData$seurat_clusters
    scData <- RunUMAP(
      scData,
      reduction = "pca.sketch",
      reduction.name = "umap.sketch",
      return.model = T,
      dims = 1:param$npcs
    )
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
    scData$seurat_clusters.projected <- factor(
      scData$seurat_clusters.projected,
      levels(scData$seurat_clusters.sketched)
    )

    # switch to full dataset
    ## TODO: do we need seurat_clusters.projected at all
    Idents(scData) <- "seurat_clusters.projected"
    scData$seurat_clusters <- Idents(scData)
    DefaultAssay(scData) <- myAssay
  }

  logStep("UMAP done")

  # get markers and annotations
  vars.to.regress = NULL
  posMarkers <- FindAllMarkers(
    object = scData,
    test.use = param$DE.method,
    only.pos = TRUE,
    latent.vars = vars.to.regress,
    min.pct = param$min.pct,
    return.thresh = param$pvalue_allMarkers,
    logfc.threshold = param$logfc.threshold
  )
  ## Significant markers
  posMarkers <- posMarkers[, c(
    "gene",
    "cluster",
    "pct.1",
    "pct.2",
    "avg_log2FC",
    "p_val_adj"
  )]
  posMarkers$cluster <- makeGroupingVariableSortedFactor(posMarkers$cluster)
  diff_pct = abs(posMarkers$pct.1 - posMarkers$pct.2)
  posMarkers$diff_pct <- diff_pct
  posMarkers <- posMarkers[order(posMarkers$diff_pct, decreasing = TRUE), ] %>%
    mutate_if(is.numeric, round, digits = 30)
  posMarkers <- posMarkers[posMarkers$p_val_adj < param$pvalue_allMarkers, ]
  rownames(posMarkers) <- NULL
  writexl::write_xlsx(posMarkers, path = "posMarkers.xlsx")
  logStep("cluster markers done")

  ## BANKSY
  lambda <- ifelse(is.null(param$lambda), 0.8, as.numeric(param$lambda))
  niche_res <- ifelse(
    is.null(param$nicheResolution),
    0.5,
    as.numeric(param$nicheResolution)
  )

  myDefAssay <- DefaultAssay(scData)
  myIdents <- Idents(scData)
  myClusters <- scData$seurat_clusters
  scData <- RunBanksy(
    scData,
    lambda = lambda,
    assay = myDefAssay,
    slot = "data",
    features = "variable",
    k_geom = 30,
    verbose = FALSE
  )
  ## RunBanksy creates the new default assay"BANKSY"
  scData <- RunPCA(
    scData,
    assay = "BANKSY",
    reduction.name = "pca.banksy",
    features = rownames(scData),
    npcs = 30,
    verbose = FALSE
  )
  scData <- FindNeighbors(
    scData,
    reduction = "pca.banksy",
    dims = 1:12,
    verbose = FALSE
  )
  scData <- findClustersFast(
    scData,
    graph.name = "BANKSY_snn",
    cluster.name = "banksy_cluster",
    resolution = niche_res,
    verbose = FALSE
  )
  logStep("BANKSY clustering done")
  # Use original assay for marker identification (not BANKSY augmented features)
  DefaultAssay(scData) <- myDefAssay
  posMarkersBanksy <- FindAllMarkers(
    scData,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  posMarkersBanksy <- posMarkersBanksy[, c(
    "gene",
    "cluster",
    "pct.1",
    "pct.2",
    "avg_log2FC",
    "p_val_adj"
  )]
  posMarkersBanksy$cluster <- makeGroupingVariableSortedFactor(
    posMarkersBanksy$cluster
  )
  if (nrow(posMarkersBanksy) > 0) {
    posMarkersBanksy$diff_pct <- abs(
      posMarkersBanksy$pct.1 - posMarkersBanksy$pct.2
    )
    posMarkersBanksy <- posMarkersBanksy[
      order(posMarkersBanksy$diff_pct, decreasing = TRUE),
    ]
  }
  writexl::write_xlsx(posMarkersBanksy, "posMarkersBanksy.xlsx")
  logStep("BANKSY markers done")

  # Reset default assay and the transcriptional clusters: on Seurat 5.5.1
  # FindClusters() overwrites seurat_clusters even with a cluster.name (see
  # the XeniumSeurat fix and tests/testthat/test_xeniumSeuratClusters.R).
  DefaultAssay(scData) <- myDefAssay
  Idents(scData) <- myIdents
  scData$seurat_clusters <- myClusters
  # }, error = function(e) {
  #   ezLog("banksy failed", e)
  #   #writexl::write_xlsx(data.frame(), "posMarkersBanksy.xlsx")
  # })

  # 7. RCTD Annotation
  ref_path <- NULL
  if (!is.null(param$rctdFile) && param$rctdFile != "") {
    ref_path <- param$rctdFile
    cat(
      paste("Using manual RCTD reference:", ref_path, "\n"),
      file = "log.txt",
      append = TRUE
    )
  } else if (
    ezIsSpecified(param$rctdReference) && param$rctdReference != "None"
  ) {
    ref_relative <- sub(" \\([^)]+\\)$", "", param$rctdReference)
    ref_path <- file.path("/srv/GT/databases/RCTD_References", ref_relative)
    cat(
      paste("Using RCTD reference from dropdown:", ref_path, "\n"),
      file = "log.txt",
      append = TRUE
    )
  }

  rctdWeights <- NULL
  rctdResultsDf <- NULL
  if (!is.null(ref_path)) {
    stopifnot(file.exists(ref_path))
    cat(
      paste("Running RCTD with reference:", ref_path, "\n"),
      file = "log.txt",
      append = TRUE
    )
    ref_obj <- ezLoadRobj(ref_path)

    # Check if it is a spacexr Reference object
    if (!inherits(ref_obj, "Reference")) {
      if (is.list(ref_obj) && "reference" %in% names(ref_obj)) {
        ref_obj <- ref_obj$reference
      } else if (inherits(ref_obj, "Seurat")) {
        # Convert Seurat object to RCTD Reference
        cat(
          "Converting Seurat object to RCTD Reference...\n",
          file = "log.txt",
          append = TRUE
        )
        ref_counts <- Seurat::GetAssayData(ref_obj, layer = "counts")
        # Try common cell type annotation columns
        celltype_col <- intersect(
          c("author_cell_type", "cell_type", "celltype", "CellType"),
          colnames(ref_obj@meta.data)
        )[1]
        if (is.na(celltype_col)) {
          warning(
            "No cell type column found in Seurat reference. Skipping RCTD."
          )
          ref_obj <- NULL
        } else {
          ref_celltypes <- ref_obj@meta.data[[celltype_col]]
          names(ref_celltypes) <- colnames(ref_obj)
          ref_celltypes <- as.factor(ref_celltypes)
          ref_obj <- spacexr::Reference(ref_counts, ref_celltypes)
          cat(
            paste(
              "Created RCTD Reference with",
              length(levels(ref_celltypes)),
              "cell types\n"
            ),
            file = "log.txt",
            append = TRUE
          )
        }
      } else {
        warning(
          "Loaded object is not a valid RCTD Reference or Seurat object. Skipping RCTD."
        )
        ref_obj <- NULL
      }
    }

    if (!is.null(ref_obj)) {
      # Prepare Query (SpatialRNA object)
      counts <- GetAssayData(scData, assay = myAssay, layer = "counts")
      ## references keyed by Ensembl IDs (e.g. the Allen brain references) need
      ## the query genes renamed; symbol-keyed references are left as they are
      counts <- rename_rctd_query_genes(
        counts,
        setNames(featInfo$gene_id, rownames(featInfo)),
        rownames(ref_obj@counts)
      )
      coords <- GetTissueCoordinates(scData)

      # Ensure coords match counts columns
      if ("x" %in% colnames(coords) && "y" %in% colnames(coords)) {
        if ("cell" %in% colnames(coords)) {
          rownames(coords) <- coords$cell
        }
        coords <- coords[, c("x", "y")]
      } else {
        colnames(coords)[1:2] <- c("x", "y")
        if ("cell" %in% colnames(coords)) {
          rownames(coords) <- coords$cell
        }
        coords <- coords[, c("x", "y")]
      }

      # Match cells
      common_cells <- intersect(colnames(counts), rownames(coords))
      if (length(common_cells) == 0) {
        warning(
          "No common cells between counts and coordinates. Skipping RCTD."
        )
        ref_obj <- NULL
      } else {
        counts <- counts[, common_cells, drop = FALSE]
        coords <- coords[common_cells, , drop = FALSE]

        umi_min <- ifelse(
          is.null(param$rctdUMImin),
          20,
          as.numeric(param$rctdUMImin)
        )
        useRctdPy <- !identical(param$rctdEngine, "spacexr") && rctdPyAvailable()
        if (!identical(param$rctdEngine, "spacexr") && !useRctdPy) {
          futile.logger::flog.warn(
            "rctd-py env %s not installed: falling back to spacexr", RCTD_PY_ENV
          )
        }
        logStep(paste("RCTD start, engine", if (useRctdPy) "rctd-py" else "spacexr"))
        if (!useRctdPy) {
          query.puck <- SpatialRNA(coords, counts, Matrix::colSums(counts))
          myRCTD <- create.RCTD(
            query.puck,
            ref_obj,
            max_cores = param$cores,
            UMI_min = umi_min
          )
          myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
          norm_weights <- normalize_weights(myRCTD@results$weights)
          results_df <- myRCTD@results$results_df
        } else {
          ## rctd-py: same model, 14-48 h -> minutes on p32810/p42441-sized runs
          ## torch.compile on a GPU or for large inputs, see rctdPyCompile
          rctdRes <- runRctdPy(counts, coords, ref_obj, umi_min,
                               cores = as.integer(param$cores))
          norm_weights <- rctdRes$weights
          results_df <- rctdRes$results_df
        }

        rctdWeights <- as.matrix(norm_weights)
        rctdResultsDf <- results_df
        # Add to Seurat metadata - Primary cell type assignment; bins that
        # RCTD rejected get NA (shown as "Unassigned" in the report)
        max_type <- rctd_main_type(norm_weights, results_df)
        scData <- AddMetaData(
          scData,
          metadata = max_type,
          col.name = "RCTD_Main"
        )

        # Add normalized weights as metadata columns
        weight_df <- as.data.frame(norm_weights)
        colnames(weight_df) <- paste0("rctd.weight.", colnames(weight_df))
        common_cells <- intersect(colnames(scData), rownames(weight_df))
        if (length(common_cells) > 0) {
          for (wt_col in colnames(weight_df)) {
            wt_vals <- weight_df[common_cells, wt_col]
            names(wt_vals) <- common_cells
            scData <- AddMetaData(
              scData,
              metadata = wt_vals,
              col.name = wt_col
            )
          }
        }

        # Add doublet/singlet classification
        if (!is.null(results_df)) {
          scData <- AddMetaData(scData, metadata = results_df)
        }

        cat("RCTD annotation completed\n", file = "log.txt", append = TRUE)
        logStep("RCTD done")
      }
    }
  }

  # 8. Niches and cell-cell communication: on extracted inputs, never on
  # the object; each step can fail without ending the job
  status <- new.env()
  hd <- tryCatch(
    visiumhd_niche_inputs(scData, myAssay, rctdWeights, rctdResultsDf, param, input),
    error = function(e) {
      ezLog("Niche inputs could not be prepared: ", conditionMessage(e))
      NULL
    }
  )
  communication <- visiumhd_step(
    "communication", param$computeCellChat, hd, visiumhd_communication(hd, param), status,
    needsRctd = FALSE
  )
  logStep("communication done")
  misty <- visiumhd_step("misty", param$computeMisty, hd, visiumhd_misty(hd, param), status)
  logStep("misty done")
  spatialCellChat <- visiumhd_step(
    "spatialCellChat", param$computeSpatialCellChat, hd, visiumhd_spatialcellchat(hd, param), status
  )
  logStep("spatialCellChat done")
  liana <- visiumhd_step("liana", param$computeLiana, hd, visiumhd_liana(hd, param), status)
  logStep("liana done")
  nicheParams <- list(micronsPerPixel = hd$micronsPerPixel, binUm = hd$binUm, nicheGridSize = param$nicheGridSize)
  remove(hd)
  gc()
  if (!is.null(communication)) {
    qs2::qs_save(lapply(communication, `[[`, "cellchat"), "cellchat.qs2", nthreads = 4)
    sheets <- list()
    for (level in names(communication)) {
      communication[[level]]$cellchat <- NULL
      sheets[[paste0("interactions_", level)]] <- communication[[level]]$interactions
      sheets[[paste0("pathways_", level)]] <- communication[[level]]$pathways
    }
    writexl::write_xlsx(sheets, path = "cellchat.xlsx")
  }
  if (!is.null(misty)) {
    writexl::write_xlsx(misty[c("improvements", "contributions", "importances")], path = "misty.xlsx")
  }
  if (!is.null(spatialCellChat)) {
    writexl::write_xlsx(spatialCellChat[c("interactions", "pathways")], path = "spatialCellChat.xlsx")
  }
  if (!is.null(liana)) {
    writexl::write_xlsx(list(liana = liana), path = "liana.xlsx")
  }
  qs2::qs_save(
    list(
      status = as.list(status),
      params = nicheParams,
      communication = communication,
      misty = misty,
      spatialCellChat = spatialCellChat,
      liana = liana
    ),
    "nicheResults.qs2",
    nthreads = 4
  )

  ## generate template for manual cluster annotation -----
  ## we only deal with one sample
  stopifnot(length(input$getNames()) == 1)
  clusterInfos <- ezFrame(
    Sample = input$getNames(),
    Cluster = levels(Idents(scData)),
    ClusterLabel = ""
  )
  nTopMarkers <- 10
  topMarkers <- posMarkers %>%
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

  #Save some results in external files
  # bulkSignalPerCluster <- AggregateExpression(scData, group.by = 'ident', assays=myAssay)[[1]]
  # bulkSignalPerCluster <- data.frame(GeneSymbol = rownames(scData), Count = bulkSignalPerCluster)
  bulkSignalPerSample <- AggregateExpression(
    scData,
    group.by = 'Sample',
    assays = myAssay
  )[[1]]
  bulkSignalPerSample <- data.frame(
    GeneSymbol = rownames(scData),
    Count = bulkSignalPerSample
  )
  writexl::write_xlsx(bulkSignalPerSample, path = "bulkSignalPerSample.xlsx")

  logStep("rendering report")
  makeRmdReport(
    param = param,
    output = output,
    input = input,
    scData = scData,
    scData_unfiltered = scData_unfiltered,
    rmdFile = "VisiumHDSeurat.Rmd",
    reportTitle = paste0(param$name, ": ", input$getNames()),
    use.qs2 = TRUE
  )

  logStep("report done")
  gc()
  return("Success")
}
