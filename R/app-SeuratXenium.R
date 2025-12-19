# App wrapper for Seurat Xenium Analysis
#
# This app processes Xenium data using Seurat, performs QC, Clustering, and RCTD Annotation.
# It follows FGCZ best practices and recommendations.
# Runs in SAMPLE mode - one job per sample.

ezMethodSeuratXenium <- function(input = NA, output = NA, param = NA, htmlFile = "00index.html") {
  ezLoadPackage("Seurat")
  ezLoadPackage("spacexr")
  ezLoadPackage("ggplot2")
  ezLoadPackage("dplyr")
  ezLoadPackage("scales")
  ezLoadPackage("patchwork")
  ezLoadPackage("data.table")
  ezLoadPackage("qs2")
  ezLoadPackage("writexl")
  ezLoadPackage("SeuratWrappers")
  ezLoadPackage("Banksy")

  # Setup directory - use sample name from input
  cwd <- getwd()
  sampleName <- input$getNames()
  setwdNew(sampleName)

  # Get input path for this sample
  xeniumPath <- file.path(param$dataRoot, input$getColumn("XeniumPath"))

  ezWrite(paste("Processing sample:", sampleName, "at", xeniumPath), "log.txt", append = TRUE)

  if (!dir.exists(xeniumPath)) {
    stop(paste("Error: Path not found:", xeniumPath))
  }

  # 1. Load Xenium Data
  sdata <- LoadXenium(data.dir = xeniumPath, fov = "fov")

  # Load additional cell metadata not loaded by Seurat (cell_area, nucleus_area)
  cells_file <- file.path(xeniumPath, "cells.csv.gz")
  if (file.exists(cells_file)) {
    cells_meta <- data.table::fread(cells_file, select = c("cell_id", "cell_area", "nucleus_area"))
    meta_idx <- match(colnames(sdata), cells_meta$cell_id)

    if (sum(!is.na(meta_idx)) > 0) {
      sdata$cell_area <- cells_meta$cell_area[meta_idx]
      sdata$nucleus_area <- cells_meta$nucleus_area[meta_idx]
      ezWrite(paste("Added cell_area and nucleus_area to", sum(!is.na(meta_idx)), "cells"), "log.txt", append = TRUE)
    }
  }

  # Rename cells with sample prefix
  sdata <- RenameCells(sdata, add.cell.id = sampleName)
  sdata$orig.ident <- sampleName
  sdata$Sample <- sampleName

  # 2. QC
  min_counts <- ifelse(is.null(param$minCounts), 10, as.numeric(param$minCounts))
  min_features <- ifelse(is.null(param$minFeatures), 5, as.numeric(param$minFeatures))

  ezWrite(paste("Filtering cells: minCounts =", min_counts, ", minFeatures =", min_features), "log.txt", append = TRUE)
  cells_before <- ncol(sdata)
  sdata <- subset(sdata, subset = nCount_Xenium > min_counts & nFeature_Xenium > min_features)
  cells_after <- ncol(sdata)
  ezWrite(paste("Cells after QC:", cells_after, "/", cells_before), "log.txt", append = TRUE)

  # 3. Normalization (Xenium-specific with median counts as scale factor)
  sdata <- NormalizeData(sdata, scale.factor = median(sdata$nCount_Xenium))

  # 4. Feature Selection & Scaling
  sdata <- FindVariableFeatures(sdata, selection.method = "vst", nfeatures = 2000)
  sdata <- ScaleData(sdata)

  # 5. Dimension Reduction
  sdata <- RunPCA(sdata, verbose = FALSE)
  sdata <- RunUMAP(sdata, dims = 1:30, verbose = FALSE)

  # 6. Clustering
  res <- ifelse(is.null(param$Cluster_resolution), 0.5, as.numeric(param$Cluster_resolution))
  sdata <- FindNeighbors(sdata, dims = 1:30, verbose = FALSE)
  sdata <- FindClusters(sdata, resolution = res, verbose = FALSE)

  # 7. RCTD Annotation
  ref_path <- NULL
  if (!is.null(param$rctdFile) && param$rctdFile != "") {
    ref_path <- param$rctdFile
    ezWrite(paste("Using manual RCTD reference:", ref_path), "log.txt", append = TRUE)
  } else if (!is.null(param$rctdReference) && param$rctdReference != "" && param$rctdReference != "None") {
    ref_relative <- sub(" \\([^)]+\\)$", "", param$rctdReference)
    ref_path <- file.path("/srv/GT/databases/RCTD_References", ref_relative)
    ezWrite(paste("Using RCTD reference from dropdown:", ref_path), "log.txt", append = TRUE)
  }

  if (!is.null(ref_path) && file.exists(ref_path)) {
    ezWrite(paste("Running RCTD with reference:", ref_path), "log.txt", append = TRUE)

    # Load Reference
    ref_obj <- readRDS(ref_path)

    # Check if it is a spacexr Reference object
    if (!inherits(ref_obj, "Reference")) {
      if (is.list(ref_obj) && "reference" %in% names(ref_obj)) {
        ref_obj <- ref_obj$reference
      } else {
        warning("Loaded object is not a valid RCTD Reference. Skipping RCTD.")
        ref_obj <- NULL
      }
    }

    if (!is.null(ref_obj)) {
      # Prepare Query (SpatialRNA object)
      counts <- GetAssayData(sdata, assay = "Xenium", layer = "counts")
      # Ensure counts is a proper sparse matrix (required for RCTD)
      if (!inherits(counts, "dgCMatrix")) {
        counts <- as(counts, "dgCMatrix")
      }
      coords <- GetTissueCoordinates(sdata)

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
      counts <- counts[, common_cells, drop = FALSE]
      coords <- coords[common_cells, , drop = FALSE]

      # Create SpatialRNA object
      query.puck <- SpatialRNA(coords, counts, colSums(counts))

      # Run RCTD
      umi_min <- ifelse(is.null(param$rctdUMImin), 100, as.numeric(param$rctdUMImin))
      myRCTD <- create.RCTD(query.puck, ref_obj, max_cores = param$cores, UMI_min = umi_min)
      myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

      # Extract results
      results <- myRCTD@results
      norm_weights <- normalize_weights(results$weights)

      # Add to Seurat metadata - Primary cell type assignment
      max_type <- colnames(norm_weights)[max.col(norm_weights, ties.method = "first")]
      names(max_type) <- rownames(norm_weights)
      sdata <- AddMetaData(sdata, metadata = max_type, col.name = "RCTD_Main")

      # Add normalized weights as metadata columns
      weight_df <- as.data.frame(norm_weights)
      colnames(weight_df) <- paste0("rctd.weight.", colnames(weight_df))
      common_cells <- intersect(colnames(sdata), rownames(weight_df))
      if (length(common_cells) > 0) {
        for (wt_col in colnames(weight_df)) {
          wt_vals <- weight_df[common_cells, wt_col]
          names(wt_vals) <- common_cells
          sdata <- AddMetaData(sdata, metadata = wt_vals, col.name = wt_col)
        }
      }

      # Add doublet/singlet classification
      if (!is.null(results$results_df)) {
        sdata <- AddMetaData(sdata, metadata = results$results_df)
      }

      ezWrite("RCTD annotation completed", "log.txt", append = TRUE)
    }
  } else if (!is.null(ref_path)) {
    ezWrite(paste("RCTD reference not found:", ref_path), "log.txt", append = TRUE)
  }

  # Assign final object
  scData <- sdata

  # Pre-compute cluster markers
  ezWrite("Computing cluster markers...", "log.txt", append = TRUE)
  Idents(scData) <- "seurat_clusters"
  posMarkers <- FindAllMarkers(scData, only.pos = TRUE, min.pct = 0.25,
                                logfc.threshold = 0.25, verbose = FALSE)
  if (nrow(posMarkers) > 0) {
    posMarkers$diff_pct <- abs(posMarkers$pct.1 - posMarkers$pct.2)
    posMarkers <- posMarkers[order(posMarkers$diff_pct, decreasing = TRUE), ]
  }
  writexl::write_xlsx(posMarkers, "posMarkers.xlsx")

  # BANKSY Spatial Niches
  ezWrite("Running BANKSY spatial analysis...", "log.txt", append = TRUE)
  tryCatch({
    lambda <- ifelse(is.null(param$lambda), 0.8, as.numeric(param$lambda))
    niche_res <- ifelse(is.null(param$Niche_resolution), 0.5, as.numeric(param$Niche_resolution))

    scData <- RunBanksy(scData, lambda = lambda, assay = "Xenium",
                        slot = "data", features = "variable", k_geom = 30,
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
    DefaultAssay(scData) <- "Xenium"
    Idents(scData) <- "seurat_clusters"
    ezWrite("BANKSY analysis completed", "log.txt", append = TRUE)
  }, error = function(e) {
    ezWrite(paste("BANKSY analysis failed:", e$message), "log.txt", append = TRUE)
    writexl::write_xlsx(data.frame(), "posMarkersBanksy.xlsx")
  })

  # Export Xenium Explorer compatible CSV files
  ezWrite("Exporting Xenium Explorer CSV files...", "log.txt", append = TRUE)

  extract_cell_id <- function(cell_names) {
    sapply(cell_names, function(x) {
      parts <- strsplit(as.character(x), "-")[[1]]
      last <- parts[length(parts)]
      if (grepl("^[0-9]+$", last)) return(last)
      nums <- regmatches(x, regexpr("[0-9]+", x))
      if (length(nums) > 0) return(nums)
      return(x)
    }, USE.NAMES = FALSE)
  }

  cell_ids <- extract_cell_id(colnames(scData))

  # Clusters for Xenium Explorer
  write.csv(data.frame(
    cell_id = cell_ids,
    group = as.character(scData$seurat_clusters)
  ), "clusters_for_explorer.csv", row.names = FALSE)

  # Niches for Xenium Explorer
  if ("banksy_cluster" %in% colnames(scData@meta.data)) {
    write.csv(data.frame(
      cell_id = cell_ids,
      group = as.character(scData$banksy_cluster)
    ), "niches_for_explorer.csv", row.names = FALSE)
  }

  # Cell types for Xenium Explorer
  if ("RCTD_Main" %in% colnames(scData@meta.data)) {
    write.csv(data.frame(
      cell_id = cell_ids,
      group = as.character(scData$RCTD_Main)
    ), "celltypes_for_explorer.csv", row.names = FALSE)
  }

  # Save Results as qs2
  qs2::qs_save(scData, "scData.qs2", nthreads = 8)
  qs2::qs_save(param, "param.qs2")
  qs2::qs_save(input, "input.qs2")

  # Generate Vitessce-optimized Zarr for fast visualization
  # This enables <5 second load times in exploreVitessceXenium
  # Disabled by default until fully tested
  if (isTRUE(param$generateVitessceZarr) && !is.null(param$generateVitessceZarr)) {
    ezWrite("Generating Vitessce Zarr for fast visualization...", con = "log.txt", append = TRUE)

    tryCatch({
      vitessce_script <- "/home/pgueguen/git/paul-scripts/Internal_Dev/SeuratXeniumApp/scripts/generate_vitessce_zarr.R"

      if (file.exists(vitessce_script)) {
        source(vitessce_script)

        result <- generate_vitessce_zarr(
          seurat_path = "scData.qs2",
          output_dir = ".",
          seurat_obj = scData,
          nthreads = 8
        )

        ezWrite(paste("Vitessce Zarr created:", result$zarr_path), con = "log.txt", append = TRUE)
      } else {
        ezWrite("Warning: Vitessce generation script not found, skipping", con = "log.txt", append = TRUE)
      }
    }, error = function(e) {
      ezWrite(paste("Warning: Vitessce Zarr generation failed:", e$message), con = "log.txt", append = TRUE)
    })
  }

  # Generate Report
  makeRmdReport(rmdFile = "SeuratXenium.Rmd", reportTitle = paste0(sampleName, " - Xenium Analysis"))

  return("Success")
}

EzAppSeuratXenium <- setRefClass("EzAppSeuratXenium",
  contains = "EzApp",
  methods = list(
    initialize = function() {
      "Initializes the application using its specific defaults."
      runMethod <<- ezMethodSeuratXenium
      name <<- "EzAppSeuratXenium"

      # Populate RCTD References
      rctd_refs <- tryCatch({
        if (dir.exists("/srv/GT/databases/RCTD_References")) {
           list.files("/srv/GT/databases/RCTD_References", pattern = ".rds$", full.names = FALSE)
        } else {
           c("Reference_Not_Found_Locally")
        }
      }, error = function(e) {
        c("Error_Listing_References")
      })

      if (length(rctd_refs) == 0) rctd_refs <- c("None")

      appDefaults <<- rbind(
        minCounts = ezFrame(Type = "numeric", DefaultValue = 10,
                            Description = "Minimum counts per cell"),
        minFeatures = ezFrame(Type = "numeric", DefaultValue = 5,
                              Description = "Minimum features per cell"),
        Cluster_resolution = ezFrame(Type = "numeric", DefaultValue = 0.5,
                             Description = "Seurat clustering resolution"),
        lambda = ezFrame(Type = "numeric", DefaultValue = 0.8,
                         Description = "BANKSY spatial weighting (0-1)"),
        Niche_resolution = ezFrame(Type = "numeric", DefaultValue = 0.5,
                                   Description = "BANKSY niche clustering resolution"),
        rctdReference = ezFrame(Type = "charVector", DefaultValue = rctd_refs[1],
                                Description = "RCTD Reference to use"),
        rctdFile = ezFrame(Type = "character", DefaultValue = "",
                           Description = "Manual override: Full path to custom RCTD reference .rds file"),
        rctdUMImin = ezFrame(Type = "numeric", DefaultValue = 100,
                             Description = "Minimum UMI count for RCTD annotation"),
        computeSC = ezFrame(Type = "logical", DefaultValue = TRUE,
                            Description = "Compute Seurat Analysis")
      )
    }
  )
)
