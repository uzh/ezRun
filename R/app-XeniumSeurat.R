# App wrapper for Xenium Seurat Analysis
#
# This app processes Xenium data using Seurat, performs QC, Clustering, and RCTD Annotation.
# It follows FGCZ best practices and recommendations.
# Runs in SAMPLE mode - one job per sample.

#' Patched RunBanksy function for Seurat v5 compatibility
#' Uses 'layer' instead of deprecated 'slot' argument (SeuratObject >= 5.0.0)
#' Issue: https://github.com/satijalab/seurat/issues/10250
RunBanksy_v5 <- function(
  object,
  lambda,
  assay = "RNA",
  layer = "data",
  use_agf = FALSE,
  dimx = NULL,
  dimy = NULL,
  dimz = NULL,
  ndim = 2,
  features = "variable",
  group = NULL,
  split.scale = TRUE,
  k_geom = 15,
  n = 2,
  sigma = 1.5,
  alpha = 0.05,
  k_spatial = 10,
  spatial_mode = "kNN_median",
  assay_name = "BANKSY",
  M = NULL,
  chunk_size = NULL,
  parallel = FALSE,
  num_cores = NULL,
  verbose = TRUE
) {
  if (lambda < 0 || lambda > 1) {
    stop("Lambda must be between 0 and 1")
  }

  # Get data using 'layer' instead of deprecated 'slot'
  if (verbose) {
    message("Fetching data from layer ", layer, " from assay ", assay)
  }
  data_own <- Seurat::GetAssayData(
    object = object,
    assay = assay,
    layer = layer
  )

  if (features[1] != "all") {
    if (verbose) {
      message("Subsetting by features")
    }
    if (features[1] == "variable") {
      feat <- Seurat::VariableFeatures(object)
      if (length(feat) == 0) {
        warning("No variable features found. Running FindVariableFeatures")
        object <- Seurat::FindVariableFeatures(object)
        feat <- Seurat::VariableFeatures(object)
      }
    } else {
      feat <- features[which(rownames(object) %in% features)]
      if (length(feat) == 0) stop("None of the specified features found")
    }
    data_own <- data_own[feat, , drop = FALSE]
  }
  data_own <- as.matrix(x = data_own)

  # Get locations
  if (!is.null(dimx) & !is.null(dimy)) {
    locs <- data.frame(
      sdimx = unlist(object[[dimx]]),
      sdimy = unlist(object[[dimy]])
    )
    rownames(locs) <- colnames(object)
    if (!is.null(dimz)) {
      locs$sdimz <- object[[dimz]]
    }
    obj_samples <- colnames(data_own)
    locs_samples <- rownames(locs)
    if (any(is.na(match(obj_samples, locs_samples)))) {
      na_id <- which(is.na(match(obj_samples, locs_samples)))
      warning(
        "No centroids found for samples: ",
        paste(obj_samples[na_id], collapse = ", ")
      )
      data_own <- data_own[, -na_id, drop = FALSE]
    }
    locs <- locs[match(obj_samples, locs_samples), , drop = FALSE]
  } else {
    coords <- Seurat::GetTissueCoordinates(object)
    if ("cell" %in% colnames(coords)) {
      rownames(coords) <- coords$cell
    }
    locs <- coords[, seq_len(ndim), drop = FALSE]
  }

  dim_names <- paste0("sdim", c("x", "y", "z"))
  colnames(locs) <- dim_names[seq_len(ncol(locs))]

  if (!is.null(group)) {
    if (verbose) {
      message("Staggering locations by ", group)
    }
    locs[, 1] <- locs[, 1] + abs(min(locs[, 1]))
    max_x <- max(locs[, 1]) * 2
    n_groups <- length(unique(unlist(object[[group]])))
    shift <- seq(from = 0, length.out = n_groups, by = max_x)
    shift_order <- match(
      unique(unlist(object[[group]])),
      names(table(object[[group]]))
    )
    locs[, 1] <- locs[, 1] + rep(shift, table(object[[group]])[shift_order])
    object <- Seurat::AddMetaData(
      object,
      metadata = locs,
      col.name = paste0("staggered_", colnames(locs))
    )
  }

  # Match cells
  common_cells <- intersect(colnames(data_own), rownames(locs))
  if (length(common_cells) == 0) {
    stop("No common cells between expression data and coordinates")
  }
  if (length(common_cells) < ncol(data_own)) {
    if (verbose) {
      message("Subsetting to ", length(common_cells), " common cells")
    }
    data_own <- data_own[, common_cells, drop = FALSE]
    locs <- locs[common_cells, , drop = FALSE]
  }

  # Compute BANKSY
  knn_list <- lapply(k_geom, function(kg) {
    Banksy:::computeNeighbors(
      locs,
      spatial_mode = spatial_mode,
      k_geom = kg,
      n = n,
      sigma = sigma,
      alpha = alpha,
      k_spatial = k_spatial,
      verbose = verbose
    )
  })

  M_seq <- seq(0, max(Banksy:::getM(use_agf, M)))
  center <- rep(TRUE, length(M_seq))
  center[1] <- FALSE

  har <- Map(
    function(knn_df, M_val, cen) {
      x <- Banksy:::computeHarmonics(
        gcm = data_own,
        knn_df = knn_df,
        M = M_val,
        center = cen,
        verbose = verbose,
        chunk_size = chunk_size,
        parallel = parallel,
        num_cores = num_cores
      )
      rownames(x) <- paste0(rownames(x), ".m", M_val)
      x
    },
    knn_list,
    M_seq,
    center
  )

  lambdas <- Banksy:::getLambdas(lambda, n_harmonics = length(har))
  if (verbose) {
    message("Creating Banksy matrix")
  }
  data_banksy <- c(list(data_own), har)

  if (verbose) {
    message(
      "Scaling BANKSY matrix. Do not call ScaleData on assay ",
      assay_name
    )
  }
  fast_scaler <- function(mat, obj, grp, split.scl, verb) {
    if (!is.null(grp) && split.scl) {
      groups <- unique(unlist(obj[[grp]]))
      scaled_list <- lapply(groups, function(g) {
        idx <- which(unlist(obj[[grp]]) == g)
        idx <- idx[idx %in% seq_len(ncol(mat))]
        if (length(idx) > 0) t(scale(t(mat[, idx, drop = FALSE]))) else NULL
      })
      scaled_list <- scaled_list[!sapply(scaled_list, is.null)]
      if (length(scaled_list) > 0) {
        do.call(cbind, scaled_list)[, colnames(mat)]
      } else {
        t(scale(t(mat)))
      }
    } else {
      t(scale(t(mat)))
    }
  }

  data_scaled <- lapply(
    data_banksy,
    fast_scaler,
    object,
    group,
    split.scale,
    verbose
  )
  data_banksy <- Map(function(lam, mat) lam * mat, lambdas, data_banksy)
  data_scaled <- Map(function(lam, mat) lam * mat, lambdas, data_scaled)
  data_banksy <- do.call(rbind, data_banksy)
  data_scaled <- do.call(rbind, data_scaled)

  # Create BANKSY assay
  if (grepl(pattern = "counts", x = layer)) {
    banksy_assay <- Seurat::CreateAssayObject(counts = data_banksy)
  } else {
    banksy_assay <- Seurat::CreateAssayObject(data = data_banksy)
  }

  if (verbose) {
    message("Setting default assay to ", assay_name)
  }
  object[[assay_name]] <- banksy_assay
  Seurat::DefaultAssay(object) <- assay_name
  object <- Seurat::SetAssayData(
    object,
    layer = "scale.data",
    new.data = data_scaled,
    assay = assay_name
  )
  if (verbose) {
    message("BANKSY assay created successfully")
  }
  return(object)
}

ezMethodXeniumSeurat <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
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

  ezWrite(
    paste("Processing sample:", sampleName, "at", xeniumPath),
    "log.txt",
    append = TRUE
  )

  if (!dir.exists(xeniumPath)) {
    stop(paste("Error: Path not found:", xeniumPath))
  }

  # 1. Load Xenium Data
  sdata <- LoadXenium(data.dir = xeniumPath, fov = "fov")

  # Load additional cell metadata not loaded by Seurat (cell_area, nucleus_area)
  cells_file <- file.path(xeniumPath, "cells.csv.gz")
  if (file.exists(cells_file)) {
    cells_meta <- data.table::fread(
      cells_file,
      select = c("cell_id", "cell_area", "nucleus_area")
    )
    meta_idx <- match(colnames(sdata), cells_meta$cell_id)

    if (sum(!is.na(meta_idx)) > 0) {
      sdata$cell_area <- cells_meta$cell_area[meta_idx]
      sdata$nucleus_area <- cells_meta$nucleus_area[meta_idx]
      ezWrite(
        paste(
          "Added cell_area and nucleus_area to",
          sum(!is.na(meta_idx)),
          "cells"
        ),
        "log.txt",
        append = TRUE
      )
    }
  }

  # Rename cells with sample prefix
  sdata <- RenameCells(sdata, add.cell.id = sampleName)
  sdata$orig.ident <- sampleName
  sdata$Sample <- sampleName

  # 2. QC - compute thresholds and flags
  min_counts <- ifelse(
    is.null(param$minCounts),
    10,
    as.numeric(param$minCounts)
  )
  min_features <- ifelse(
    is.null(param$minFeatures),
    5,
    as.numeric(param$minFeatures)
  )

  # Add QC flags (matching VisiumHD style for consistent reporting)
  sdata$qc.lib <- sdata$nCount_Xenium < min_counts
  sdata$qc.nexprs <- sdata$nFeature_Xenium < min_features
  sdata$useCell <- !(sdata$qc.lib | sdata$qc.nexprs)
  sdata$discard <- !sdata$useCell

  # Store unfiltered object for QC reporting (with QC flags)
  ezWrite(
    "Saving unfiltered object for QC reporting...",
    "log.txt",
    append = TRUE
  )
  qs2::qs_save(sdata, "scData.unfiltered.qs2", nthreads = 8)

  ezWrite(
    paste(
      "Filtering cells: minCounts =",
      min_counts,
      ", minFeatures =",
      min_features
    ),
    "log.txt",
    append = TRUE
  )
  cells_before <- ncol(sdata)
  sdata <- subset(sdata, cells = which(sdata$useCell))
  cells_after <- ncol(sdata)
  ezWrite(
    paste("Cells after QC:", cells_after, "/", cells_before),
    "log.txt",
    append = TRUE
  )

  # 3. Normalization (Xenium-specific with median counts as scale factor)
  sdata <- NormalizeData(sdata, scale.factor = median(sdata$nCount_Xenium))

  # 4. Feature Selection & Scaling
  sdata <- FindVariableFeatures(
    sdata,
    selection.method = "vst",
    nfeatures = 2000
  )
  sdata <- ScaleData(sdata)

  # 5. Dimension Reduction
  sdata <- RunPCA(sdata, verbose = FALSE)
  sdata <- RunUMAP(sdata, dims = 1:30, verbose = FALSE)

  # 6. Clustering
  res <- ifelse(
    is.null(param$clusterResolution),
    0.5,
    as.numeric(param$clusterResolution)
  )
  sdata <- FindNeighbors(sdata, dims = 1:30, verbose = FALSE)
  sdata <- FindClusters(sdata, resolution = res, verbose = FALSE)

  # 7. RCTD Annotation
  ref_path <- NULL
  if (!is.null(param$rctdFile) && param$rctdFile != "") {
    ref_path <- param$rctdFile
    ezWrite(
      paste("Using manual RCTD reference:", ref_path),
      "log.txt",
      append = TRUE
    )
  } else if (
    ezIsSpecified(param$rctdReference) && param$rctdReference != "None"
  ) {
    ref_relative <- sub(" \\([^)]+\\)$", "", param$rctdReference)
    ref_path <- file.path("/srv/GT/databases/RCTD_References", ref_relative)
    ezWrite(
      paste("Using RCTD reference from dropdown:", ref_path),
      "log.txt",
      append = TRUE
    )
  }

  if (!is.null(ref_path)) {
    stopifnot(file.exists(ref_path))
    ezWrite(
      paste("Running RCTD with reference:", ref_path),
      "log.txt",
      append = TRUE
    )
    ref_obj <- ezLoadRobj(ref_path)

    # Check if it is a spacexr Reference object
    if (!inherits(ref_obj, "Reference")) {
      if (is.list(ref_obj) && "reference" %in% names(ref_obj)) {
        ref_obj <- ref_obj$reference
      } else if (inherits(ref_obj, "Seurat")) {
        # Convert Seurat object to RCTD Reference
        ezWrite(
          "Converting Seurat object to RCTD Reference...",
          "log.txt",
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
          # RCTD requires cell_types to be a factor
          ref_celltypes <- as.factor(ref_celltypes)
          ref_obj <- spacexr::Reference(ref_counts, ref_celltypes)
          ezWrite(
            paste(
              "Created RCTD Reference with",
              length(levels(ref_celltypes)),
              "cell types"
            ),
            "log.txt",
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
      counts <- GetAssayData(sdata, assay = "Xenium", layer = "counts")
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
      if (length(common_cells) == 0) {
        warning(
          "No common cells between counts and coordinates. Skipping RCTD."
        )
        ref_obj <- NULL
      } else {
        counts <- counts[, common_cells, drop = FALSE]
        coords <- coords[common_cells, , drop = FALSE]

        # Create SpatialRNA object
        query.puck <- SpatialRNA(coords, counts, Matrix::colSums(counts))

        # Run RCTD
        umi_min <- ifelse(
          is.null(param$rctdUMImin),
          20,
          as.numeric(param$rctdUMImin)
        )
        myRCTD <- create.RCTD(
          query.puck,
          ref_obj,
          max_cores = param$cores,
          UMI_min = umi_min
        )
        myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

        # Extract results
        results <- myRCTD@results
        norm_weights <- normalize_weights(results$weights)

        # Add to Seurat metadata - Primary cell type assignment
        max_type <- colnames(norm_weights)[max.col(
          norm_weights,
          ties.method = "first"
        )]
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
    }
  } else if (!is.null(ref_path)) {
    ezWrite(
      paste("RCTD reference not found:", ref_path),
      "log.txt",
      append = TRUE
    )
  }

  # Assign final object
  scData <- sdata

  # Pre-compute cluster markers
  ezWrite("Computing cluster markers...", "log.txt", append = TRUE)
  Idents(scData) <- "seurat_clusters"
  posMarkers <- FindAllMarkers(
    scData,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  if (nrow(posMarkers) > 0) {
    posMarkers$diff_pct <- abs(posMarkers$pct.1 - posMarkers$pct.2)
    posMarkers <- posMarkers[order(posMarkers$diff_pct, decreasing = TRUE), ]
  }
  writexl::write_xlsx(posMarkers, "posMarkers.xlsx")

  # BANKSY Spatial Niches
  ezWrite("Running BANKSY spatial analysis...", "log.txt", append = TRUE)
  tryCatch(
    {
      lambda <- ifelse(is.null(param$lambda), 0.8, as.numeric(param$lambda))
      niche_res <- ifelse(
        is.null(param$nicheResolution),
        0.5,
        as.numeric(param$nicheResolution)
      )

      # Use embedded RunBanksy_v5 for Seurat v5 compatibility
      # SeuratWrappers::RunBanksy uses deprecated 'slot' argument that fails with SeuratObject >= 5.0
      scData <- RunBanksy_v5(
        scData,
        lambda = lambda,
        assay = "Xenium",
        layer = "data",
        features = "variable",
        k_geom = 30,
        verbose = FALSE
      )
      DefaultAssay(scData) <- "BANKSY"
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
      scData <- FindClusters(
        scData,
        cluster.name = "banksy_cluster",
        resolution = niche_res,
        verbose = FALSE
      )

      # Niche markers
      Idents(scData) <- "banksy_cluster"
      posMarkersBanksy <- FindAllMarkers(
        scData,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.25,
        verbose = FALSE
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

      # Reset default assay
      DefaultAssay(scData) <- "Xenium"
      Idents(scData) <- "seurat_clusters"
      ezWrite("BANKSY analysis completed", "log.txt", append = TRUE)
    },
    error = function(e) {
      ezWrite(
        paste("BANKSY analysis failed:", e$message),
        "log.txt",
        append = TRUE
      )
      writexl::write_xlsx(data.frame(), "posMarkersBanksy.xlsx")
    }
  )

  # Export Xenium Explorer compatible CSV files
  # Cell IDs in Seurat match original Xenium format (e.g., "aaaddlda-1")
  ezWrite("Exporting Xenium Explorer CSV files...", "log.txt", append = TRUE)

  cell_ids <- colnames(scData)

  # Clusters for Xenium Explorer
  write.csv(
    data.frame(
      cell_id = cell_ids,
      group = as.character(scData$seurat_clusters)
    ),
    "clusters_for_explorer.csv",
    row.names = FALSE
  )

  # Niches for Xenium Explorer
  if ("banksy_cluster" %in% colnames(scData@meta.data)) {
    write.csv(
      data.frame(
        cell_id = cell_ids,
        group = as.character(scData$banksy_cluster)
      ),
      "niches_for_explorer.csv",
      row.names = FALSE
    )
  }

  # Cell types for Xenium Explorer
  if ("RCTD_Main" %in% colnames(scData@meta.data)) {
    write.csv(
      data.frame(
        cell_id = cell_ids,
        group = as.character(scData$RCTD_Main)
      ),
      "celltypes_for_explorer.csv",
      row.names = FALSE
    )
  }

  # Generate Vitessce-optimized Zarr for fast visualization
  # This enables <5 second load times in exploreVitessceXenium
  # Disabled by default until fully tested
  if (
    isTRUE(param$generateVitessceZarr) && !is.null(param$generateVitessceZarr)
  ) {
    ezWrite(
      "Generating Vitessce Zarr for fast visualization...",
      con = "log.txt",
      append = TRUE
    )

    tryCatch(
      {
        # Try ezRun bundled script first, fall back to development location
        vitessce_script <- system.file(
          "R",
          "vitessce-utils.R",
          package = "ezRun"
        )
        if (!file.exists(vitessce_script) || vitessce_script == "") {
          vitessce_script <- "/home/pgueguen/git/paul-scripts/Internal_Dev/SeuratXeniumApp/scripts/generate_vitessce_zarr.R"
        }

        if (file.exists(vitessce_script)) {
          source(vitessce_script)

          result <- generate_vitessce_zarr(
            seurat_path = "scData.qs2",
            output_dir = ".",
            seurat_obj = scData,
            nthreads = 8
          )

          ezWrite(
            paste("Vitessce Zarr created:", result$zarr_path),
            con = "log.txt",
            append = TRUE
          )
        } else {
          ezWrite(
            "Warning: Vitessce generation script not found, skipping",
            con = "log.txt",
            append = TRUE
          )
        }
      },
      error = function(e) {
        ezWrite(
          paste("Warning: Vitessce Zarr generation failed:", e$message),
          con = "log.txt",
          append = TRUE
        )
      }
    )
  }

  # Save final analyzed object for Rmd report
  ezWrite("Saving final analyzed object...", "log.txt", append = TRUE)
  qs2::qs_save(scData, "scData.qs2", nthreads = param$cores)

  # Generate Report
  makeRmdReport(
    param = param,
    scData = scData,
    input = input,
    use.qs = TRUE,
    rmdFile = "XeniumSeurat.Rmd",
    reportTitle = paste0(sampleName, " - Xenium Analysis")
  )

  return("Success")
}

EzAppXeniumSeurat <- setRefClass(
  "EzAppXeniumSeurat",
  contains = "EzApp",
  methods = list(
    initialize = function() {
      "Initializes the application using its specific defaults."
      runMethod <<- ezMethodXeniumSeurat
      name <<- "EzAppXeniumSeurat"

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
        minCounts = ezFrame(
          Type = "numeric",
          DefaultValue = 10,
          Description = "Minimum counts per cell"
        ),
        minFeatures = ezFrame(
          Type = "numeric",
          DefaultValue = 5,
          Description = "Minimum features per cell"
        ),
        clusterResolution = ezFrame(
          Type = "numeric",
          DefaultValue = 0.5,
          Description = "Seurat clustering resolution"
        ),
        lambda = ezFrame(
          Type = "numeric",
          DefaultValue = 0.8,
          Description = "BANKSY spatial weighting (0-1)"
        ),
        nicheResolution = ezFrame(
          Type = "numeric",
          DefaultValue = 0.5,
          Description = "BANKSY niche clustering resolution"
        ),
        rctdReference = ezFrame(
          Type = "charVector",
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
        computeSC = ezFrame(
          Type = "logical",
          DefaultValue = TRUE,
          Description = "Compute Seurat Analysis"
        )
      )
    }
  )
)
