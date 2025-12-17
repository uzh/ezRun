# App wrapper for Seurat Xenium Analysis
# 
# This app processes Xenium data using Seurat, performs QC, Clustering, and RCTD Annotation.
# It follows FGCZ best practices and recommendations.

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
  
  # Setup directory
  cwd <- getwd()
  setwdNew("SeuratXenium")
  
  # Get input paths
  dataset <- input$meta
  xeniumPaths <- input$getColumn("XeniumPath")
  names(xeniumPaths) <- rownames(dataset)
  
  # Initialize results
  seurat_objects <- list()
  
  # Process each sample
  for (sampleName in names(xeniumPaths)) {
    xeniumPath <- file.path(param$dataRoot, xeniumPaths[sampleName])
    
    ezWrite(paste("Processing sample:", sampleName, "at", xeniumPath), "log.txt", append=TRUE)
    
    if (!dir.exists(xeniumPath)) {
      ezWrite(paste("Error: Path not found:", xeniumPath), "log.txt", append=TRUE)
      next
    }
    
    # 1. Load Xenium Data
    # Seurat v5 LoadXenium
    tryCatch({
      # Assuming standard Xenium output structure
      # We usually look for 'transcripts.csv.gz' or 'transcripts.parquet' and 'cell_feature_matrix'
      # LoadXenium handles the complexity
      sdata <- LoadXenium(data.dir = xeniumPath, fov = "fov")

      # Load additional cell metadata not loaded by Seurat (cell_area, nucleus_area)
      cells_file <- file.path(xeniumPath, "cells.csv.gz")
      if (file.exists(cells_file)) {
        cells_meta <- data.table::fread(cells_file, select = c("cell_id", "cell_area", "nucleus_area"))
        # Match to Seurat cell names (before RenameCells, cell names match cell_id directly)
        meta_idx <- match(colnames(sdata), cells_meta$cell_id)

        if (sum(!is.na(meta_idx)) > 0) {
          sdata$cell_area <- cells_meta$cell_area[meta_idx]
          sdata$nucleus_area <- cells_meta$nucleus_area[meta_idx]
          ezWrite(paste("Added cell_area and nucleus_area to", sum(!is.na(meta_idx)), "cells"), "log.txt", append = TRUE)
        }
      }

      # Rename cells to avoid conflicts if merging later (though here we might analyze per sample or merge)
      sdata <- RenameCells(sdata, add.cell.id = sampleName)
      sdata$orig.ident <- sampleName
      
      # 2. QC
      # Basic QC: Filter based on counts and features
      # Using params if available, else defaults
      min_counts <- ifelse(is.null(param$minCounts), 10, param$minCounts)
      min_features <- ifelse(is.null(param$minFeatures), 5, param$minFeatures)
      
      sdata <- subset(sdata, subset = nCount_Xenium > min_counts & nFeature_Xenium > min_features)
      
      # 3. Normalization (Xenium-specific with median counts as scale factor)
      sdata <- NormalizeData(sdata, scale.factor = median(sdata$nCount_Xenium))
      
      # 4. Feature Selection & Scaling
      sdata <- FindVariableFeatures(sdata, selection.method = "vst", nfeatures = 2000)
      sdata <- ScaleData(sdata)
      
      # 5. Dimension Reduction
      sdata <- RunPCA(sdata, verbose = FALSE)
      sdata <- RunUMAP(sdata, dims = 1:30, verbose = FALSE)
      
      # 6. Clustering
      res <- ifelse(is.null(param$Cluster_resolution), 0.5, param$Cluster_resolution)
      sdata <- FindNeighbors(sdata, dims = 1:30, verbose = FALSE)
      sdata <- FindClusters(sdata, resolution = res, verbose = FALSE)
      
      # 7. RCTD Annotation
      if (!is.null(param$rctdReference) && param$rctdReference != "" && param$rctdReference != "None") {
        # Extract folder name (strip species annotation)
        ref_folder <- sub(" \\([^)]+\\)$", "", param$rctdReference)
        ref_dir <- file.path("/srv/GT/databases/RCTD_References", ref_folder)

        # Find reference file
        if (dir.exists(ref_dir)) {
          available_files <- list.files(ref_dir, pattern = "\\.rds$", full.names = FALSE)
          ezWrite(paste("Available RCTD files in", ref_folder, ":",
                        paste(available_files, collapse = ", ")), "log.txt", append = TRUE)

          # Use specified file or first available
          if (!is.null(param$rctdFile) && param$rctdFile != "" &&
              param$rctdFile %in% available_files) {
            ref_path <- file.path(ref_dir, param$rctdFile)
          } else if (length(available_files) > 0) {
            ref_path <- file.path(ref_dir, available_files[1])
            ezWrite(paste("Using default file:", available_files[1]), "log.txt", append = TRUE)
          } else {
            ref_path <- NULL
            ezWrite("No .rds files found in reference folder", "log.txt", append = TRUE)
          }
        } else {
          ref_path <- NULL
          ezWrite(paste("Reference folder not found:", ref_dir), "log.txt", append = TRUE)
        }

        if (!is.null(ref_path) && file.exists(ref_path)) {
          ezWrite(paste("Running RCTD with reference:", ref_path), "log.txt", append=TRUE)
          
          # Load Reference
          ref_obj <- readRDS(ref_path)
          
          # Check if it is a spacexr Reference object
          if (!inherits(ref_obj, "Reference")) {
             # Try to convert if it's a list or something else (based on previous exploration)
             if (is.list(ref_obj) && "reference" %in% names(ref_obj)) {
               ref_obj <- ref_obj$reference
             } else {
               warning("Loaded object is not a valid RCTD Reference. Skipping RCTD.")
               ref_obj <- NULL
             }
          }
          
          if (!is.null(ref_obj)) {
            # Prepare Query (SpatialRNA object)
            # spacexr requires specific format.
            # We extract counts from Seurat.
            # Use layer instead of slot for Seurat v5 compatibility
            counts <- GetAssayData(sdata, assay = "Xenium", layer = "counts")
            coords <- GetTissueCoordinates(sdata)
            
            # Ensure coords match counts columns
            # Seurat v5 GetTissueCoordinates returns a dataframe with x, y (and maybe cell names as rownames)
            # spacexr expects a dataframe with x, y and rownames as cell names
            
            # Check for coord column names (x, y)
            if ("x" %in% colnames(coords) && "y" %in% colnames(coords)) {
               # Check if cell names are in a column "cell"
               if ("cell" %in% colnames(coords)) {
                   rownames(coords) <- coords$cell
               }
               coords <- coords[, c("x", "y")]
            } else {
               # Fallback for some seurat versions
               colnames(coords)[1:2] <- c("x", "y")
               if ("cell" %in% colnames(coords)) {
                   rownames(coords) <- coords$cell
               }
               coords <- coords[, c("x", "y")]
            }
            
            # Match cells
            common_cells <- intersect(colnames(counts), rownames(coords))
            counts <- counts[, common_cells]
            coords <- coords[common_cells, ]
            
            # Create SpatialRNA object
            query.puck <- SpatialRNA(coords, counts, colSums(counts))
            
            # Run RCTD
            # doublet_mode = 'doublet' is standard for RCTD, but 'full' is also option.
            # Best_practices.md suggests "RCTD specifically for doublet detection".
            # Get UMI_min parameter (default 100)
            umi_min <- ifelse(is.null(param$rctdUMImin), 100, as.numeric(param$rctdUMImin))
            myRCTD <- create.RCTD(query.puck, ref_obj, max_cores = param$cores, UMI_min = umi_min)
            myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
            
            # Extract results
            results <- myRCTD@results
            # Primary assignment
            norm_weights <- normalize_weights(results$weights) 
            cell_type_names <- myRCTD@cell_type_info$info[[2]]
            spatialRNA <- myRCTD@spatialRNA
            
            # Add to Seurat metadata
            # Primary cell type assignment (max weight type)
            max_type <- colnames(norm_weights)[max.col(norm_weights, ties.method = "first")]
            names(max_type) <- rownames(norm_weights)
            sdata <- AddMetaData(sdata, metadata = max_type, col.name = "RCTD_Main")

            # Add normalized weights as metadata columns (for visualization)
            weight_df <- as.data.frame(norm_weights)
            colnames(weight_df) <- paste0("rctd.weight.", colnames(weight_df))
            # Match cells
            common_cells <- intersect(colnames(sdata), rownames(weight_df))
            if (length(common_cells) > 0) {
              for (wt_col in colnames(weight_df)) {
                wt_vals <- weight_df[common_cells, wt_col]
                names(wt_vals) <- common_cells
                sdata <- AddMetaData(sdata, metadata = wt_vals, col.name = wt_col)
              }
            }

            # Add doublet/singlet classification (spot_class)
            if (!is.null(results$results_df)) {
               sdata <- AddMetaData(sdata, metadata = results$results_df)
            }
          }
        }
      }
      
      seurat_objects[[sampleName]] <- sdata
      
    }, error = function(e) {
      ezWrite(paste("Error processing sample", sampleName, ":", e$message), "log.txt", append=TRUE)
    })
  }
  
  if (length(seurat_objects) == 0) {
    stop("No samples processed successfully.")
  }
  
  # Merge if multiple samples?
  # The requirement says "similar to XeniumQC", which aggregates metrics but keeps reports per sample or combined.
  # Usually standard Seurat Apps process one dataset (which might be merged) or return a list.
  # For this app, let's assume we output one merged object if feasible, or just save the list.
  # FGCZ apps typically produce one report.
  
  # Let's merge if > 1 sample
  if (length(seurat_objects) > 1) {
    scData <- merge(seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = names(seurat_objects))
    # Re-process merged object
    scData <- JoinLayers(scData)
    scData <- NormalizeData(scData, scale.factor = median(scData$nCount_Xenium))
    scData <- FindVariableFeatures(scData)
    scData <- ScaleData(scData)
    scData <- RunPCA(scData)
    scData <- RunUMAP(scData, dims = 1:30)
    scData <- FindNeighbors(scData, dims = 1:30)
    scData <- FindClusters(scData, resolution = res)
  } else {
    scData <- seurat_objects[[1]]
  }
  
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
    lambda <- ifelse(is.null(param$lambda), 0.8, param$lambda)
    niche_res <- ifelse(is.null(param$Niche_resolution), 0.5, param$Niche_resolution)

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
  }, error = function(e) {
    ezWrite(paste("BANKSY analysis failed:", e$message), "log.txt", append = TRUE)
    # Create empty markers file
    writexl::write_xlsx(data.frame(), "posMarkersBanksy.xlsx")
  })

  # Export Xenium Explorer compatible CSV files
  # Format: cell_id,group (numeric cell ID extracted from cell names)
  ezWrite("Exporting Xenium Explorer CSV files...", "log.txt", append = TRUE)

  # Extract numeric cell ID from cell names
  # Cell names like "aaaaaklk-1" → extract numeric part after last dash
  # Or names like "123" → keep as is
  extract_cell_id <- function(cell_names) {
    sapply(cell_names, function(x) {
      # Try to extract number after last dash
      parts <- strsplit(as.character(x), "-")[[1]]
      last <- parts[length(parts)]
      if (grepl("^[0-9]+$", last)) return(last)
      # Otherwise try to extract any number
      nums <- regmatches(x, regexpr("[0-9]+", x))
      if (length(nums) > 0) return(nums)
      # Fallback to original
      return(x)
    }, USE.NAMES = FALSE)
  }

  cell_ids <- extract_cell_id(colnames(scData))

  # Clusters for Xenium Explorer
  write.csv(data.frame(
    cell_id = cell_ids,
    group = as.character(scData$seurat_clusters)
  ), "clusters_for_explorer.csv", row.names = FALSE)

  # Niches for Xenium Explorer (if BANKSY was run)
  if ("banksy_cluster" %in% colnames(scData@meta.data)) {
    write.csv(data.frame(
      cell_id = cell_ids,
      group = as.character(scData$banksy_cluster)
    ), "niches_for_explorer.csv", row.names = FALSE)
  }

  # Cell types for Xenium Explorer (if RCTD was run)
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

  # Generate Report
  makeRmdReport(rmdFile = "SeuratXenium.Rmd", reportTitle = paste0(param$name, " Report"))

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
        # We assume the path exists on the server where this runs.
        # If we are local and it doesn't exist, we provide a placeholder.
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
        computeSC = ezFrame(Type = "logical", DefaultValue = TRUE,
                            Description = "Compute Seurat Analysis")
      )
    }
  )
)
