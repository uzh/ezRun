#!/usr/bin/env Rscript
# ==============================================================================
# Generate Vitessce-optimized Zarr from Seurat Object
# ==============================================================================
#
# This script converts a Seurat object to a Vitessce-optimized Zarr store.
# It can be called from:
#   1. ezRun XeniumSeuratApp (post-processing step)
#   2. Command line: Rscript vitessce-utils.R <seurat_path> <output_dir>
#   3. Sourced as a library function
#
# Requirements:
#   R: anndataR, Seurat, qs2, Matrix
#   Python: anndata, zarr, scipy (via miniforge or pixi)
#
# ==============================================================================

# Helper for string repetition (avoid overwriting base `*`)
`%x%` <- function(a, b) {
  if (is.character(a) && is.numeric(b)) {
    return(paste(rep(a, b), collapse = ""))
  }
  stop("Invalid operands for %x%")
}

# Find Python executable - check multiple possible locations
find_python <- function() {
  candidates <- c(
    # Docker (Shiny app)
    "/opt/pixi-env/bin/python",
    # FGCZ compute cluster - dedicated vitessce env (if created)
    "/usr/local/ngseq/miniforge3/envs/vitessce/bin/python",
    # FGCZ compute cluster - environments with anndata
    "/usr/local/ngseq/miniforge3/envs/gi_cellrank2/bin/python",
    "/usr/local/ngseq/miniforge3/envs/gi_bin2cell/bin/python",
    "/usr/local/ngseq/miniforge3/envs/gi_cellbender_0.3.2/bin/python",
    # FGCZ compute cluster - general purpose env
    "/usr/local/ngseq/miniforge3/envs/ps_pgueguen/bin/python",
    # System Python
    "/usr/bin/python3",
    Sys.which("python3"),
    Sys.which("python")
  )

  for (path in candidates) {
    if (file.exists(path)) {
      # Verify it has anndata, zarr, and scipy (all required)
      check <- suppressWarnings(system(
        paste(shQuote(path), "-c 'import anndata; import zarr; import scipy'"),
        intern = TRUE,
        ignore.stderr = TRUE
      ))
      exit_code <- attr(check, "status")
      if (is.null(exit_code) || exit_code == 0) {
        return(path)
      }
    }
  }

  stop(
    "Could not find Python with anndata installed. Tried: ",
    paste(candidates, collapse = ", ")
  )
}

# Get Python script path (for finding seurat_to_vitessce_zarr.py)
get_vitessce_python_script <- function() {
  # First check: ezRun package installation
  ezrun_path <- system.file(
    "python",
    "seurat_to_vitessce_zarr.py",
    package = "ezRun"
  )
  if (nzchar(ezrun_path) && file.exists(ezrun_path)) {
    return(ezrun_path)
  }

  # Second check: Development location (paul-scripts)
  dev_path <- "/home/pgueguen/git/paul-scripts/Internal_Dev/SeuratXeniumApp/scripts/seurat_to_vitessce_zarr.py"
  if (file.exists(dev_path)) {
    return(dev_path)
  }

  # Third check: Same directory as this script (when sourced)
  if (sys.nframe() > 0) {
    frame_files <- sapply(seq_len(sys.nframe()), function(i) sys.frame(i)$ofile)
    frame_files <- frame_files[!sapply(frame_files, is.null)]
    if (length(frame_files) > 0) {
      script_dir <- dirname(frame_files[[length(frame_files)]])
      candidate <- file.path(script_dir, "seurat_to_vitessce_zarr.py")
      if (file.exists(candidate)) {
        return(candidate)
      }
    }
  }

  # Fourth check: Command line script location
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("--file=", "", file_arg)))
    candidate <- file.path(script_dir, "seurat_to_vitessce_zarr.py")
    if (file.exists(candidate)) {
      return(candidate)
    }
  }

  stop(
    "Could not find seurat_to_vitessce_zarr.py. ",
    "Ensure ezRun is installed or the script is in the expected location."
  )
}

#' Convert Seurat object to h5ad format
#'
#' @param seurat_obj Seurat object
#' @param h5ad_path Output path for h5ad file
#' @param nthreads Number of threads for parallel operations
#' @return Path to created h5ad file
seurat_to_h5ad <- function(seurat_obj, h5ad_path, nthreads = 8) {
  library(anndataR)
  library(Matrix)
  library(Seurat)

  # Get counts matrix - handle Seurat v5
  assay_name <- DefaultAssay(seurat_obj)
  assay_obj <- seurat_obj@assays[[assay_name]]

  ezLog("  Extracting expression matrix from assay: ", assay_name)

  if (inherits(assay_obj, "Assay5")) {
    layer_names <- names(assay_obj@layers)
    count_layers <- grep("^counts", layer_names, value = TRUE)

    if (length(count_layers) > 1) {
      ezLog("  Joining ", length(count_layers), " counts layers...")
      seurat_obj <- JoinLayers(seurat_obj, assay = assay_name)
    }
  }
  # Use layer parameter (slot is defunct in SeuratObject 5.0.0+)
  expr_matrix <- GetAssayData(seurat_obj, assay = assay_name, layer = "counts")

  ezLog(
    "  Matrix dimensions: ",
    nrow(expr_matrix),
    " genes x ",
    ncol(expr_matrix),
    " cells"
  )

  # Transpose for AnnData (cells x genes)
  ezLog("  Transposing to cells x genes format...")
  expr_matrix_t <- Matrix::t(expr_matrix)

  # Get metadata
  obs_df <- seurat_obj@meta.data
  ezLog("  Metadata columns: ", ncol(obs_df))

  # Get embeddings
  obsm_list <- list()

  # UMAP
  if ("umap" %in% names(seurat_obj@reductions)) {
    obsm_list[["X_umap"]] <- Embeddings(seurat_obj, "umap")
    ezLog("  Added UMAP embeddings")
  }

  # PCA
  if ("pca" %in% names(seurat_obj@reductions)) {
    obsm_list[["X_pca"]] <- Embeddings(seurat_obj, "pca")[,
      1:min(50, ncol(Embeddings(seurat_obj, "pca")))
    ]
    ezLog("  Added PCA embeddings (first 50 PCs)")
  }

  # Spatial coordinates from FOV/Xenium images
  if (length(seurat_obj@images) > 0) {
    tryCatch(
      {
        coords <- GetTissueCoordinates(seurat_obj)
        obsm_list[["spatial"]] <- as.matrix(coords[, c("x", "y")])
        ezLog("  Added spatial coordinates")
      },
      error = function(e) {
        ezLog("  Warning: Could not extract spatial coordinates: ", e$message)
      }
    )
  }

  # Create AnnData
  ezLog("  Creating AnnData object...")
  adata <- AnnData(
    X = expr_matrix_t,
    obs = obs_df,
    obsm = if (length(obsm_list) > 0) obsm_list else NULL
  )

  # Write to h5ad
  ezLog("  Writing h5ad: ", h5ad_path)
  # Use mode = "w" to allow overwriting existing files
  write_h5ad(adata, h5ad_path, mode = "w")

  return(h5ad_path)
}

#' Load Seurat object from various file formats
#'
#' @param seurat_path Path to Seurat file (.qs2, .qs, .rds)
#' @param nthreads Number of threads for qs/qs2 reading
#' @return Seurat object
load_seurat_object <- function(seurat_path, nthreads = 8) {
  if (!file.exists(seurat_path)) {
    stop("Seurat file not found: ", seurat_path)
  }

  ezLog("Loading Seurat object: ", seurat_path)

  fileExtension <- tolower(tools::file_ext(seurat_path))

  if (fileExtension == "qs2") {
    library(qs2)
    obj <- qs_read(seurat_path, nthreads = nthreads)
  } else if (fileExtension == "qs") {
    library(qs)
    obj <- qread(seurat_path, nthreads = nthreads)
  } else if (fileExtension == "rds") {
    obj <- readRDS(seurat_path)
  } else {
    stop("Unsupported file format: ", fileExtension)
  }

  ezLog("  Loaded: ", ncol(obj), " cells, ", nrow(obj), " features")
  return(obj)
}

#' Generate Vitessce-optimized Zarr from Seurat object
#'
#' This is the main function to call. It:
#' 1. Loads the Seurat object
#' 2. Converts to h5ad (intermediate format)
#' 3. Calls Python script to create optimized Zarr
#'
#' @param seurat_path Path to Seurat file (.qs2, .qs, .rds)
#' @param output_dir Output directory (zarr will be at output_dir/vitessce/data.zarr)
#' @param seurat_obj Optional: pre-loaded Seurat object (skips loading)
#' @param nthreads Number of threads for parallel operations
#' @param cleanup_h5ad Whether to delete intermediate h5ad file
#' @return List with zarr_path and metadata_path
#' @export
generate_vitessce_zarr <- function(
  seurat_path,
  output_dir,
  seurat_obj = NULL,
  nthreads = 8,
  cleanup_h5ad = TRUE
) {
  ezLog("=" %x% 60)
  ezLog("Generating Vitessce-optimized Zarr")
  ezLog("=" %x% 60)
  ezLog("Input: ", seurat_path)
  ezLog("Output: ", output_dir)

  # Find Python with anndata
  python_cmd <- find_python()
  ezLog("Using Python: ", python_cmd)

  # Find Python optimization script
  python_script <- get_vitessce_python_script()
  ezLog("Using Python script: ", python_script)

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Step 1: Load Seurat object if not provided
  if (is.null(seurat_obj)) {
    ezLog("\nStep 1/3: Loading Seurat object...")
    seurat_obj <- load_seurat_object(seurat_path, nthreads = nthreads)
  } else {
    ezLog("\nStep 1/3: Using provided Seurat object")
    ezLog("  ", ncol(seurat_obj), " cells, ", nrow(seurat_obj), " features")
  }

  # Step 2: Convert to h5ad
  ezLog("\nStep 2/3: Converting to h5ad...")
  h5ad_path <- file.path(output_dir, "temp_conversion.h5ad")
  seurat_to_h5ad(seurat_obj, h5ad_path, nthreads = nthreads)

  # Free memory before Python conversion
  rm(seurat_obj)
  gc(verbose = FALSE)

  # Step 3: Call Python for optimized Zarr
  ezLog("\nStep 3/3: Creating optimized Zarr...")
  cmd <- sprintf(
    "%s %s %s %s",
    shQuote(python_cmd),
    shQuote(python_script),
    shQuote(h5ad_path),
    shQuote(output_dir)
  )

  result <- system(cmd, intern = TRUE)
  cat(paste(result, collapse = "\n"), "\n")

  # Check for Python script success
  zarr_path <- file.path(output_dir, "vitessce", "data.zarr")
  metadata_path <- file.path(output_dir, "vitessce", "metadata.json")

  if (!dir.exists(zarr_path)) {
    stop("Zarr creation failed - directory not found: ", zarr_path)
  }

  if (!file.exists(metadata_path)) {
    stop("Metadata file not created: ", metadata_path)
  }

  # Cleanup intermediate h5ad
  if (cleanup_h5ad && file.exists(h5ad_path)) {
    unlink(h5ad_path)
    ezLog("Cleaned up intermediate h5ad file")
  }

  ezLog("\n", "=" %x% 60)
  ezLog("Vitessce Zarr generation complete!")
  ezLog("  Zarr: ", zarr_path)
  ezLog("  Metadata: ", metadata_path)
  ezLog("=" %x% 60)

  return(list(
    zarr_path = zarr_path,
    metadata_path = metadata_path
  ))
}

# ==============================================================================
# COMMAND LINE INTERFACE
# ==============================================================================

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 2) {
  args <- commandArgs(trailingOnly = TRUE)
  seurat_path <- args[1]
  output_dir <- args[2]

  ezLog("Running from command line")
  ezLog("  Seurat: ", seurat_path)
  ezLog("  Output: ", output_dir)

  result <- generate_vitessce_zarr(seurat_path, output_dir)

  ezLog("\nSuccess! Zarr at: ", result$zarr_path)
}
