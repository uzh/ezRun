###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## Shared BANKSY helpers (used by XeniumSeurat and SpatialSeurat).

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
    ezLog("Fetching data from layer ", layer, " from assay ", assay)
  }
  data_own <- Seurat::GetAssayData(
    object = object,
    assay = assay,
    layer = layer
  )

  if (features[1] != "all") {
    if (verbose) {
      ezLog("Subsetting by features")
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
      ezLog("Staggering locations by ", group)
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
      ezLog("Subsetting to ", length(common_cells), " common cells")
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
    ezLog("Creating Banksy matrix")
  }
  data_banksy <- c(list(data_own), har)

  if (verbose) {
    ezLog(
      "Scaling BANKSY matrix. Do not call ScaleData on assay ",
      assay_name
    )
  }
  fast_scaler <- function(mat, obj, grp, split.scl, verb) {
    # scale() divides by sd, so a zero-variance gene (a panel gene detected in
    # no cell - common, since every panel gene is a variable feature here)
    # yields a whole NaN row. Seurat's own ScaleData returns 0 in that case;
    # match it, otherwise the BANKSY assay carries NaN features into PCA.
    zeroVarSafe <- function(m) {
      s <- t(scale(t(m)))
      s[!is.finite(s)] <- 0
      s
    }
    if (!is.null(grp) && split.scl) {
      groups <- unique(unlist(obj[[grp]]))
      scaled_list <- lapply(groups, function(g) {
        idx <- which(unlist(obj[[grp]]) == g)
        idx <- idx[idx %in% seq_len(ncol(mat))]
        if (length(idx) > 0) zeroVarSafe(mat[, idx, drop = FALSE]) else NULL
      })
      scaled_list <- scaled_list[!sapply(scaled_list, is.null)]
      if (length(scaled_list) > 0) {
        do.call(cbind, scaled_list)[, colnames(mat)]
      } else {
        zeroVarSafe(mat)
      }
    } else {
      zeroVarSafe(mat)
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
    ezLog("Setting default assay to ", assay_name)
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
    ezLog("BANKSY assay created successfully")
  }
  return(object)
}
