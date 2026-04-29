###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Locate the CellRanger filtered feature H5 from a CountMatrix path.
##' @description The dataset.tsv `CountMatrix` column may point either at the matrix
##'   mtx directory (`*/sample_filtered_feature_bc_matrix/`, Read10X layout) or at
##'   a directory that already contains the H5 file. This helper returns the first
##'   matching `*filtered_feature_bc_matrix.h5` found in either the path itself or
##'   its parent directory.
##' @param countMatrixPath Path to the CountMatrix value.
##' @return Character path to the H5 file, or character(0) if none found.
##' @keywords internal
findFilteredH5 <- function(countMatrixPath) {
  candidates <- c(countMatrixPath, dirname(countMatrixPath))
  for (d in candidates) {
    if (!dir.exists(d)) next
    h5 <- list.files(d, pattern = "filtered_feature_bc_matrix\\.h5$",
                     full.names = TRUE, recursive = FALSE)
    if (length(h5) > 0) return(h5[1])
  }
  character(0)
}

##' @title Detect modalities present in a single-cell input directory.
##' @description Inspects a CellRanger / CellRanger Multi / CellRanger ARC output
##'   directory and reports which modalities are available for downstream
##'   processing. Robust to the CountMatrix path being either the matrix mtx
##'   directory (Read10X layout) or a count parent directory.
##' @param countMatrixPath Path to a directory containing or adjacent to
##'   CellRanger count outputs (the value of the dataset.tsv `CountMatrix` column).
##' @return Named list with logical flags: hasRNA, hasADT, hasVDJ_T, hasVDJ_B, hasATAC.
##' @export
detectModalities <- function(countMatrixPath) {
  stopifnot(dir.exists(countMatrixPath))

  h5 <- findFilteredH5(countMatrixPath)
  hasRNA <- length(h5) > 0

  hasADT <- hasRNA && h5HasAntibodyCapture(h5)

  # VDJ siblings can live one or two levels up from the CountMatrix.
  # Plain CellRanger:        <sample>/<count_or_outs>/vdj_t/...
  # CellRanger Multi (FGCZ): per_sample_outs/<sample>-cellRanger/{count,vdj_t}/...
  vdj_candidates <- unique(c(dirname(countMatrixPath),
                             dirname(dirname(countMatrixPath))))
  hasVDJ_T <- any(vapply(vdj_candidates, function(p) {
    file.exists(file.path(p, "vdj_t", "filtered_contig_annotations.csv"))
  }, logical(1)))
  hasVDJ_B <- any(vapply(vdj_candidates, function(p) {
    file.exists(file.path(p, "vdj_b", "filtered_contig_annotations.csv"))
  }, logical(1)))

  # ATAC siblings live alongside or just above the count matrix.
  atac_candidates <- unique(c(countMatrixPath, dirname(countMatrixPath)))
  hasATAC <- any(vapply(atac_candidates, function(p) {
    file.exists(file.path(p, "atac_fragments.tsv.gz")) &&
      file.exists(file.path(p, "atac_peaks.bed"))
  }, logical(1)))

  list(hasRNA = hasRNA, hasADT = hasADT,
       hasVDJ_T = hasVDJ_T, hasVDJ_B = hasVDJ_B,
       hasATAC = hasATAC)
}

##' @title Extract Antibody Capture (ADT) counts from a CellRanger H5.
##' @description Reads a multi-modal CellRanger filtered_feature_bc_matrix.h5 and
##'   returns the Antibody Capture submatrix with `<sample>_<barcode>` colnames so
##'   it aligns with a Seurat object built upstream (which prefixes barcodes by
##'   sample name). When the H5 has only a single feature type, that matrix is
##'   returned unchanged when it is Antibody Capture, otherwise an empty matrix.
##' @param h5path Path to the CellRanger HDF5.
##' @param sampleName Sample identifier used to prefix barcodes (matches
##'   `obj$Sample` / `obj$cellBarcode`).
##' @return Sparse matrix with ADT features as rows and prefixed cell barcodes
##'   as columns, or NULL if no ADT features are present.
##' @export
readADTCounts <- function(h5path, sampleName) {
  stopifnot(file.exists(h5path))
  cts <- Seurat::Read10X_h5(h5path, use.names = TRUE, unique.features = TRUE)
  if (is.list(cts)) {
    if (!"Antibody Capture" %in% names(cts)) return(NULL)
    adt <- cts[["Antibody Capture"]]
  } else {
    # Single feature type - inspect the H5 directly to confirm it's ADT.
    if (!h5HasAntibodyCapture(h5path)) return(NULL)
    adt <- cts
  }
  if (!is.null(adt) && ncol(adt) > 0) {
    colnames(adt) <- paste0(sampleName, "_", colnames(adt))
  }
  adt
}

##' @title Find the directory holding a CellRanger Multi sample's vdj_* sibling.
##' @description Returns the directory that contains `vdj_t/` and/or `vdj_b/`
##'   relative to the CountMatrix path. Mirrors the search logic in
##'   `detectModalities()`.
##' @param countMatrixPath Path to the CountMatrix value.
##' @return Character path of the parent directory holding vdj_*/, or character(0).
##' @keywords internal
findVDJSiblingDir <- function(countMatrixPath) {
  candidates <- unique(c(dirname(countMatrixPath), dirname(dirname(countMatrixPath))))
  for (p in candidates) {
    if (dir.exists(file.path(p, "vdj_t")) || dir.exists(file.path(p, "vdj_b"))) {
      return(p)
    }
  }
  character(0)
}

##' @title Locate a CellRanger VDJ filtered_contig_annotations.csv.
##' @param countMatrixPath Path to the CountMatrix value.
##' @param chain "T" for vdj_t, "B" for vdj_b.
##' @return Character path to the CSV, or character(0) if absent.
##' @keywords internal
findVDJContigCsv <- function(countMatrixPath, chain = c("T", "B")) {
  chain <- match.arg(chain)
  parent <- findVDJSiblingDir(countMatrixPath)
  if (length(parent) == 0) return(character(0))
  sub <- if (chain == "T") "vdj_t" else "vdj_b"
  csv <- file.path(parent, sub, "filtered_contig_annotations.csv")
  if (file.exists(csv)) csv else character(0)
}

##' @title Prefix VDJ contig barcodes with a sample name to match Seurat colnames.
##' @description CellRanger VDJ contig CSVs report bare cell barcodes
##'   (`AAACCTGCATTCTCAT-1`); ScSeurat prefixes barcodes with the sample
##'   (`<Sample>_<barcode>`). This helper rewrites the `barcode` column so
##'   downstream `combineTCR/BCR` outputs join cleanly with the Seurat object.
##' @param df data.frame loaded from filtered_contig_annotations.csv.
##' @param sample Character sample name to prepend.
##' @return data.frame with `barcode` and (if present) `contig_id` rewritten.
##' @keywords internal
prefixVDJBarcodes <- function(df, sample) {
  stopifnot("barcode" %in% colnames(df))
  if (any(grepl(paste0("^", sample, "_"), df$barcode))) {
    return(df)  # already prefixed
  }
  df$barcode <- paste0(sample, "_", df$barcode)
  if ("contig_id" %in% colnames(df)) {
    df$contig_id <- paste0(sample, "_", df$contig_id)
  }
  df
}

##' @title Attach VDJ-T / VDJ-B clonal annotations to a Seurat object.
##' @description Wraps `scRepertoire::combineTCR` / `combineBCR` and
##'   `combineExpression` so `obj@meta.data` gains CTaa / CTgene / cloneSize /
##'   Frequency columns. Returns the (possibly subset) Seurat object plus the
##'   raw combined clones list as an attribute.
##' @param obj Seurat object (already annotated).
##' @param vdjTPath Optional path to vdj_t/filtered_contig_annotations.csv.
##' @param vdjBPath Optional path to vdj_b/filtered_contig_annotations.csv.
##' @param sampleName Sample identifier matching obj$Sample / column prefix.
##' @return Seurat object with VDJ metadata. The combined clones list is
##'   attached as `attr(obj, "vdjCombined")`; saving as `qs2` strips
##'   attributes, so callers may also export the raw list to TSV.
##' @export
processVDJ <- function(obj, vdjTPath = NULL, vdjBPath = NULL, sampleName) {
  if (!requireNamespace("scRepertoire", quietly = TRUE)) {
    stop("scRepertoire is required for processVDJ()")
  }
  contigs <- list()
  chains <- character()
  if (!is.null(vdjTPath) && nzchar(vdjTPath) && file.exists(vdjTPath)) {
    contigs[[sampleName]] <- prefixVDJBarcodes(
      utils::read.csv(vdjTPath, stringsAsFactors = FALSE), sampleName
    )
    chains <- "T"
  } else if (!is.null(vdjBPath) && nzchar(vdjBPath) && file.exists(vdjBPath)) {
    contigs[[sampleName]] <- prefixVDJBarcodes(
      utils::read.csv(vdjBPath, stringsAsFactors = FALSE), sampleName
    )
    chains <- "B"
  } else {
    warning("processVDJ called with no usable contig CSV; returning obj unchanged.")
    return(obj)
  }

  combined <- if (chains == "T") {
    scRepertoire::combineTCR(contigs, samples = sampleName,
                             removeNA = FALSE, removeMulti = FALSE,
                             filterMulti = FALSE)
  } else {
    scRepertoire::combineBCR(contigs, samples = sampleName,
                             removeNA = FALSE, removeMulti = FALSE,
                             filterMulti = FALSE)
  }

  # combineTCR/BCR prefixes barcodes again as <sample>_<barcode>; collapse the
  # double prefix back to the Seurat-style single prefix so combineExpression
  # finds the cells.
  for (i in seq_along(combined)) {
    bc <- combined[[i]]$barcode
    dup <- paste0(sampleName, "_", sampleName, "_")
    bc <- sub(dup, paste0(sampleName, "_"), bc, fixed = TRUE)
    combined[[i]]$barcode <- bc
  }

  obj <- scRepertoire::combineExpression(
    combined, obj,
    cloneCall = "strict",
    group.by = "sample",
    proportion = FALSE,
    cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500)
  )
  # exploreSC crashes on NA factor values when group.by/color uses cloneSize
  # (server-import.R hits `if (grepl(NA, ...))`). Convert NAs to an explicit
  # "No clonotype" level so the saved scMultiData.qs2 stays Shiny-safe.
  if ("cloneSize" %in% colnames(obj@meta.data)) {
    pal_levels <- names(cloneSizePalette())
    cs <- as.character(obj$cloneSize)
    cs[is.na(cs)] <- "No clonotype"
    obj$cloneSize <- factor(cs, levels = pal_levels)
  }
  attr(obj, "vdjCombined") <- combined
  attr(obj, "vdjChain")    <- chains
  obj
}

##' @title Pick the best cell-type / annotation metadata column on a Seurat object.
##' @description Returns the first column matching a known annotation pattern,
##'   or NULL if none found. Used by the report to pick a UMAP `group.by` column
##'   for "by cell type" panels.
##' @param obj Seurat object.
##' @return Character column name or NULL.
##' @export
pickCellTypeColumn <- function(obj) {
  # Priority requested by FGCZ:
  # CyteType > cellxgene > AzimuthPanHuman > Azimuth > scType.
  # Other manual / experimental annotations come last so an upstream
  # CyteType column always wins when present.
  candidates <- c(
    # CyteType (Nygen Analytics CyteTypeR)
    "CyteTypeR.annotation", "CyteType_annotation", "cytetype", "CyteType",
    # cellxgene mapping (label-transfer from a cellxgene reference dataset)
    "cellxgene.labels", "cellxgene_labels", "cellxgene", "cellxgene_celltype",
    # Azimuth Pan-Human
    "azimuth_pan_human", "AzimuthPanHuman.labels", "panhuman.predicted.celltype",
    "PanHuman.celltype",
    # Azimuth tissue references
    "predicted.celltype.l2", "predicted.celltype.l1",
    "Azimuth.labels", "azimuth.labels",
    # scType
    "scType.labels", "scType_label", "sctype", "scType",
    # Manual / experimental fallbacks
    "manual_celltype", "celltype", "CellType", "cell_type",
    "Cell_Type_Experimental",
    "SingleR.labels"
  )
  hit <- intersect(candidates, colnames(obj@meta.data))
  if (length(hit) > 0) hit[1] else NULL
}

##' @title scRepertoire-style cloneSize palette for Seurat DimPlots.
##' @description Returns a 6-colour vector matching scRepertoire's default
##'   yellow-to-black ramp (Hyperexpanded -> None) so DimPlot(group.by =
##'   "cloneSize") and scRepertoire bar plots use consistent colours.
##' @keywords internal
cloneSizePalette <- function() {
  c("Hyperexpanded (100 < X <= 500)" = "#FFFE9E",
    "Large (20 < X <= 100)"          = "#F7AA2D",
    "Medium (5 < X <= 20)"           = "#E45358",
    "Small (1 < X <= 5)"             = "#A01975",
    "Single (0 < X <= 1)"            = "#4C1258",
    "None ( < X <= 0)"               = "#040404",
    "No clonotype"                   = "#BDBDBD")
}

##' @title Run WNN integration when 2+ dimensional modalities are present.
##' @description Wraps `Seurat::FindMultiModalNeighbors` selecting from
##'   {`pca`, `adt.pca`, `lsi`} based on which reductions exist on the object.
##'   Computes a `wnn.umap` and adds per-modality weights to metadata.
##' @param obj Seurat object.
##' @param dims_rna,dims_adt,dims_atac Dimension ranges per modality.
##' @return Seurat object with `reductions$wnn.umap`, `wsnn` graph, and
##'   per-modality `*.weight` columns. Returns obj unchanged when fewer than
##'   2 dimensional modalities exist.
##' @export
runWNN <- function(obj, dims_rna = 1:20, dims_adt = NULL, dims_atac = NULL,
                   resolution = 0.5) {
  reds <- names(obj@reductions)
  components <- list()
  if ("pca" %in% reds)     components$rna  <- list(red = "pca",     dims = dims_rna)
  if ("adt.pca" %in% reds) components$adt  <- list(red = "adt.pca", dims = dims_adt %||% seq_len(min(18L, ncol(obj[["adt.pca"]]))))
  if ("lsi" %in% reds)     components$atac <- list(red = "lsi",     dims = dims_atac %||% 2:min(30L, ncol(obj[["lsi"]])))
  if (length(components) < 2L) {
    message("WNN skipped: fewer than 2 dimensional modalities (", length(components), ").")
    return(obj)
  }

  reductions <- vapply(components, `[[`, "red", FUN.VALUE = character(1))
  dims_list  <- lapply(components, `[[`, "dims")
  message("Running WNN on: ", paste(reductions, collapse = " + "),
          " | resolution = ", resolution)

  obj <- Seurat::FindMultiModalNeighbors(
    obj,
    reduction.list = unname(as.list(reductions)),
    dims.list      = unname(dims_list),
    knn.graph.name = "wknn",
    snn.graph.name = "wsnn",
    weighted.nn.name = "weighted.nn",
    verbose = FALSE
  )
  n_neighbors <- min(30L, max(2L, ncol(obj) - 1L))
  obj <- Seurat::RunUMAP(obj, nn.name = "weighted.nn",
                         reduction.name = "wnn.umap",
                         reduction.key  = "wnnUMAP_",
                         n.neighbors = n_neighbors, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, graph.name = "wsnn", algorithm = 3,
                              resolution = resolution, verbose = FALSE)
  cluster_col <- paste0("wsnn_res.", resolution)
  cl <- obj[[cluster_col, drop = TRUE]]
  if (is.factor(cl) && all(!is.na(suppressWarnings(as.integer(levels(cl)))))) {
    sorted_levels <- as.character(sort(as.integer(levels(cl))))
    cl <- factor(cl, levels = sorted_levels)
    obj[[cluster_col]] <- cl
    obj$seurat_clusters <- cl
    Seurat::Idents(obj) <- cluster_col
  }
  obj@misc$wnn_cluster_col <- cluster_col
  obj@misc$wnn_resolution  <- resolution
  obj
}

##' @title clonalOccupy with all cloneSize bins in the legend.
##' @description Wraps `scRepertoire::clonalOccupy()` so empty cloneSize bins
##'   still appear (with their colour) in the legend. scRepertoire drops rows
##'   with `n == 0` before plotting, which strips empty levels from the
##'   geom data and prevents ggplot2 from emitting their fill swatches even
##'   with `drop = FALSE`. We pull the data via `exportTable = TRUE`, inject
##'   `n = 0L` phantom rows for missing levels, then rebuild the bar plot
##'   with the canonical `cloneSizePalette()`.
##' @param obj Seurat object with scRepertoire metadata.
##' @param x.axis Metadata column for the x axis.
##' @param palette Named colour vector; defaults to `cloneSizePalette()`.
##' @return ggplot object.
##' @export
cloneOccupyFull <- function(obj, x.axis, palette = cloneSizePalette(),
                            proportion = FALSE) {
  df <- scRepertoire::clonalOccupy(obj, x.axis = x.axis, exportTable = TRUE)
  if (!is.factor(df[[x.axis]])) df[[x.axis]] <- factor(df[[x.axis]])
  df$cloneSize <- factor(df$cloneSize, levels = names(palette))
  empty <- setdiff(names(palette), as.character(df$cloneSize))
  if (length(empty) > 0L) {
    pad <- data.frame(matrix(NA, nrow = length(empty), ncol = 0))
    pad[[x.axis]] <- factor(rep(levels(df[[x.axis]])[1], length(empty)),
                            levels = levels(df[[x.axis]]))
    pad$cloneSize <- factor(empty, levels = names(palette))
    pad$n <- 0L
    df <- rbind(df, pad[, colnames(df), drop = FALSE])
  }
  if (isTRUE(proportion)) {
    totals <- stats::ave(df$n, df[[x.axis]], FUN = sum)
    df$pct <- ifelse(totals > 0, 100 * df$n / totals, 0)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x.axis]], y = .data[["pct"]],
                                          fill = .data[["cloneSize"]])) +
      ggplot2::geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
      ggplot2::geom_text(data = df[df$n > 0, ],
                         ggplot2::aes(label = sprintf("%.0f%%", .data[["pct"]])),
                         position = ggplot2::position_stack(vjust = 0.5),
                         size = 3) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%"),
                                  expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::ylab("% of cells per group")
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x.axis]], y = .data[["n"]],
                                          fill = .data[["cloneSize"]])) +
      ggplot2::geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
      ggplot2::geom_text(data = df[df$n > 0, ],
                         ggplot2::aes(label = .data[["n"]]),
                         position = ggplot2::position_stack(vjust = 0.5)) +
      ggplot2::ylab("Single Cells")
  }
  p +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE,
                               limits = names(palette)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.title.x = ggplot2::element_blank())
}

##' @title DimPlot of cloneSize with all 6 swatches in the legend.
##' @description Seurat's `DimPlot(group.by = "cloneSize")` only emits legend
##'   keys for levels actually drawn, so empty bins (Hyperexpanded, Small,
##'   Rare, None) silently lose their colour swatch. We rebuild the plot
##'   from the embedding + metadata, append one phantom point per missing
##'   level so ggplot2 generates a key for it, then hide the phantoms via
##'   `override.aes(alpha = 1)`. The visible points are unchanged.
##' @param obj Seurat object with scRepertoire metadata + reduction.
##' @param reduction Reduction name (e.g. "umap", "adt.umap", "wnn.umap").
##' @param palette Named colour vector; defaults to `cloneSizePalette()`.
##' @param point_size Point size for visible cells.
##' @param title Plot title (optional).
##' @return ggplot object, or NULL when the reduction or cloneSize column
##'   is missing.
##' @export
cloneDimPlot <- function(obj, reduction, palette = cloneSizePalette(),
                         point_size = 0.5, title = NULL) {
  if (!reduction %in% names(obj@reductions)) return(NULL)
  if (!"cloneSize" %in% colnames(obj@meta.data)) return(NULL)
  emb <- as.data.frame(Seurat::Embeddings(obj, reduction))
  axis_x <- colnames(emb)[1]; axis_y <- colnames(emb)[2]
  emb$cloneSize <- factor(obj$cloneSize, levels = names(palette))
  empty <- setdiff(names(palette), as.character(unique(emb$cloneSize)))
  phantom <- if (length(empty) > 0) {
    data.frame(setNames(list(rep(mean(emb[[axis_x]], na.rm = TRUE), length(empty)),
                              rep(mean(emb[[axis_y]], na.rm = TRUE), length(empty))),
                         c(axis_x, axis_y)),
               cloneSize = factor(empty, levels = names(palette)))
  } else NULL
  emb_ord <- emb[order(!is.na(emb$cloneSize)), ]
  p <- ggplot2::ggplot(emb_ord,
                       ggplot2::aes(x = .data[[axis_x]], y = .data[[axis_y]],
                                    color = .data[["cloneSize"]])) +
    ggplot2::geom_point(size = point_size)
  if (!is.null(phantom)) {
    p <- p + ggplot2::geom_point(data = phantom, alpha = 0)
  }
  p <- p +
    ggplot2::scale_color_manual(values = palette, na.value = "grey85",
                                drop = FALSE, limits = names(palette)) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 6, alpha = 1, stroke = 0.5))) +
    ggplot2::theme_classic() +
    ggplot2::theme(aspect.ratio = 1,
                   legend.key.size = ggplot2::unit(1.2, "lines"),
                   legend.text = ggplot2::element_text(size = 10))
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}

##' @title Load a BD Rhapsody pre-built Seurat object.
##' @description BD Rhapsody pipelines produce a `<sample>_Seurat.rds` next to
##'   the per-sample MEX folder, already populated with RNA and (when AbSeq is
##'   present) ADT assays plus tSNE/UMAP. This loader reads that file, prefixes
##'   barcodes `<sample>_<barcode>`, and aligns metadata to ScMultiOmics
##'   conventions. If the RDS is absent, returns NULL so callers can fall
##'   back.
##' @param resultDir Path to the BD Rhapsody result directory (one sample).
##' @param sampleName Sample name to prefix barcodes; ScSeurat-style.
##' @return Seurat object or NULL.
##' @export
loadBDRhapsody <- function(resultDir, sampleName) {
  rds <- list.files(resultDir, pattern = "_Seurat\\.rds$", full.names = TRUE,
                    recursive = FALSE)
  if (length(rds) == 0) {
    message("No *_Seurat.rds found in ", resultDir)
    return(NULL)
  }
  obj <- readRDS(rds[1])
  if (!inherits(obj, "Seurat")) return(NULL)
  obj <- Seurat::RenameCells(obj,
                             new.names = paste0(sampleName, "_", colnames(obj)))
  obj$Sample <- sampleName
  obj$cellBarcode <- sub(paste0("^", sampleName, "_"), "", colnames(obj))
  Seurat::DefaultAssay(obj) <- "RNA"
  # BD's RDS ships UMAP/tSNE but no `pca` reduction. Compute one so downstream
  # WNN integration can see RNA + ADT as two dimensional modalities.
  if (!"pca" %in% names(obj@reductions)) {
    message("BD Seurat has no PCA - computing RNA PCA + neighbors + clusters.")
    if (length(Seurat::VariableFeatures(obj)) == 0) {
      obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
    }
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- Seurat::FindNeighbors(obj, dims = 1:30, verbose = FALSE)
    obj <- Seurat::FindClusters(obj, resolution = 0.5, verbose = FALSE)
  }
  obj
}

##' @title Locate the CellRanger ARC ATAC fragments + peaks for a CountMatrix.
##' @param countMatrixPath Path to the CountMatrix value.
##' @return Named list with `fragments` and `peaks` paths, or NULL if absent.
##' @keywords internal
findATACFiles <- function(countMatrixPath) {
  candidates <- unique(c(countMatrixPath, dirname(countMatrixPath)))
  for (p in candidates) {
    fr <- file.path(p, "atac_fragments.tsv.gz")
    pk <- file.path(p, "atac_peaks.bed")
    if (file.exists(fr) && file.exists(pk)) {
      return(list(fragments = fr, peaks = pk))
    }
  }
  NULL
}

##' @title Resolve gene annotation GRanges for a Signac chromatin assay.
##' @description Tries the EnsDb package matching the param$refBuild species; the
##'   caller may also supply a pre-loaded annotation via param$annotation.
##' @param refBuild ezRun refBuild string (or species name).
##' @return GRanges with seqlevelsStyle "UCSC", or NULL if no EnsDb available.
##' @keywords internal
getATACAnnotation <- function(refBuild = NULL) {
  if (is.null(refBuild)) refBuild <- ""
  human <- grepl("Homo_sapiens|GRCh38|hg38|human", refBuild, ignore.case = TRUE)
  mouse <- grepl("Mus_musculus|GRCm|mm10|mm39|mouse", refBuild, ignore.case = TRUE)
  ensdb_pkg <- if (human) "EnsDb.Hsapiens.v86" else if (mouse) "EnsDb.Mmusculus.v79" else NULL
  if (is.null(ensdb_pkg) || !requireNamespace(ensdb_pkg, quietly = TRUE)) {
    return(NULL)
  }
  ensdb <- getExportedValue(ensdb_pkg, ensdb_pkg)
  ann <- Signac::GetGRangesFromEnsDb(ensdb = ensdb)
  GenomeInfoDb::seqlevelsStyle(ann) <- "UCSC"
  ann
}

##' @title Attach an ATAC assay to a multiome Seurat object.
##' @description Reads the CellRanger ARC peak matrix from the H5, builds a
##'   `Signac::ChromatinAssay` (with fragments + annotation), runs TF-IDF + SVD,
##'   computes an ATAC UMAP, and adds a `GeneActivity` assay.
##' @param obj Seurat object with RNA assay (cells already prefixed `<Sample>_<barcode>`).
##' @param fragmentsPath Path to atac_fragments.tsv.gz.
##' @param peaksPath Path to atac_peaks.bed (currently informational; the H5
##'   matrix carries the called peak set).
##' @param refBuild ezRun refBuild string used to resolve EnsDb annotation.
##' @param sampleName Sample name to prefix barcodes; must match `obj$Sample`.
##' @return Seurat object with `assays$ATAC`, `reductions$lsi`, `reductions$atac.umap`,
##'   and `assays$GeneActivity`.
##' @export
processATAC <- function(obj, fragmentsPath, peaksPath, refBuild = NULL,
                        sampleName) {
  if (!requireNamespace("Signac", quietly = TRUE)) {
    stop("Signac is required for processATAC()")
  }
  stopifnot(file.exists(fragmentsPath))

  cm_dir <- dirname(fragmentsPath)
  h5 <- list.files(cm_dir, pattern = "filtered_feature_bc_matrix\\.h5$", full.names = TRUE)
  if (length(h5) == 0) stop("Could not find filtered_feature_bc_matrix.h5 next to ", fragmentsPath)
  cts <- Seurat::Read10X_h5(h5[1], use.names = TRUE, unique.features = TRUE)
  if (!is.list(cts) || !"Peaks" %in% names(cts)) {
    stop("Expected a 'Peaks' submatrix in ", h5[1], "; got: ", paste(names(cts), collapse = ", "))
  }
  peak_counts <- cts[["Peaks"]]
  colnames(peak_counts) <- paste0(sampleName, "_", colnames(peak_counts))

  shared <- intersect(colnames(peak_counts), colnames(obj))
  if (length(shared) < ncol(obj) * 0.5) {
    warning(sprintf("Only %d/%d cells overlap ATAC barcodes - check sample prefixing.",
                    length(shared), ncol(obj)))
  }
  obj <- obj[, shared]
  peak_counts <- peak_counts[, shared, drop = FALSE]

  ann <- getATACAnnotation(refBuild)

  # Build a fragments object with bare-barcode -> sample-prefixed-barcode map
  bare_bcs <- sub(paste0("^", sampleName, "_"), "", shared)
  cell_map <- setNames(bare_bcs, shared)
  frag_obj <- Signac::CreateFragmentObject(path = fragmentsPath, cells = cell_map,
                                           validate.fragments = FALSE)

  chrom <- Signac::CreateChromatinAssay(
    counts = peak_counts,
    sep = c(":", "-"),
    fragments = frag_obj,
    annotation = ann,
    min.cells = 0,
    min.features = 0
  )
  obj[["ATAC"]] <- chrom
  Seurat::DefaultAssay(obj) <- "ATAC"

  obj <- Signac::FindTopFeatures(obj, min.cutoff = "q5")
  obj <- Signac::RunTFIDF(obj)
  obj <- Signac::RunSVD(obj)

  n_lsi <- min(30L, ncol(Seurat::Embeddings(obj, "lsi")))
  n_neighbors <- min(30L, max(2L, ncol(obj) - 1L))
  obj <- Seurat::RunUMAP(obj, reduction = "lsi", dims = 2:n_lsi,
                         reduction.name = "atac.umap",
                         n.neighbors = n_neighbors, verbose = FALSE)

  if (!is.null(ann)) {
    ga <- tryCatch(Signac::GeneActivity(obj), error = function(e) {
      message("GeneActivity skipped: ", conditionMessage(e)); NULL
    })
    if (!is.null(ga)) {
      obj[["GeneActivity"]] <- Seurat::CreateAssayObject(counts = ga)
      obj <- Seurat::NormalizeData(obj, assay = "GeneActivity",
                                   normalization.method = "LogNormalize",
                                   scale.factor = median(obj$nCount_GeneActivity),
                                   verbose = FALSE)
    }
  }

  Seurat::DefaultAssay(obj) <- "RNA"
  obj
}

##' @title Add an ADT assay to a Seurat object with normalization, PCA and UMAP.
##' @description Attaches a Seurat v5 ADT assay built from raw antibody counts,
##'   normalizes (ADTnorm or CLR), runs PCA on the ADT features and a UMAP
##'   embedding stored as `adt.umap`. Cells absent from `adtCounts` are dropped.
##' @param obj Seurat object (already has RNA assay + clusters).
##' @param adtCounts Sparse matrix or matrix of ADT raw counts (features x cells).
##'   Cell barcodes must match `colnames(obj)`.
##' @param normMethod "ADTnorm" (default) or "CLR".
##' @param npcs Number of ADT PCs (clipped to nrow(adt) - 1).
##' @return Seurat object with `assays$ADT`, `reductions$adt.pca`, `reductions$adt.umap`.
##'   `DefaultAssay` is restored to "RNA" before return.
##' @export
processADT <- function(obj, adtCounts, normMethod = c("ADTnorm", "CLR"), npcs = 18) {
  normMethod <- match.arg(normMethod)
  shared <- intersect(colnames(adtCounts), colnames(obj))
  if (length(shared) < ncol(obj) * 0.5) {
    warning(sprintf("Only %d/%d cells overlap ADT barcodes - check sample prefixing.",
                    length(shared), ncol(obj)))
  }
  obj <- obj[, shared]
  adt_aligned <- adtCounts[, shared, drop = FALSE]

  obj[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = adt_aligned)
  Seurat::DefaultAssay(obj) <- "ADT"

  if (normMethod == "ADTnorm" && requireNamespace("ADTnorm", quietly = TRUE)) {
    norm_mat <- ADTnorm::ADTnorm(
      cell_x_adt = t(as.matrix(adt_aligned)),
      cell_x_feature = data.frame(sample = rep("all", ncol(adt_aligned))),
      marker_to_process = rownames(adt_aligned),
      save_outpath = tempdir(),
      study_name = "ezrun",
      save_intermediate_fig = FALSE
    )
    obj <- Seurat::SetAssayData(obj, assay = "ADT", layer = "data",
                                new.data = t(norm_mat))
  } else {
    obj <- Seurat::NormalizeData(obj, assay = "ADT",
                                 normalization.method = "CLR", margin = 2,
                                 verbose = FALSE)
  }

  npcs <- min(npcs, nrow(adt_aligned) - 1)
  obj <- Seurat::ScaleData(obj, assay = "ADT", verbose = FALSE)
  obj <- Seurat::RunPCA(obj, assay = "ADT", reduction.name = "adt.pca",
                        npcs = npcs, features = rownames(obj[["ADT"]]),
                        approx = FALSE, verbose = FALSE)
  n_neighbors <- min(30L, max(2L, ncol(obj) - 1L))
  obj <- Seurat::RunUMAP(obj, assay = "ADT", reduction = "adt.pca",
                         dims = seq_len(npcs), reduction.name = "adt.umap",
                         n.neighbors = n_neighbors, verbose = FALSE)

  Seurat::DefaultAssay(obj) <- "RNA"
  obj
}

##' @title Check whether a CellRanger HDF5 has any Antibody Capture features.
##' @param h5path Path to a CellRanger filtered_feature_bc_matrix.h5 file.
##' @return TRUE if any feature has type "Antibody Capture", FALSE otherwise
##'   (or if the file is missing / empty / unreadable).
##' @keywords internal
h5HasAntibodyCapture <- function(h5path) {
  if (!file.exists(h5path)) return(FALSE)
  if (file.size(h5path) == 0) return(FALSE)
  tryCatch({
    hf <- hdf5r::H5File$new(h5path, mode = "r")
    on.exit(hf$close_all(), add = TRUE)
    if (!hf$exists("matrix/features/feature_type")) return(FALSE)
    ft <- hf[["matrix/features/feature_type"]]$read()
    any(ft == "Antibody Capture")
  }, error = function(e) FALSE)
}
