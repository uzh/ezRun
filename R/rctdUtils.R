###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## RCTD helpers for spot-based data: spacexr (R) and rctd-py (PyTorch, same
## model). XeniumSeurat and VisiumHDSeurat keep their own inline reference
## code; these functions follow the same rules (reference coercion,
## gene-overlap check, UMI_min_sigma fallback, fixed seed).

RCTD_REFERENCE_DIR <- "/srv/GT/databases/RCTD_References"

##' @title Resolve the RCTD reference path from the app parameters
##' @description \code{rctdFile} (full path) wins over \code{rctdReference},
##'   which is relative to \code{RCTD_REFERENCE_DIR} and may carry a
##'   " (description)" suffix as used in the SUSHI dropdowns.
##' @return a path or NULL when no reference was requested.
resolve_rctd_reference_path <- function(param) {
  if (ezIsSpecified(param$rctdFile) && param$rctdFile != "") {
    return(param$rctdFile)
  }
  if (ezIsSpecified(param$rctdReference) && param$rctdReference != "None") {
    refRelative <- sub(" \\([^)]+\\)$", "", param$rctdReference)
    return(file.path(RCTD_REFERENCE_DIR, refRelative))
  }
  return(NULL)
}

##' @title Load an RCTD reference
##' @description Accepts a spacexr \code{Reference}, a list holding one in
##'   \code{$reference}, or a Seurat object with raw counts and a cell-type
##'   column.
##' @param refPath path to an .rds or .qs2 file.
##' @param cellTypeColumn cell-type column used for Seurat references.
##' @return a spacexr \code{Reference}.
load_rctd_reference <- function(refPath, cellTypeColumn = NULL) {
  if (!file.exists(refPath)) {
    stop("RCTD reference not found: ", refPath)
  }
  refObj <- ezLoadRobj(refPath)
  if (inherits(refObj, "Reference")) {
    return(refObj)
  }
  if (is.list(refObj) && "reference" %in% names(refObj)) {
    return(refObj$reference)
  }
  if (inherits(refObj, "Seurat")) {
    candidates <- c(
      cellTypeColumn,
      "author_cell_type",
      "cell_type",
      "celltype",
      "CellType",
      "class",
      "subclass"
    )
    ctCol <- intersect(candidates, colnames(refObj@meta.data))[1]
    if (is.na(ctCol)) {
      stop("No cell-type column found in the Seurat reference ", refPath)
    }
    cellTypes <- as.factor(refObj@meta.data[[ctCol]])
    names(cellTypes) <- colnames(refObj)
    counts <- Seurat::GetAssayData(refObj, layer = "counts")
    return(spacexr::Reference(counts, cellTypes, Matrix::colSums(counts)))
  }
  stop("Not a spacexr Reference, list(reference=) or Seurat object: ", refPath)
}

##' @title Name the query genes like the reference
##' @description References built from atlases often use Ensembl IDs (e.g. the
##'   Allen brain references) while the Seurat object carries (uniquified)
##'   symbols. Renames the query rows to gene IDs when that gives the larger
##'   overlap; all genes are kept, so per-spot UMI totals do not change.
##' @param counts gene x spot count matrix with symbol row names.
##' @param geneIds gene IDs named by the row names of \code{counts}.
##' @param refGenes genes of the reference.
##' @return the count matrix, with gene IDs as row names if those match better.
rename_rctd_query_genes <- function(counts, geneIds, refGenes) {
  bySymbol <- sum(rownames(counts) %in% refGenes)
  ids <- unname(geneIds[rownames(counts)])
  byId <- sum(!is.na(ids) & ids %in% refGenes)
  if (byId > bySymbol) {
    keep <- !is.na(ids) & !duplicated(ids)
    counts <- counts[keep, , drop = FALSE]
    rownames(counts) <- ids[keep]
  }
  nShared <- max(byId, bySymbol)
  if (nShared < 100) {
    stop(sprintf(
      "RCTD reference shares only %d genes with the data (wrong species or ID type?)",
      nShared
    ))
  }
  ezLog(sprintf(
    "RCTD: %d genes shared with the reference (matched by %s)",
    nShared,
    if (byId > bySymbol) "gene ID" else "symbol"
  ))
  return(counts)
}

##' @title Match query gene names to the reference
##' @description \code{rename_rctd_query_genes} followed by restriction to the
##'   shared genes.
##' @inheritParams rename_rctd_query_genes
##' @return the count matrix restricted to genes shared with the reference.
match_rctd_genes <- function(counts, geneIds, refGenes) {
  counts <- rename_rctd_query_genes(counts, geneIds, refGenes)
  shared <- intersect(rownames(counts), refGenes)
  return(counts[shared, , drop = FALSE])
}

##' @title Main cell type per bin from a doublet-mode RCTD run
##' @description The cell type with the highest weight, except for bins that
##'   RCTD classified as \code{"reject"} (it could not fit them confidently),
##'   which get NA, as in XeniumSeurat.
##' @param weights bin x cell type weight matrix.
##' @param resultsDf RCTD results_df with a \code{spot_class} column, row
##'   names = bins; NULL keeps every bin.
##' @return named character vector, one entry per row of \code{weights}.
rctd_main_type <- function(weights, resultsDf = NULL) {
  mainType <- colnames(weights)[max.col(weights, ties.method = "first")]
  names(mainType) <- rownames(weights)
  if (!is.null(resultsDf) && "spot_class" %in% colnames(resultsDf)) {
    spotClass <- as.character(resultsDf[rownames(weights), "spot_class"])
    mainType[spotClass %in% "reject"] <- NA
  }
  return(mainType)
}

##' @title UMI threshold of the spots RCTD uses to fit sigma
##' @description 300 as in spacexr, or the 75th UMI quantile when fewer than
##'   100 spots reach 300 (shallow or probe-based data).
rctd_umi_min_sigma <- function(nUmi) {
  umiMinSigma <- 300
  if (sum(nUmi > umiMinSigma) < 100) {
    umiMinSigma <- max(1, floor(stats::quantile(nUmi, 0.75, names = FALSE)))
    ezLog("RCTD: few spots above 300 UMIs; UMI_min_sigma set to ", umiMinSigma)
  }
  return(umiMinSigma)
}

##' @title Run RCTD on spots
##' @param counts gene x spot raw counts (genes already matched).
##' @param coords data.frame with columns x and y, row names = spots.
##' @param reference spacexr Reference.
##' @param cores number of cores.
##' @param umiMin minimal UMIs per spot.
##' @param doubletMode "full" for multi-cell spots (Visium SD).
##' @param seed random seed (spacexr samples spots to fit sigma).
##' @return list with \code{weights} (spot x cell type, rows sum to 1) and
##'   \code{resultsDf} (NULL in full mode).
run_rctd_spots <- function(
  counts,
  coords,
  reference,
  cores = 1,
  umiMin = 100,
  doubletMode = "full",
  seed = 42
) {
  coords <- coords[colnames(counts), c("x", "y")]
  nUmi <- Matrix::colSums(counts)
  puck <- spacexr::SpatialRNA(coords, counts, nUmi)
  umiMinSigma <- rctd_umi_min_sigma(nUmi)
  ## spacexr starts fresh worker processes that do not inherit the BLAS
  ## thread pinning of this session: without these variables every worker
  ## uses all cores of the node (observed 8 workers x 17 threads, 4x slower)
  threadVars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
  oldVars <- Sys.getenv(threadVars, unset = NA)
  do.call(Sys.setenv, as.list(stats::setNames(rep("1", length(threadVars)), threadVars)))
  on.exit(
    {
      for (v in threadVars) {
        if (is.na(oldVars[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, stats::setNames(list(oldVars[[v]]), v))
      }
    },
    add = TRUE
  )
  set.seed(seed)
  rctd <- spacexr::create.RCTD(
    puck,
    reference,
    max_cores = cores,
    UMI_min = umiMin,
    UMI_min_sigma = umiMinSigma
  )
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = doubletMode)
  weights <- as.matrix(spacexr::normalize_weights(rctd@results$weights))
  return(list(weights = weights, resultsDf = rctd@results$results_df))
}

## rctd-py through its CLI, in the gi_rctd-py env. Same model as spacexr:
## p32810-sized Visium HD runs went from 14.7 h (spacexr) to 13 min on a GPU
## and 38 min on 8 CPU cores, with 99.7% agreement on the main type.
## RCTD_PY_BIN (a directory holding `rctd`) overrides the env, for testing.
RCTD_PY_ENV <- "gi_rctd-py_0.3.8"
RCTD_PY_CONDA <- "/usr/local/ngseq/miniforge3"

## TRUE when rctd-py can run: the test override is set or the env exists.
rctdPyAvailable <- function() {
  nzchar(Sys.getenv("RCTD_PY_BIN")) ||
    file.exists(file.path(RCTD_PY_CONDA, "envs", RCTD_PY_ENV, "bin", "rctd"))
}

## TRUE when SLURM gave the job a GPU.
rctdPyHasGpu <- function() {
  gpus <- Sys.getenv("CUDA_VISIBLE_DEVICES")
  nzchar(gpus) && gpus != "NoDevFiles"
}

## torch.compile costs ~13 min up front on CPU and makes every later step
## faster. Measured on 8 CPU cores: p37810 Visium SD Rep1 (2,319 spots, full
## mode) 14.9 min compiled vs 1.3 min not; Visium HD 16 um FFPE (99,739
## bins, doublet mode) 31 min compiled vs 102 min not. Break-even is around
## 15,000 spots.
RCTD_PY_COMPILE_MIN_SPOTS <- 15000

## Compile on a GPU, and on CPU only for large inputs.
rctdPyCompile <- function(nSpots) {
  rctdPyHasGpu() || nSpots >= RCTD_PY_COMPILE_MIN_SPOTS
}

##' @title Run RCTD with rctd-py
##' @param counts gene x spot raw counts (genes already matched).
##' @param coords data.frame with columns x and y, row names = spots.
##' @param ref spacexr Reference.
##' @param umiMin minimal UMIs per spot.
##' @param mode "doublet" (Visium HD bins, cells) or "full" (multi-cell
##'   Visium SD spots).
##' @param device "auto" uses a GPU when the job has one; "cpu" never does.
##' @param umiMinSigma UMI threshold for the sigma fit; NULL keeps 300.
##' @param cores torch CPU threads; NULL leaves torch's default (all cores
##'   it can see, which oversubscribes a SLURM allocation).
##' @param compile FALSE passes --no-compile; NULL decides by rctdPyCompile
##'   (always on a GPU, on CPU only for large inputs).
##' @param workDir where the h5ad files go.
##' @return what spacexr's path yields: row-normalised \code{weights}
##'   (spots x types, only the spots RCTD fitted) and, in doublet mode, a
##'   \code{results_df} with spot_class / first_type / second_type (NULL in
##'   full mode).
runRctdPy <- function(
  counts,
  coords,
  ref,
  umiMin,
  mode = "doublet",
  device = "auto",
  umiMinSigma = NULL,
  cores = NULL,
  compile = NULL,
  env = RCTD_PY_ENV,
  workDir = file.path(getwd(), "rctd_py")
) {
  mode <- match.arg(mode, c("doublet", "full"))
  if (is.null(compile)) {
    compile <- rctdPyCompile(ncol(counts))
  }
  dir.create(workDir, showWarnings = FALSE)
  qFile <- file.path(workDir, "query.h5ad")
  rFile <- file.path(workDir, "reference.h5ad")
  outFile <- file.path(workDir, "result.h5ad")
  anndataR::write_h5ad(
    anndataR::AnnData(
      X = Matrix::t(counts),
      obs = data.frame(row.names = colnames(counts)),
      var = data.frame(row.names = rownames(counts)),
      obsm = list(spatial = as.matrix(coords[colnames(counts), c("x", "y")]))
    ),
    qFile, mode = "w"
  )
  anndataR::write_h5ad(
    anndataR::AnnData(
      X = Matrix::t(ref@counts),
      obs = data.frame(cell_type = as.character(ref@cell_types),
                       row.names = colnames(ref@counts)),
      var = data.frame(row.names = rownames(ref@counts))
    ),
    rFile, mode = "w"
  )
  binDir <- Sys.getenv("RCTD_PY_BIN")
  if (nzchar(binDir)) {
    withr::local_path(binDir, action = "prefix")
  } else {
    Herper::local_CondaEnv(env, pathToMiniConda = RCTD_PY_CONDA)
  }
  if (!is.null(cores)) {
    withr::local_envvar(c(
      OMP_NUM_THREADS = cores, MKL_NUM_THREADS = cores,
      OPENBLAS_NUM_THREADS = cores
    ))
  }
  ezSystem(paste(
    "rctd run", qFile, rFile, "--mode", mode, "--umi-min", umiMin,
    if (!is.null(umiMinSigma)) paste("--umi-min-sigma", umiMinSigma),
    "--device", device, if (!compile) "--no-compile", "-o", outFile
  ))
  res <- readRctdPyResult(outFile, mode = mode)
  ## the h5ad copies (~1.5 GB for a Visium HD section) do not belong in the
  ## result folder
  unlink(c(qFile, rFile, outFile))
  cells <- rctdResultCells(res$obsNames, colnames(counts))
  kept <- res$kept
  weights <- res$weights[kept, , drop = FALSE]
  dimnames(weights) <- list(cells[kept], res$cellTypes)
  weights <- weights / rowSums(weights)
  results_df <- NULL
  if (mode == "doublet") {
    results_df <- data.frame(
      spot_class = factor(res$spot_class[kept],
                          levels = c("reject", "singlet", "doublet_certain",
                                     "doublet_uncertain")),
      first_type = res$first_type[kept],
      second_type = res$second_type[kept],
      row.names = cells[kept]
    )
  }
  list(weights = weights, results_df = results_df)
}

## The fields of an rctd-py result h5ad, read with rhdf5. Not anndataR:
## anndata >= 0.13 writes obs/_index as a nullable-string-array, which
## anndataR 1.2 cannot decode, and it then drops the whole obs table.
## `kept` marks the spots RCTD fitted: doublet mode labels the others
## "filtered", full mode leaves their weights NaN.
readRctdPyResult <- function(file, mode = "doublet") {
  rd <- function(path) rhdf5::h5read(file, path)
  categorical <- function(col) {
    codes <- as.integer(rd(paste0("obs/", col, "/codes")))
    cats <- as.character(rd(paste0("obs/", col, "/categories")))
    cats[replace(codes + 1L, codes < 0L, NA)]
  }
  obsNames <- rd("obs/_index")
  if (is.list(obsNames)) {
    obsNames <- obsNames$values # nullable-string-array: values + mask
  }
  weights <- t(rd("obsm/rctd_weights")) # h5 row-major -> R column-major
  res <- list(
    obsNames = as.character(obsNames),
    weights = weights,
    cellTypes = as.character(rd("uns/rctd_cell_type_names"))
  )
  if (mode == "doublet") {
    res$spot_class <- categorical("rctd_spot_class")
    res$first_type <- categorical("rctd_first_type")
    res$second_type <- categorical("rctd_second_type")
    res$kept <- res$spot_class != "filtered"
  } else {
    res$kept <- rowSums(is.na(weights)) == 0
  }
  res
}

## rctd-py keeps the query's row order; refuse a result that does not.
rctdResultCells <- function(obsNames, queryCells) {
  if (length(obsNames) != length(queryCells)) {
    stop(sprintf("rctd-py returned %d rows for %d query bins",
                 length(obsNames), length(queryCells)))
  }
  if (!identical(as.character(obsNames), queryCells)) {
    stop("rctd-py result rows are not in query order")
  }
  queryCells
}
