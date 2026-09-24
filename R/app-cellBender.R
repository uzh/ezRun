###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppCellBender <-
  setRefClass(
    "EzAppCellBender",
    contains = "EzApp",
    methods = list(
      ## CellBender unconditional. DropletUtils gated on either the raw .h5 needing
      ## to be created, or on multimodal (ATAC "Peaks") input needing filtering.
      ## rhdf5 unconditional (feature-type inspection runs for every sample).
      citation = function() {
        c(
          "Fleming, S.J. et al. Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender. Nature Methods 20, 1323-1335 (2023). https://doi.org/10.1038/s41592-023-01943-7",
          "Lun, A.T.L. et al. EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome Biology 20, 63 (2019). https://doi.org/10.1186/s13059-019-1662-y",
          "Fischer, B. & Pau, G. rhdf5: R Interface to HDF5. R package version 2.56.0. https://doi.org/10.18129/B9.bioc.rhdf5"
        )
      },
      ## Defaults quoted here were read from `cellbender remove-background --help` and
      ## cellbender/remove_background/consts.py in conda env gi_cellbender_0.3.2
      ## (the package installed there is CellBender 0.3.0).
      methods_facts = function() {
        c(
          ## ezMethodCellBender input path: UnfilteredCountMatrix, else raw_feature_bc_matrix of the multi output (app-cellBender.R:105-126)
          "CellBender was run on the unfiltered (raw) droplet count matrix including empty droplets (for CellRanger Multi input without an UnfilteredCountMatrix column, the library-level raw matrix of the pool), not on the cell-filtered matrix.",
          ## dropPeaksFromH5 (app-cellBender.R:54-91), called at app-cellBender.R:159
          "For multiome inputs, ATAC peak features were removed before CellBender; all other feature types (gene expression, antibody capture, multiplexing capture) were kept and processed together.",
          ## command line built at app-cellBender.R:161-176
          "The app passed only the input, the output and either --cuda (when gpu is above 0) or --cpu-threads, plus any flags given in cmdOptions; every other setting was the CellBender default unless cmdOptions overrode it.",
          ## CellBender defaults from remove-background --help
          "Unless overridden in cmdOptions, CellBender used the full model (ambient RNA plus barcode swapping), 150 training epochs, a learning rate of 1e-4 with a one-cycle schedule, a target false positive rate (--fpr) of 0.01, a 64-dimensional latent space with one 512-unit encoder layer, 90% of droplets for training, and excluded droplets with fewer than 5 UMIs.",
          ## app never sets --expected-cells / --total-droplets-included (app-cellBender.R:161-170)
          "Unless given in cmdOptions, the expected number of cells and the number of droplets included in the analysis were estimated by CellBender's own heuristic from the ranked UMI-count curve.",
          ## consts.py RANDOM_SEED = 1234, applied in run.py (pyro.util.set_rng_seed); the app sets no seed
          "CellBender used its fixed internal random seed (1234); the app set no seed of its own.",
          ## --estimator default mckp (--help); consts.py CELL_PROB_CUTOFF = 0.5; outputs kept at app-cellBender.R:180-187
          "Denoised counts were computed with the MCKP estimator (CellBender default), and two matrices were kept: the full matrix with every input barcode, and the droplets with a posterior cell probability above 0.5 (CellBender's filtered output)."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodCellBender
        name <<- "EzAppCellBender"
        appDefaults <<- rbind(
          cmdOptions = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "for -expected-cells and --total-droplets-included"
          ),
          gpu = ezFrame(
            Type = "numeric",
            DefaultValue = 0,
            Description = "defines the number of gpu to run it with cuda option"
          )
        )
      }
    )
  )

#' Drop ATAC Peaks from a multimodal 10x h5 before CellBender inference.
#'
#' CellBender models ambient RNA (and also antibody/hashtag contamination
#' for CITE-seq / CellPlex inputs), but it is NOT designed for ATAC
#' fragments. On ARC/multiome inputs the ~130-180k Peaks features inflate
#' runtime by ~10x because MCKP posterior estimation scales with
#' features x cells. We drop only the "Peaks" feature_type and keep
#' Gene Expression, Antibody Capture, Multiplexing Capture, etc.
#'
#' Returns the path to the filtered h5, or the original inputFile if no
#' Peaks were present.
dropPeaksFromH5 <- function(inputFile, sampleName) {
  featureTypes <- as.character(rhdf5::h5read(inputFile, "matrix/features/feature_type"))
  typeCounts <- table(featureTypes)

  if (!"Peaks" %in% names(typeCounts)) {
    return(inputFile)
  }

  ezLog(paste0(
    "Multimodal input detected: ",
    paste(names(typeCounts), typeCounts, sep = "=", collapse = ", ")
  ))
  ezLog("Dropping Peaks (ATAC) for CellBender; keeping all other modalities.")

  sce <- DropletUtils::read10xCounts(inputFile, col.names = TRUE)
  keepMask <- rowData(sce)$Type != "Peaks"
  sceKeep <- sce[keepMask, ]

  ezLog(paste0(
    "Keeping ", sum(keepMask), " non-Peaks features, dropping ",
    sum(!keepMask), " Peaks."
  ))

  outH5 <- paste0(sampleName, "_nopeaks.h5")
  if (file.exists(outH5)) file.remove(outH5)

  DropletUtils::write10xCounts(
    outH5,
    counts(sceKeep),
    gene.id = rowData(sceKeep)$ID,
    gene.symbol = rowData(sceKeep)$Symbol,
    type = "HDF5",
    genome = unique(as.character(rhdf5::h5read(inputFile, "matrix/features/genome"))),
    version = "3"
  )

  return(outH5)
}

ezMethodCellBender <- function(input = NA, output = NA, param = NA) {
  require(DropletUtils)

  sampleName = input$getNames()
  setwdNew(sampleName)

  # Initialize cmDir before tryCatch
  cmDir <- NULL

  # Try to get path for single modality case first
  tryCatch(
    {
      if ("UnfilteredCountMatrix" %in% input$colNames) {
        cmDir <- input$getFullPaths("UnfilteredCountMatrix")
      } else {
        # Multi-modal case - construct path from ResultDir
        resultPath <- input$getColumn("ResultDir")
        cmDir <- file.path(
          param$dataRoot,
          resultPath,
          "multi/count/raw_feature_bc_matrix"
        )
        ## file.exists(), NOT exists(): exists() tests for an R OBJECT of that name, so it is
        ## always FALSE for a path string and this branch always overwrote the run-level path
        ## above with one that does not exist for any `cellranger multi` output. Measured:
        ## file.exists(p) TRUE while exists(p) FALSE, which made every Multi CellBender job
        ## fail on a bad path unless the dataset happened to carry UnfilteredCountMatrix.
        if (!file.exists(cmDir)) {
          cmDir <- file.path(
            param$dataRoot,
            resultPath,
            "count/sample_raw_feature_bc_matrix"
          )
        }
      }
    },
    error = function(e) {
      stop(sprintf("Failed to construct path: %s", e$message))
    }
  )

  if (is.null(cmDir)) {
    stop("Failed to get valid path for count matrix")
  }

  inputFile <- paste0(cmDir, '.h5')
  if (!file.exists(inputFile)) {
    warning('RawCountMatrix missing! Creating it instead')

    sce <- read10xCounts(cmDir, col.names = TRUE)

    # Save as h5 file using write10xCounts
    ezLog("Saving to h5 format...")
    inputFile <- paste0(sampleName, '.h5')
    write10xCounts(
      inputFile,
      counts(sce),
      type = "HDF5",
      genome = param$ezRef@refFeatureFile,
      version = "3",
      chemistry = input$getColumn("SCDataOrigin")
    )

    ezLog("Created h5 file: ", inputFile)
  }

  inputFile <- dropPeaksFromH5(inputFile, sampleName)

  cmd <- paste(
    "cellbender remove-background",
    "--input",
    inputFile,
    "--output cellbender.h5"
  )

  if (param$cmdOptions != '') {
    cmd <- paste(cmd, param$cmdOptions)
  }

  if (param$gpu > 0) {
    cmd <- paste(cmd, "--cuda")
  } else {
    cmd <- paste(cmd, '--cpu-threads', param$cores)
  }
  system(cmd)

  ##Post processing for Seurat
  cmd <- paste(
    "ptrepack --complevel 5 cellbender_filtered.h5:/matrix cellbender_filtered_seurat.h5:/matrix"
  )
  system(cmd)
  cmd <- paste(
    "ptrepack --complevel 5 cellbender.h5:/matrix cellbender_raw_seurat.h5:/matrix"
  )
  system(cmd)
  ##Clean Up:
  system(
    'rm ckpt.tar.gz cellbender_posterior.h5 cellbender.h5 cellbender_filtered.h5'
  )
  if (paste0(cmDir, '.h5') != inputFile) {
    system(paste("rm", inputFile))
  }

  return("Success")
}
