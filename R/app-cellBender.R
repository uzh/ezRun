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
        if (!exists(cmDir)) {
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
