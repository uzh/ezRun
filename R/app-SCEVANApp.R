###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title SCEVAN CNV Analysis Method
##' @description Main processing function for SCEVAN CNV inference from scRNA-seq data
##' @param input an object of class EzDataset with input metadata
##' @param output an object of class EzDataset with output metadata
##' @param param a list of parameters including:
##'   \itemize{
##'     \item{name: }{Output name for the analysis}
##'     \item{dataRoot: }{Root directory for data files}
##'     \item{cellTypeColumn: }{Metadata column containing cell type annotations}
##'     \item{referenceCellTypes: }{Cell types to use as normal reference}
##'     \item{sampleColumn: }{Metadata column containing sample identifiers}
##'     \item{organism: }{Organism: "human" or "mouse"}
##'     \item{cores: }{Number of cores for parallel processing}
##'     \item{subclones: }{Enable subclone analysis}
##'     \item{betaVega: }{Segmentation granularity parameter}
##'   }
##' @param htmlFile the name of the HTML output file (default: "00index.html")
##' @return Returns "Success" if analysis completes successfully
ezMethodSCEVANApp <- function(input = NA, output = NA, param = NA,
                              htmlFile = "00index.html") {
  # ============================================================================
  # 1. LOAD REQUIRED PACKAGES
  # ============================================================================
  ezLoadPackage('SCEVAN')
  ezLoadPackage('Seurat')
  ezLoadPackage('qs2')
  ezLoadPackage('dplyr')
  ezLoadPackage('tidyr')

  # ============================================================================
  # 2. SET WORKING DIRECTORY
  # ============================================================================
  setwdNew('SCEVAN')

  # ============================================================================
  # 3. GET INPUT METADATA AND VALIDATE
  # ============================================================================
  dataset <- input$meta
  nSamples <- nrow(dataset)

  message("SCEVAN CNV Analysis")
  message("===================")
  message("Number of input files: ", nSamples)

  # Get input paths - handle both column names for compatibility
  # with ScSeuratCombine (SeuratObject) and ScSeurat (SC Seurat) apps
  seuratPaths <- if (input$hasColumn("SeuratObject")) {
    input$getColumn("SeuratObject")
  } else if (input$hasColumn("SC Seurat")) {
    input$getColumn("SC Seurat")
  } else {
    stop("No Seurat object column found. Expected 'SeuratObject' or 'SC Seurat'")
  }

  # Validate organism parameter
  if (!param$organism %in% c("human", "mouse")) {
    stop("Organism must be 'human' or 'mouse', got: ", param$organism)
  }

  # ============================================================================
  # 4. LOAD SEURAT OBJECT
  # ============================================================================
  message("Loading Seurat object...")

  # Get full path to Seurat object
  seuratPath <- file.path(param$dataRoot, seuratPaths[1])

  if (!file.exists(seuratPath)) {
    stop("Seurat object not found: ", seuratPath)
  }

  # Load Seurat object (handles both qs2 and RDS)
  seuratObj <- ezLoadRobj(seuratPath, nthreads = param$cores)

  # Validate it's a Seurat object
  if (!inherits(seuratObj, "Seurat")) {
    stop("Input file is not a Seurat object: ", class(seuratObj))
  }

  message("Loaded Seurat object with ", ncol(seuratObj), " cells and ",
          nrow(seuratObj), " genes")

  # ============================================================================
  # 5. VALIDATE AND PARSE PARAMETERS
  # ============================================================================
  metadata <- seuratObj@meta.data

  # Validate cellTypeColumn exists
  if (!param$cellTypeColumn %in% colnames(metadata)) {
    stop("Cell type column '", param$cellTypeColumn,
         "' not found in metadata. Available columns: ",
         paste(colnames(metadata), collapse = ", "))
  }

  # Validate sampleColumn exists
  if (!param$sampleColumn %in% colnames(metadata)) {
    stop("Sample column '", param$sampleColumn,
         "' not found in metadata. Available columns: ",
         paste(colnames(metadata), collapse = ", "))
  }

  # Parse reference cell types (handle comma or newline separation)
  referenceCellTypes <- param$referenceCellTypes %>%
    strsplit(",|\n") %>%
    unlist() %>%
    trimws() %>%
    .[. != ""]

  if (length(referenceCellTypes) == 0) {
    stop("No reference cell types specified. Please provide at least one reference cell type.")
  }

  message("Reference cell types: ", paste(referenceCellTypes, collapse = ", "))

  # Identify normal cells
  normCells <- rownames(metadata)[metadata[[param$cellTypeColumn]] %in% referenceCellTypes]

  if (length(normCells) == 0) {
    stop("No cells matched the reference cell types. ",
         "Please check that cell type names match exactly. ",
         "Available cell types: ",
         paste(unique(metadata[[param$cellTypeColumn]]), collapse = ", "))
  }

  normPct <- 100 * length(normCells) / ncol(seuratObj)
  message("Identified ", length(normCells), " normal cells (",
          round(normPct, 2), "% of total)")

  if (normPct < 5) {
    warning("Only ", round(normPct, 2),
            "% of cells are reference cells. This may affect CNV calling quality. ",
            "Consider adding more reference cell types.")
  }

  # ============================================================================
  # 6. PREPARE SAMPLE INFORMATION
  # ============================================================================
  # Get unique samples
  samples <- unique(metadata[[param$sampleColumn]])
  nSamplesInData <- length(samples)

  message("Found ", nSamplesInData, " unique samples in the data")

  if (nSamplesInData == 0) {
    stop("No samples found in column '", param$sampleColumn, "'")
  }

  # ============================================================================
  # 7. RUN SCEVAN FOR EACH SAMPLE
  # ============================================================================
  message("\n=== Running SCEVAN pipelineCNA ===\n")

  resultsList <- list()
  failedSamples <- c()

  for (i in seq_along(samples)) {
    sampleName <- samples[i]
    message("\n[", i, "/", nSamplesInData, "] Processing sample: ", sampleName)

    # Subset Seurat object for this sample
    cellsToKeep <- colnames(seuratObj)[metadata[[param$sampleColumn]] == sampleName]
    message("  Number of cells: ", length(cellsToKeep))

    seuratSub <- subset(seuratObj, cells = cellsToKeep)

    # Extract count matrix (handle Seurat v5 layers)
    if (length(seuratSub@assays$RNA@layers) > 0 &&
        "counts" %in% names(seuratSub@assays$RNA@layers)) {
      # Seurat v5 with layers
      countMtx <- as.matrix(seuratSub@assays$RNA@layers$counts)
      rownames(countMtx) <- rownames(seuratSub@assays$RNA)
      colnames(countMtx) <- colnames(seuratSub)
    } else if (!is.null(seuratSub@assays$RNA@counts) &&
               length(seuratSub@assays$RNA@counts) > 0) {
      # Seurat v3/v4 with counts slot
      countMtx <- as.matrix(GetAssayData(seuratSub, slot = "counts"))
    } else {
      warning("  No counts found for sample ", sampleName, ", skipping")
      failedSamples <- c(failedSamples, sampleName)
      next
    }

    # Identify normal cells in this sample
    currentNormCells <- intersect(normCells, colnames(countMtx))
    message("  Normal cells in this sample: ", length(currentNormCells))

    if (length(currentNormCells) < 10) {
      warning("  Only ", length(currentNormCells),
              " normal cells in sample ", sampleName,
              ". CNV calling may be unreliable. Skipping.")
      failedSamples <- c(failedSamples, sampleName)
      next
    }

    # Run SCEVAN pipelineCNA
    message("  Running pipelineCNA...")
    scevanResult <- tryCatch({
      pipelineCNA(
        count_mtx = countMtx,
        sample = sampleName,
        norm_cell = currentNormCells,
        par_cores = param$cores,
        SUBCLONES = param$subclones,
        beta_vega = param$betaVega,
        ClonalCN = TRUE,
        plotTree = FALSE,
        organism = param$organism
      )
    }, error = function(e) {
      warning("  SCEVAN failed for sample ", sampleName, ": ", e$message)
      return(NULL)
    })

    if (!is.null(scevanResult)) {
      resultsList[[sampleName]] <- scevanResult
      message("  ✓ Sample ", sampleName, " completed successfully")
    } else {
      failedSamples <- c(failedSamples, sampleName)
    }
  }

  # ============================================================================
  # 8. CHECK RESULTS
  # ============================================================================
  nSuccessful <- length(resultsList)
  nFailed <- length(failedSamples)

  message("\n=== SCEVAN Processing Summary ===")
  message("Successful: ", nSuccessful, "/", nSamplesInData)
  if (nFailed > 0) {
    message("Failed: ", nFailed, " (", paste(failedSamples, collapse = ", "), ")")
  }

  if (nSuccessful == 0) {
    stop("SCEVAN failed for all samples. Check error messages above.")
  }

  # ============================================================================
  # 9. RUN MULTI-SAMPLE COMPARISON (if >1 sample)
  # ============================================================================
  multiSampleResults <- NULL

  if (nSuccessful > 1) {
    message("\n=== Running multi-sample comparison ===\n")

    # Prepare count matrices list
    listCountMtx <- list()
    listNormCells <- list()

    for (sampleName in names(resultsList)) {
      cellsToKeep <- colnames(seuratObj)[metadata[[param$sampleColumn]] == sampleName]
      seuratSub <- subset(seuratObj, cells = cellsToKeep)

      # Extract count matrix
      if (length(seuratSub@assays$RNA@layers) > 0 &&
          "counts" %in% names(seuratSub@assays$RNA@layers)) {
        countMtx <- as.matrix(seuratSub@assays$RNA@layers$counts)
        rownames(countMtx) <- rownames(seuratSub@assays$RNA)
        colnames(countMtx) <- colnames(seuratSub)
      } else {
        countMtx <- as.matrix(GetAssayData(seuratSub, slot = "counts"))
      }

      listCountMtx[[sampleName]] <- countMtx
      listNormCells[[sampleName]] <- intersect(normCells, colnames(countMtx))
    }

    # Run multi-sample comparison
    multiSampleResults <- tryCatch({
      multiSampleComparisonClonalCN(
        listCountMtx = listCountMtx,
        listNormCells = listNormCells,
        analysisName = param$name,
        organism = param$organism,
        par_cores = param$cores
      )
    }, error = function(e) {
      warning("Multi-sample comparison failed: ", e$message)
      return(NULL)
    })

    if (!is.null(multiSampleResults)) {
      message("✓ Multi-sample comparison completed")
    }
  }

  # ============================================================================
  # 10. SAVE RESULTS
  # ============================================================================
  message("\n=== Saving results ===\n")

  # Save SCEVAN results
  qs2::qs_save(resultsList, "scevan_results_list.qs2", nthreads = param$cores)
  message("✓ Saved SCEVAN results list")

  if (!is.null(multiSampleResults)) {
    qs2::qs_save(multiSampleResults, "scevan_multi_sample_results.qs2",
                 nthreads = param$cores)
    message("✓ Saved multi-sample comparison results")
  }

  # Save Seurat object for report
  qs2::qs_save(seuratObj, "seurat_object.qs2", nthreads = param$cores)
  message("✓ Saved Seurat object")

  # Save parameters for Rmd
  paramForRmd <- param
  paramForRmd$referenceCellTypes <- referenceCellTypes
  paramForRmd$normCells <- normCells
  paramForRmd$samples <- samples
  paramForRmd$nSuccessful <- nSuccessful
  paramForRmd$failedSamples <- failedSamples
  write_rds(paramForRmd, 'param.rds')
  message("✓ Saved parameters")

  # ============================================================================
  # 11. GENERATE HTML REPORT
  # ============================================================================
  message("\n=== Generating HTML report ===\n")

  reportTitle <- paste('SCEVAN CNV Analysis -', param$name)

  makeRmdReport(
    rmdFile = "SCEVAN.Rmd",
    reportTitle = reportTitle
  )

  # ============================================================================
  # 12. RETURN SUCCESS
  # ============================================================================
  message("\n✓ SCEVAN analysis completed successfully!\n")
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSCEVANApp(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run SCEVAN CNV analysis
##' @section Functions:
##' \describe{
##'   \item{initialize}{Initializes the application using its specific defaults.}
##' }
##' @section Slots:
##' \describe{
##'   \item{runMethod}{The method to run (ezMethodSCEVANApp).}
##'   \item{name}{The name of the application (EzAppSCEVANApp).}
##'   \item{appDefaults}{Default parameters for the application.}
##' }
EzAppSCEVANApp <-
  setRefClass("EzAppSCEVANApp",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSCEVANApp
        name <<- "EzAppSCEVANApp"

        # Define application-specific default parameters
        appDefaults <<- rbind(
          # Cell type column containing annotations
          cellTypeColumn = ezFrame(
            Type = "character",
            DefaultValue = "predicted.id",
            Description = "Metadata column containing cell type annotations (e.g., predicted.id, cellTypeIntegrated, Azimuth.celltype.l2)"
          ),

          # Reference cell types (normal cells)
          referenceCellTypes = ezFrame(
            Type = "character",
            DefaultValue = "CD4-positive, alpha-beta T cell,CD8-positive, alpha-beta T cell,B cell,natural killer cell,plasma cell",
            Description = "Comma or newline separated list of cell types to use as normal reference for CNV calling"
          ),

          # Sample column
          sampleColumn = ezFrame(
            Type = "character",
            DefaultValue = "Sample",
            Description = "Metadata column containing sample identifiers"
          ),

          # Organism
          organism = ezFrame(
            Type = "charVector",
            DefaultValue = "human",
            Description = "Organism: human or mouse",
            Values = c("human", "mouse")
          ),

          # Number of cores
          cores = ezFrame(
            Type = "numeric",
            DefaultValue = 48,
            Description = "Number of cores for parallel processing (SCEVAN is computationally intensive)"
          ),

          # Enable subclones analysis
          subclones = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Enable subclone identification and analysis"
          ),

          # Beta vega parameter
          betaVega = ezFrame(
            Type = "numeric",
            DefaultValue = 0.5,
            Description = "Segmentation granularity parameter (0-1). Lower values = more segments."
          )
        )
      }
    )
  )
