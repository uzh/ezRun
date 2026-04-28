###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title EzAppScMultiOmics
##' @description Multi-omics extension layered on top of an annotated ScSeurat
##'   object. Reads `scData.qs2` from a previous ScSeurat run, attaches
##'   modalities discovered next to the original `CountMatrix` (ADT today;
##'   VDJ + ATAC + WNN in later phases) and writes `scMultiData.qs2` plus an
##'   HTML report.
##' @export
EzAppScMultiOmics <-
  setRefClass(
    "EzAppScMultiOmics",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodScMultiOmics
        name <<- "EzAppScMultiOmics"
        appDefaults <<- rbind(
          runWNN = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "Run WNN when 2+ dimensional modalities are present (Phase 3)."
          ),
          adtNorm = ezFrame(
            Type = "character",
            DefaultValue = "ADTnorm",
            Description = "ADT normalization method: 'ADTnorm' or 'CLR'."
          ),
          vdjChain = ezFrame(
            Type = "character",
            DefaultValue = "auto",
            Description = "VDJ chain selection: TCR | BCR | both | auto (Phase 2)."
          ),
          npcsADT = ezFrame(
            Type = "numeric",
            DefaultValue = 18,
            Description = "Number of PCs for the ADT assay."
          )
        )
      }
    )
  )

##' @title Resolve the path to the upstream ScSeurat scData.qs2.
##' @description ScSeurat exposes its output via different column names depending
##'   on the SUSHI dataset shape. Try the most specific link first, then fall
##'   back to the report directory.
##' @param input EzDataset.
##' @return Character path to scData.qs2 (existence not guaranteed).
##' @keywords internal
findScDataPath <- function(input) {
  cols <- input$colNames
  if ("SC Seurat" %in% cols) {
    p <- input$getFullPaths("SC Seurat")
    if (file.exists(p)) return(p)
  }
  for (cand in c("SC Cluster Report", "Report", "Static Report")) {
    if (cand %in% cols) {
      p <- input$getFullPaths(cand)
      d <- if (dir.exists(p)) p else dirname(p)
      hit <- file.path(d, "scData.qs2")
      if (file.exists(hit)) return(hit)
    }
  }
  ""
}

##' @title ezMethodScMultiOmics runtime.
##' @description Loads an annotated `scData.qs2`, detects modalities next to the
##'   original `CountMatrix`, attaches an ADT assay (Phase 1), saves
##'   `scMultiData.qs2`, and renders the multi-omics HTML report.
##' @param input EzDataset (rows of the dataset.tsv).
##' @param output EzDataset (output rows).
##' @param param Job parameters.
##' @param htmlFile Output HTML basename.
##' @return "Success" on success.
##' @export
ezMethodScMultiOmics <- function(input = NA, output = NA, param = NA,
                                 htmlFile = "00index.html") {
  library(Seurat)
  library(qs2)

  cwd <- getwd()
  reportName <- if ("Report" %in% output$colNames) {
    basename(output$getColumn("Report"))
  } else {
    paste0(input$getNames(), "_ScMultiOmicsReport")
  }
  setwdNew(reportName)
  on.exit(setwd(cwd), add = TRUE)

  ## 1. Locate inputs ---------------------------------------------------------
  scDataPath <- findScDataPath(input)
  if (!nzchar(scDataPath) || !file.exists(scDataPath)) {
    stop("scData.qs2 not found in input dataset; run ScSeurat first.")
  }
  if (!"CountMatrix" %in% input$colNames) {
    stop("CountMatrix column missing from input dataset; required to locate ADT/VDJ/ATAC siblings.")
  }
  countMatrixPath <- input$getFullPaths("CountMatrix")

  ## 2. Load annotated RNA object --------------------------------------------
  message("Loading annotated scData from: ", scDataPath)
  obj <- qs2::qs_read(scDataPath, nthreads = max(1L, param$cores %||% 4L))

  ## 3. Detect modalities -----------------------------------------------------
  mod <- detectModalities(countMatrixPath)
  message("Detected modalities: ",
          paste(names(mod)[unlist(mod)], collapse = ", "))
  saveRDS(mod, "modalities.rds")

  ## 4. Per-modality processing (Phase 1: ADT only) --------------------------
  if (isTRUE(mod$hasADT)) {
    h5 <- findFilteredH5(countMatrixPath)
    sampleName <- input$getNames()
    adt <- readADTCounts(h5, sampleName = sampleName)
    if (!is.null(adt) && ncol(adt) > 0) {
      message("Adding ADT assay: ", nrow(adt), " features x ", ncol(adt), " cells.")
      obj <- processADT(obj, adt,
                        normMethod = param$adtNorm %||% "ADTnorm",
                        npcs = param$npcsADT %||% 18)
    } else {
      warning("hasADT was TRUE but readADTCounts returned no data; skipping.")
      mod$hasADT <- FALSE
    }
  }

  ## 5. Save and render -------------------------------------------------------
  makeRmdReport(
    scMultiData = obj,
    param = param,
    modalities = mod,
    rmdFile = "ScMultiOmics.Rmd",
    reportTitle = paste0("ScMultiOmics: ", input$getNames()),
    use.qs2 = TRUE,
    nthreads = max(1L, param$cores %||% 4L)
  )

  return("Success")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
