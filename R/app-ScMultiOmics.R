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
  scDataPath      <- findScDataPath(input)
  countMatrixPath <- if ("CountMatrix" %in% input$colNames)
    input$getFullPaths("CountMatrix") else ""
  bdResultDir     <- if ("BDRhapsodyPath" %in% input$colNames)
    input$getFullPaths("BDRhapsodyPath") else ""
  scDataOrigin    <- if ("SCDataOrigin" %in% input$colNames)
    input$getColumn("SCDataOrigin") else ""

  isBD <- nzchar(bdResultDir) || identical(unname(scDataOrigin), "BDRhapsody")

  ## 2. Load annotated RNA object --------------------------------------------
  if (isBD) {
    bd_root <- if (nzchar(bdResultDir)) bdResultDir else
      if (nzchar(countMatrixPath)) dirname(countMatrixPath) else
      stop("BD Rhapsody mode requires BDRhapsodyPath or CountMatrix.")
    message("Loading BD Rhapsody Seurat from: ", bd_root)
    obj <- loadBDRhapsody(bd_root, sampleName = input$getNames())
    if (is.null(obj)) stop("BD Rhapsody Seurat not found in ", bd_root)
  } else {
    if (!nzchar(scDataPath) || !file.exists(scDataPath)) {
      stop("scData.qs2 not found in input dataset; run ScSeurat first.")
    }
    if (!nzchar(countMatrixPath)) {
      stop("CountMatrix column missing from input dataset; required to locate ADT/VDJ/ATAC siblings.")
    }
    message("Loading annotated scData from: ", scDataPath)
    obj <- qs2::qs_read(scDataPath, nthreads = max(1L, param$cores %||% 4L))
  }

  ## 3. Detect modalities -----------------------------------------------------
  if (isBD) {
    mod <- list(hasRNA  = TRUE,
                hasADT  = "ADT" %in% Seurat::Assays(obj),
                hasVDJ_T = "VDJTPath" %in% input$colNames,
                hasVDJ_B = "VDJBPath" %in% input$colNames,
                hasATAC  = FALSE)
  } else {
    mod <- detectModalities(countMatrixPath)
  }
  # Standalone VDJ columns override auto-discovery
  if ("VDJTPath" %in% input$colNames &&
      nzchar(input$getColumn("VDJTPath")[[1]])) mod$hasVDJ_T <- TRUE
  if ("VDJBPath" %in% input$colNames &&
      nzchar(input$getColumn("VDJBPath")[[1]])) mod$hasVDJ_B <- TRUE
  message("Detected modalities: ",
          paste(names(mod)[unlist(mod)], collapse = ", "))
  saveRDS(mod, "modalities.rds")

  sampleName <- input$getNames()

  ## 4. Per-modality processing ----------------------------------------------
  if (isTRUE(mod$hasADT) && !isBD) {
    h5 <- findFilteredH5(countMatrixPath)
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
  } else if (isTRUE(mod$hasADT) && isBD) {
    # BD Rhapsody Seurat already carries the ADT assay; build a quick PCA + UMAP
    # so the ADT tab and (later) WNN have something to consume.
    if (!"adt.pca" %in% names(obj@reductions)) {
      message("Computing ADT PCA + UMAP on BD Rhapsody ADT assay.")
      Seurat::DefaultAssay(obj) <- "ADT"
      obj <- Seurat::NormalizeData(obj, assay = "ADT",
                                   normalization.method = "CLR", margin = 2,
                                   verbose = FALSE)
      obj <- Seurat::ScaleData(obj, assay = "ADT", verbose = FALSE)
      npcs <- min(param$npcsADT %||% 18, nrow(obj[["ADT"]]) - 1)
      obj <- Seurat::RunPCA(obj, assay = "ADT", reduction.name = "adt.pca",
                            npcs = npcs, features = rownames(obj[["ADT"]]),
                            approx = FALSE, verbose = FALSE)
      n_neighbors <- min(30L, max(2L, ncol(obj) - 1L))
      obj <- Seurat::RunUMAP(obj, assay = "ADT", reduction = "adt.pca",
                             dims = seq_len(npcs), reduction.name = "adt.umap",
                             n.neighbors = n_neighbors, verbose = FALSE)
      Seurat::DefaultAssay(obj) <- "RNA"
    }
  }

  if (isTRUE(mod$hasATAC)) {
    atac_files <- findATACFiles(countMatrixPath)
    if (!is.null(atac_files)) {
      message("Adding ATAC assay from: ", basename(atac_files$fragments))
      refBuild <- param$refBuild
      if (is.null(refBuild) && "refBuild" %in% input$colNames) {
        refBuild <- input$getColumn("refBuild")[[1]]
      }
      obj <- processATAC(obj,
                         fragmentsPath = atac_files$fragments,
                         peaksPath = atac_files$peaks,
                         refBuild = refBuild,
                         sampleName = sampleName)
    } else {
      warning("hasATAC was TRUE but ATAC files not found; skipping.")
      mod$hasATAC <- FALSE
    }
  }

  vdjChain <- param$vdjChain %||% "auto"
  wantVDJ_T <- (vdjChain %in% c("auto", "TCR", "both")) && isTRUE(mod$hasVDJ_T)
  wantVDJ_B <- (vdjChain %in% c("auto", "BCR", "both")) && isTRUE(mod$hasVDJ_B)
  if (wantVDJ_T || wantVDJ_B) {
    # Standalone-VDJ columns take priority over auto-discovery.
    vdjTPath <- NULL; vdjBPath <- NULL
    if (wantVDJ_T) {
      vdjTPath <- if ("VDJTPath" %in% input$colNames &&
                      nzchar(input$getColumn("VDJTPath")[[1]]))
        input$getFullPaths("VDJTPath") else
        findVDJContigCsv(countMatrixPath, "T")
      if (is.character(vdjTPath) && length(vdjTPath) > 0 && dir.exists(vdjTPath)) {
        vdjTPath <- file.path(vdjTPath, "filtered_contig_annotations.csv")
      }
    }
    if (wantVDJ_B) {
      vdjBPath <- if ("VDJBPath" %in% input$colNames &&
                      nzchar(input$getColumn("VDJBPath")[[1]]))
        input$getFullPaths("VDJBPath") else
        findVDJContigCsv(countMatrixPath, "B")
      if (is.character(vdjBPath) && length(vdjBPath) > 0 && dir.exists(vdjBPath)) {
        vdjBPath <- file.path(vdjBPath, "filtered_contig_annotations.csv")
      }
    }
    message("Attaching VDJ clones: ",
            if (wantVDJ_T) paste("T (", basename(vdjTPath), ")", sep = "") else "",
            if (wantVDJ_T && wantVDJ_B) " + " else "",
            if (wantVDJ_B) paste("B (", basename(vdjBPath), ")", sep = "") else "")
    obj <- processVDJ(obj, vdjTPath = vdjTPath, vdjBPath = vdjBPath,
                      sampleName = sampleName)
    # Persist the raw clones list as TSV alongside the report
    cb <- attr(obj, "vdjCombined")
    if (!is.null(cb)) {
      tsv <- file.path(getwd(), paste0("clones_", sampleName, ".tsv"))
      utils::write.table(cb[[1]], tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  ## 5. WNN integration -------------------------------------------------------
  ranWNN <- FALSE
  wnn_resolution <- param$wnnResolution %||% 0.5
  if (isTRUE(param$runWNN %||% TRUE)) {
    n_dim_mod <- sum(c("pca", "adt.pca", "lsi") %in% names(obj@reductions))
    if (n_dim_mod >= 2L) {
      obj <- runWNN(obj, resolution = wnn_resolution)
      ranWNN <- "wnn.umap" %in% names(obj@reductions)
    }
  }
  mod$ranWNN <- ranWNN
  mod$wnnResolution <- wnn_resolution

  ## 6. Save and render -------------------------------------------------------
  makeRmdReport(
    scMultiData = obj,
    param = param,
    output = output,
    modalities = mod,
    rmdFile = "ScMultiOmics.Rmd",
    reportTitle = paste0("ScMultiOmics: ", input$getNames()),
    use.qs2 = TRUE,
    nthreads = max(1L, param$cores %||% 4L)
  )

  # Also expose the object as scData.qs2 so downstream apps that read the
  # standard ScSeurat output (e.g. ScSeuratCombine) work without a separate
  # branch. SUSHI's next_dataset advertises the file via 'SC Seurat [Link]'.
  if (file.exists("scMultiData.qs2") && !file.exists("scData.qs2")) {
    file.symlink("scMultiData.qs2", "scData.qs2")
  }

  return("Success")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
