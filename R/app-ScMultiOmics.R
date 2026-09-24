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
      methods_facts = function() {
        c(
          ## ezMethodScMultiOmics loads scData.qs2 as is (app-ScMultiOmics.R:137-147); attachUpstreamAnnotations (multiOmicsUtils.R:500-590); pickCellTypeColumn (multiOmicsUtils.R:599-685)
          "RNA normalization, RNA PCA, RNA clusters and cell-type labels were taken unchanged from the upstream ScSeurat object (Azimuth, SingleR fine-label and cellxgene results saved next to it were re-attached); no new cell-type annotation was run, and the cell-type labels shown were the first available of CyteTypeR, cellxgene, Azimuth Pan-Human, Azimuth tissue reference, scType, SingleR and manual labels, in that order.",
          ## detectModalities (multiOmicsUtils.R:55-64); readADTCounts (multiOmicsUtils.R:186-195); processADT zero-count filter (multiOmicsUtils.R:1199-1206)
          "ADT counts were the Antibody Capture features of the CellRanger filtered matrix; features listed as hashtag_ids in the CellRanger Multi configuration were removed (a sample with only hashtag antibodies was treated as RNA-only), and cells with zero total ADT counts were removed from the object.",
          ## processADT ADTnorm branch (multiOmicsUtils.R:1212-1293), CLR branch (:1294-1303), non-finite -> 0 (:1307-1332); ADTnorm/Seurat defaults checked on the R 4.6 system library
          "When adtNorm is ADTnorm, ADTnorm was run separately for each marker with all cells treated as a single batch (so no cross-sample landmark alignment) and exclude_zeroes = TRUE, other ADTnorm settings at their defaults; any marker on which ADTnorm failed kept its CLR value. When adtNorm is CLR, all markers were CLR-normalized per cell (Seurat NormalizeData, margin = 2). Non-finite normalized or scaled values were set to 0.",
          ## processADT ScaleData/RunPCA/RunUMAP (multiOmicsUtils.R:1317-1339)
          "ADT data were scaled and an exact PCA (approx = FALSE) was computed on all ADT features, with npcsADT components capped at the number of ADT features minus one; the ADT UMAP used all of these components with 30 neighbours (Seurat RunUMAP: uwot, cosine metric, seed 42).",
          ## app-ScMultiOmics.R:193-214 (BD ADT); loadBDRhapsody (multiOmicsUtils.R:1009-1019)
          "For BD Rhapsody input (SCDataOrigin = BDRhapsody) the object from the BD pipeline was used, its ADT assay was always CLR-normalized (margin = 2) regardless of adtNorm, and when it had no RNA PCA one was built by LogNormalize, 2000 vst variable features, 30 PCs, a k = 20 shared-nearest-neighbour graph on those 30 PCs and Louvain clustering at a fixed resolution of 0.5.",
          ## processATAC (multiOmicsUtils.R:1120-1156); Signac 1.16.0 defaults checked
          "ATAC peak counts came from the CellRanger ARC Peaks matrix, and the object was restricted to cells present in it; peaks were not filtered (min.cells = 0), Signac FindTopFeatures used min.cutoff q5, TF-IDF used method 1 (scale factor 10000), SVD computed 50 LSI components, and the ATAC UMAP used LSI components 2 to 30 (component 1 excluded).",
          ## getATACAnnotation (multiOmicsUtils.R:1068-1080); processATAC GeneActivity (:1158-1169); Signac GeneActivity/GetGRangesFromEnsDb defaults checked
          "Gene annotation for ATAC was EnsDb.Hsapiens.v86 when refBuild names a human build and EnsDb.Mmusculus.v79 when it names a mouse build (none otherwise, and then no gene activity was computed); gene activity was Signac GeneActivity (fragments over protein-coding gene bodies extended 2 kb upstream), log-normalized with a scale factor equal to the median total gene-activity count per cell.",
          ## runWNN (multiOmicsUtils.R:773-818); FindMultiModalNeighbors / RunUMAP / FindClusters defaults checked on Seurat 5.5.1
          "When runWNN is true and at least two of RNA PCA, ADT PCA and ATAC LSI exist, Seurat FindMultiModalNeighbors (k.nn = 20) combined RNA PCs 1-20 (fixed, independent of upstream npcs), ADT PCs 1 to at most 18 and LSI components 2-30; the WNN UMAP was built from the 20 weighted nearest neighbours (seed 42), and WNN clusters were found on the weighted SNN graph with the SLM algorithm (FindClusters algorithm = 3, random.seed 0) at wnnResolution. These WNN clusters replaced seurat_clusters and are the clusters used throughout the report.",
          ## _scMultiOmics_wnn.Rmd:43-61 (FindAllMarkers); _scMultiOmics_adt.Rmd:214-262 (ADT-RNA correlation); FindAllMarkers defaults checked
          "WNN cluster markers were found with Seurat FindAllMarkers (Wilcoxon test, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) separately on the RNA, ADT and gene-activity assays, keeping markers with unadjusted p below 0.01 (FindAllMarkers default) and reporting Bonferroni-adjusted p-values; ADT-RNA agreement was the Spearman correlation of per-cluster average ADT and matched-gene RNA values, with isotype controls excluded.",
          ## processVDJ (multiOmicsUtils.R:381-411); loadBDContigs (multiOmicsUtils.R:274-289); scRepertoire 2.8.0 defaults checked
          "VDJ contigs were combined with scRepertoire combineTCR or combineBCR keeping cells with missing or extra chains (removeNA, removeMulti and filterMulti all FALSE) and dropping non-productive contigs (scRepertoire default); only one receptor type was attached per run, the TCR when TCR contigs were available and otherwise the BCR, even when vdjChain is both. For BD Rhapsody the dominant-contig AIRR table was used.",
          ## processVDJ combineBCR / clonalCluster / effectiveClone (multiOmicsUtils.R:407-451); scRepertoire 2.8.0 combineBCR and clonalCluster defaults checked
          "BCR clones were merged by combineBCR on IGH CDR3 nucleotide sequences with length-normalized Levenshtein similarity at bcrSimilarityThreshold, requiring the same V and J genes; clone identity in combineExpression then followed cloneCallTCR for both TCR and BCR data. When tcrSimilarityMerge is true, clonalCluster grouped TRB CDR3 amino-acid sequences sharing a V gene at tcrSimilarityThreshold, but if the log warns that no TRB_cluster column was produced, clones were defined by cloneCallTCR instead.",
          ## processVDJ combineExpression (multiOmicsUtils.R:453-468); _scMultiOmics_vdj.Rmd:117, 142, 160
          "Clone sizes were counted within each sample (combineExpression, proportion = FALSE) and binned as Single (1 cell), Small (2-5), Medium (6-20), Large (21-100) and Hyperexpanded (101-500), with cells lacking a clonotype labelled No clonotype; clonal overlap (Morisita-Horn) and clonal homeostasis in the report always used the scRepertoire strict clone definition (cloneCall = strict), regardless of cloneCallTCR."
        )
      },
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
          ),
          cloneCallTCR = ezFrame(
            Type = "character",
            DefaultValue = "strict",
            Description = "How to merge TCR clones in combineExpression: 'strict' (V/J + CDR3 nt) or 'aa' (CDR3 amino acid)."
          ),
          tcrSimilarityMerge = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "Run scRepertoire::clonalCluster on TCR before combineExpression to collapse near-identical CDR3 sequences."
          ),
          tcrSimilarityThreshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.85,
            Description = "Normalized similarity threshold for TCR clonalCluster (0-1). Only used when tcrSimilarityMerge = TRUE."
          ),
          bcrSimilarityThreshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.85,
            Description = "Threshold passed to combineBCR for SHM-aware clone merging."
          )
        )
      }
    )
  )

##' @title Resolve the path to the upstream ScSeurat scData.qs2.
##' @description The ScMultiOmics input dataset gates on `SC Seurat` (see
##'   ScMultiOmicsApp.rb#@required_columns), so the column is guaranteed when
##'   ScSeurat is the upstream. BD Rhapsody bypasses ScSeurat and is handled
##'   separately by the caller.
##' @param input EzDataset.
##' @return Character path to scData.qs2, or "" if not present.
##' @keywords internal
findScDataPath <- function(input) {
  if (!"SC Seurat" %in% input$colNames) return("")
  input$getFullPaths("SC Seurat")
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
      stop("'SC Seurat' column missing or scData.qs2 not found; ",
           "run ScSeurat upstream of ScMultiOmics.")
    }
    if (!nzchar(countMatrixPath)) {
      # ScSeurat propagates CountMatrix in its next_dataset, so this branch
      # only fires for hand-built dataset.tsv files that didn't include it.
      stop("CountMatrix column missing from input dataset; required to locate ",
           "ADT/VDJ/ATAC siblings. Re-run ScSeurat (which propagates ",
           "CountMatrix automatically) or hand-add the column.")
    }
    message("Loading annotated scData from: ", scDataPath)
    obj <- qs2::qs_read(scDataPath, nthreads = max(1L, param$cores %||% 4L))
    # Re-attach sibling cell-type annotations that ScSeurat saves alongside
    # scData.qs2 (Azimuth tissue refs, SingleR per-reference labels,
    # cellxgene label-transfer). ScSeurat doesn't write these onto
    # scData@meta.data directly — they live in aziResults.qs2,
    # singler.results.qs2/.rds, cellxgeneResults.qs2/.rds in the same dir.
    # Without this, ScMultiOmics's pickCellTypeColumn() only sees scType +
    # CyteTypeR + AzimuthPanHuman (the schemes that *do* write directly).
    obj <- attachUpstreamAnnotations(obj, dirname(scDataPath))
  }

  ## 3. Detect modalities -----------------------------------------------------
  if (isBD) {
    bd_assays <- Seurat::Assays(obj)
    bd_meta_cols <- colnames(obj@meta.data)
    # BD encodes VDJ as wide metadata columns on the prebuilt Seurat object
    # (BCR_Heavy_*_Dominant, TCR_Alpha_Gamma_*_Dominant, ...). Detect either
    # chain family by looking for a CDR3 nucleotide column.
    bdHasTCR <- any(grepl("^TCR_.*_CDR3_Nucleotide_Dominant$", bd_meta_cols))
    bdHasBCR <- any(grepl("^BCR_.*_CDR3_Nucleotide_Dominant$", bd_meta_cols))
    mod <- list(hasRNA   = "RNA" %in% bd_assays,
                hasADT   = "ADT" %in% bd_assays,
                hasVDJ_T = bdHasTCR,
                hasVDJ_B = bdHasBCR,
                hasATAC  = any(c("ATAC", "peaks") %in% bd_assays))
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
    adt <- readADTCounts(h5, sampleName = sampleName,
                         objBarcodes = colnames(obj))
    if (!is.null(adt) && ncol(adt) > 0) {
      message("Adding ADT assay: ", nrow(adt), " features x ", ncol(adt), " cells.")
      obj <- processADT(obj, adt,
                        normMethod = param$adtNorm %||% "ADTnorm",
                        npcs = param$npcsADT %||% 18,
                        cores = max(1L, as.integer(param$cores %||% 4L)))
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
    contigsT <- NULL; contigsB <- NULL
    if (isBD) {
      # BD ships an AIRR contig table; scRepertoire's BD format reader
      # parses it once; we filter to T/B and feed processVDJ() pre-loaded
      # data frames.
      bdContigs <- loadBDContigs(bd_root)
      contigsT <- if (wantVDJ_T) bdContigs$T else NULL
      contigsB <- if (wantVDJ_B) bdContigs$B else NULL
      if (wantVDJ_T && (is.null(contigsT) || nrow(contigsT) == 0)) {
        mod$hasVDJ_T <- FALSE; wantVDJ_T <- FALSE
      }
      if (wantVDJ_B && (is.null(contigsB) || nrow(contigsB) == 0)) {
        mod$hasVDJ_B <- FALSE; wantVDJ_B <- FALSE
      }
    } else {
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
    }
    src_T <- if (!is.null(contigsT)) sprintf("BD AIRR (%d contigs)", nrow(contigsT))
             else if (!is.null(vdjTPath)) basename(vdjTPath) else NULL
    src_B <- if (!is.null(contigsB)) sprintf("BD AIRR (%d contigs)", nrow(contigsB))
             else if (!is.null(vdjBPath)) basename(vdjBPath) else NULL
    message("Attaching VDJ clones: ",
            if (wantVDJ_T && !is.null(src_T)) paste0("T(", src_T, ")") else "",
            if (wantVDJ_T && wantVDJ_B) " + " else "",
            if (wantVDJ_B && !is.null(src_B)) paste0("B(", src_B, ")") else "")
    obj <- processVDJ(obj, vdjTPath = vdjTPath, vdjBPath = vdjBPath,
                      contigsT = contigsT, contigsB = contigsB,
                      sampleName = sampleName,
                      cloneCallTCR           = param$cloneCallTCR           %||% "strict",
                      tcrSimilarityMerge     = isTRUE(param$tcrSimilarityMerge),
                      tcrSimilarityThreshold = as.numeric(param$tcrSimilarityThreshold %||% 0.85),
                      bcrSimilarityThreshold = as.numeric(param$bcrSimilarityThreshold %||% 0.85))
    # Persist the raw clones list as TSV alongside the report
    cb <- attr(obj, "vdjCombined")
    if (!is.null(cb)) {
      tsv <- file.path(getwd(), paste0("clones_", sampleName, ".tsv"))
      utils::write.table(cb[[1]], tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  ## 5. WNN integration -------------------------------------------------------
  ranWNN <- FALSE
  # SUSHI passes all params as quoted strings via run_RApp; coerce to numeric
  # so FindMultiModalNeighbors / FindClusters get a real double.
  wnn_resolution <- suppressWarnings(as.numeric(param$wnnResolution %||% 0.5))
  if (is.na(wnn_resolution)) wnn_resolution <- 0.5
  # param$runWNN comes in as 'true'/'false'/'TRUE' string from SUSHI; treat
  # missing as default-on. Same pattern used elsewhere in ezRun for boolean
  # SUSHI params.
  runWNNFlag <- if (is.null(param$runWNN)) TRUE
                else isTRUE(as.logical(param$runWNN))
  if (runWNNFlag) {
    n_dim_mod <- sum(c("pca", "adt.pca", "lsi") %in% names(obj@reductions))
    if (n_dim_mod >= 2L) {
      obj <- runWNN(obj, resolution = wnn_resolution)
      ranWNN <- "wnn.umap" %in% names(obj@reductions)
    }
  }
  mod$ranWNN <- ranWNN
  mod$wnnResolution <- wnn_resolution

  ## 6. Save and render -------------------------------------------------------
  # Derive the exploreSC URL from the gstore-relative output Report path so
  # the Rmd can render a working "Single Cell Explorer" link without the
  # cwd-stripping heuristic (which doesn't recognize SLURM scratch dirs like
  # /scratch/p31662_..._<sample>_temp<pid>/<report_dir>).
  if ("Report" %in% output$colNames && (is.null(param$exploreSCUrl) ||
                                        !nzchar(param$exploreSCUrl %||% ""))) {
    report_rel <- output$getColumn("Report")
    if (length(report_rel) >= 1L && nzchar(report_rel[[1]])) {
      param$exploreSCUrl <- paste0(
        "https://fgcz-shiny.uzh.ch/app/exploreSC/?data=",
        sub("/+$", "", report_rel[[1]]),
        "/scMultiData.qs2"
      )
    }
  }

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
