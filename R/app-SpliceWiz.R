###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## SpliceWiz: differential alternative splicing (SE / MXE / A5SS / A3SS / IR /
## AFE / ALE) from STAR BAMs, for a two-group comparison.

##' @title Prepare (build once, cache) a SpliceWiz reference
##' @description Builds a SpliceWiz reference from the genome FASTA + GTF that the
##'   \code{EzRef} points at, and caches it next to the GTF keyed by the reference
##'   build (and the transcript-type filter) -- so it is shared across projects and
##'   built only once. Mirrors the lock/sentinel protocol of
##'   \code{\link{getSTARReference}}. If \code{param$transcriptTypes} restricts the
##'   biotypes, a filtered GTF is written into \code{workDir} and used for the build,
##'   and the cache directory is suffixed by the selection so different filters do
##'   not collide. If the reference directory beside the GTF is not writable, the
##'   reference is built into \code{workDir} for this run instead.
##' @param param the SUSHI/ezRun parameter list; \code{param$ezRef} must resolve
##'   \code{refFeatureFile} (GTF) and \code{refFastaFile} (genome FASTA).
##' @param workDir a scratch directory (never delivered to gstore) for the filtered
##'   GTF and, if needed, the reference itself.
##' @return the path to a complete SpliceWiz reference directory.
prepareSpliceWizRef <- function(param, workDir = getwd()) {
  gtfFile <- param$ezRef["refFeatureFile"]
  fastaFile <- param$ezRef["refFastaFile"]
  if (!ezIsSpecified(gtfFile)) {
    stop("refFeatureFile not defined")
  }

  ## -- optional transcript-type prefilter (as in FeatureCounts/DESeq2) ---------
  ## Restrict the biotypes fed into event discovery. A filtered GTF is written to
  ## the scratch workDir; the cache dir is suffixed by the (sorted) selection so a
  ## protein_coding reference and a full reference coexist beside the shared GTF.
  buildGtf <- gtfFile
  refSuffix <- "_SpliceWizRef"
  if (ezIsSpecified(param$transcriptTypes)) {
    ## SUSHI collapses a multi_selection param to a comma-joined string; split it.
    types <- sort(unique(trimws(unlist(strsplit(as.character(param$transcriptTypes), ",")))))
    annoFile <- param$ezRef["refAnnotationFile"]
    if (ezIsSpecified(annoFile) && file.exists(annoFile)) {
      tag <- gsub("[^A-Za-z0-9]+", "-", paste(types, collapse = "-"))
      if (nchar(tag) > 40) tag <- substr(digest::digest(types), 1, 12)
      refSuffix <- paste0("_SpliceWizRef_", tag)
      buildGtf <- file.path(workDir, "genes_filtered.gtf")
      if (!file.exists(buildGtf)) {
        rtracklayer::export(gtfByTxTypes(param, types), buildGtf)
      }
    } else {
      ezLog("SpliceWiz: refAnnotationFile missing -- ignoring transcriptTypes, using full GTF")
    }
  }

  ## cache dir beside the GTF, keyed by the reference build (as getSTARReference does)
  refDir <- sub("\\.gtf$", refSuffix, gtfFile)
  ## the reference is complete iff SpliceWiz.ref.gz sits inside it (per buildRef docs)
  sentinel <- function(d) file.exists(file.path(d, "SpliceWiz.ref.gz"))

  ## if the shared location is not writable, build into scratch for this run
  if (!sentinel(refDir) && file.access(dirname(refDir), mode = 2) != 0) {
    refDir <- file.path(workDir, basename(refDir))
  }

  ## random sleep to avoid parallel ref building
  Sys.sleep(runif(1, max = 20))

  lockFile <- paste0(refDir, ".lock")
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ## somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (i >= INDEX_BUILD_TIMEOUT) {
    stop(paste(
      "SpliceWiz reference building still in progress after",
      INDEX_BUILD_TIMEOUT,
      "min"
    ))
  }
  if (sentinel(refDir)) {
    ## already built and complete
    return(refDir)
  }

  ## no lock file and no reference, so we build it
  ezWrite(Sys.info(), con = lockFile)
  on.exit(if (file.exists(lockFile)) file.remove(lockFile), add = TRUE)
  job <- ezJobStart("SpliceWiz buildRef")
  ## genome_type left empty: build purely from the provided FASTA + GTF, so any
  ## FGCZ reference works (no reliance on hard-coded hg38/mm10 resources).
  buildRef(
    reference_path = refDir,
    fasta = fastaFile,
    gtf = buildGtf,
    genome_type = "",
    verbose = FALSE
  )
  ezWriteElapsed(job, "done")
  if (!sentinel(refDir)) {
    stop("SpliceWiz buildRef did not produce SpliceWiz.ref.gz in ", refDir)
  }
  file.remove(lockFile)
  return(refDir)
}

##' @title Add a gene_name column to a SpliceWiz ASE result
##' @description The ASE result table carries no gene symbol. Join through
##'   \code{gene_id} (\code{rowData(se)} -> \code{ref(se)$geneList}); fall back to
##'   parsing the symbol out of \code{EventName} for any residual NAs.
##' @return the \code{res} data.frame with a \code{gene_name} column.
addGeneName <- function(res, se) {
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  gl <- ref(se)$geneList
  res$gene_name <- NA_character_
  if (!is.null(gl) && all(c("gene_id", "gene_name") %in% colnames(gl)) &&
    "gene_id" %in% colnames(rd) && "EventName" %in% colnames(rd)) {
    gmap <- setNames(gl$gene_name, gl$gene_id)
    gid <- rd$gene_id[match(res$EventName, rd$EventName)]
    res$gene_name <- unname(gmap[gid])
  }
  ## fallback: symbol is the token before the first '-' or '/' (after the TYPE:)
  na <- is.na(res$gene_name) | res$gene_name == ""
  if (any(na)) {
    res$gene_name[na] <- vapply(
      res$EventName[na],
      function(e) sub("[-/].*$", "", sub(";.*$", "", sub("^[A-Za-z0-9]+:", "", e))),
      character(1)
    )
  }
  res
}

ezMethodSpliceWiz <- function(input = NA, output = NA, param = NA) {
  library(SpliceWiz)
  library(SummarizedExperiment)

  cwd <- getwd() # the SUSHI $SCRATCH_DIR; only the report dir is shipped to gstore
  on.exit(setwd(cwd), add = TRUE)
  reportName <- basename(output$getColumn("Report"))
  stopifnot(param$sampleGroup != param$refGroup)

  ## -- two-group resolution (app-JunctionSeq.R pattern) ----------------------
  param$grouping <- sub("\\ \\[.*", "", param$grouping) # strip Ruby [Factor] suffix
  condition <- make.names(input$getColumn(param$grouping))
  names(condition) <- input$getNames()
  param$sampleGroup <- make.names(param$sampleGroup)
  param$refGroup <- make.names(param$refGroup)
  samplesUse <- condition %in% c(param$sampleGroup, param$refGroup)
  input <- input$subset(samplesUse)
  condition <- condition[samplesUse]

  samples <- input$getNames()
  bamFiles <- input$getFullPaths("BAM")

  ## -- scratch working dir: ALL heavy intermediates live here, NOT in the report
  ## dir. SUSHI ships the whole report dir to gstore recursively and rm -rf's the
  ## rest of the scratch dir, so anything here is cleaned and never delivered. --
  workDir <- file.path(cwd, "sw_work")
  dir.create(workDir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workDir, recursive = TRUE, force = TRUE), add = TRUE)

  ## -- stage BAMs (+ .bai) into workDir (gstore automount can time out mid-job) -
  ## getBamLocally copies to getwd() and returns a bare basename, so run it from
  ## workDir and rebuild absolute paths.
  setwd(workDir)
  localBams <- vapply(bamFiles, getBamLocally, character(1))
  setwd(cwd)
  copied <- localBams != bamFiles
  localBams[copied] <- file.path(workDir, basename(localBams[copied]))

  ## -- reference (build once, cache; optional transcriptTypes prefilter) -------
  refDir <- prepareSpliceWizRef(param, workDir)

  nCores <- as.integer(param$cores)

  ## -- per-sample processBAM --------------------------------------------------
  pbDir <- file.path(workDir, "sw_pb")
  processBAM(
    bamfiles = unname(localBams),
    sample_names = samples,
    reference_path = refDir,
    output_path = pbDir,
    n_threads = nCores,
    overwrite = TRUE,
    verbose = FALSE
  )

  ## -- collate into an NxtSE --------------------------------------------------
  expr <- findSpliceWizOutput(pbDir)
  nxtDir <- file.path(workDir, "sw_nxt")
  ## 'both' (unstranded) libraries must be collated strand-agnostically; for
  ## stranded libraries SpliceWiz auto-detects the direction from the data.
  forceStrandAgnostic <- ezIsSpecified(param$strandMode) &&
    param$strandMode == "both"
  collateData(
    Experiment = expr,
    reference_path = refDir,
    output_path = nxtDir,
    IRMode = "SpliceOver",
    forceStrandAgnostic = forceStrandAgnostic,
    n_threads = nCores,
    overwrite = TRUE
  )

  ## -- build filtered SE with the condition factor in colData -----------------
  ## realize = TRUE pulls the HDF5-backed assays into memory so the NxtSE can be
  ## serialized (qs2) for the report and the saved object.
  se <- makeSE(nxtDir, realize = TRUE, verbose = FALSE)
  colData(se)$condition <- factor(
    condition[colnames(se)],
    levels = c(param$refGroup, param$sampleGroup)
  )
  ## applyFilters() returns a logical row mask, not an SE -- subset with it.
  se <- se[applyFilters(se, getDefaultFilters()), ]

  ## -- link COV files for plotCoverage ----------------------------------------
  ## makeSE does not auto-link the per-sample .cov files, and SpliceWiz stores
  ## cov paths RELATIVE to sourcePath via .make_path_relative(), which does a
  ## buggy *character-level* common-prefix strip -- so a cov dir that shares any
  ## leading characters with the collate dir (e.g. sw_pb vs sw_nxt) is mangled
  ## (-> "../pb/..."). Copying the .cov files INTO the collate dir (sourcePath)
  ## makes the stored relative path a bare basename, sidestepping the bug, and
  ## keeps them beside the NxtSE so plotCoverage resolves them in the separate
  ## report-render process (which runs while workDir still exists). Symlinks do
  ## NOT work: the setter normalizePath()s them back to the original sibling dir.
  covSrc <- file.path(pbDir, paste0(colnames(se), ".cov"))
  covDst <- file.path(nxtDir, basename(covSrc))
  haveCov <- file.exists(covSrc)
  file.copy(covSrc[haveCov], covDst[haveCov], overwrite = TRUE)
  if (all(file.exists(covDst))) {
    covfile(se) <- covDst
  }

  ## -- differential ASE test (engine chosen by param$aseMethod) ---------------
  aseFun <- switch(
    param$aseMethod,
    limma = ASE_limma,
    DESeq = ASE_DESeq,
    edgeR = ASE_edgeR,
    ASE_limma
  )
  aseArgs <- list(
    se = se,
    test_factor = "condition",
    test_nom = param$sampleGroup,
    test_denom = param$refGroup,
    IRmode = ifelse(ezIsSpecified(param$IRmode), param$IRmode, "all")
  )
  ## optional secondary co-variate (~ grouping2 + condition), consistent with the
  ## DESeq2/EdgeR apps' grouping2 parameter.
  if (ezIsSpecified(param$grouping2)) {
    grouping2Col <- sub("\\ \\[.*", "", param$grouping2)
    if (input$hasColumn(grouping2Col)) {
      colData(se)$grouping2 <- factor(make.names(input$getColumn(grouping2Col)[colnames(se)]))
      aseArgs$se <- se
      aseArgs$batch1 <- "grouping2"
    }
  }
  res <- as.data.frame(do.call(aseFun, aseArgs))

  ## -- standardize the FDR / p-value columns across engines -------------------
  fdrCol <- switch(param$aseMethod,
    limma = "adj.P.Val", edgeR = "FDR", DESeq = "padj", "adj.P.Val")
  pvalCol <- switch(param$aseMethod,
    limma = "P.Value", edgeR = "PValue", DESeq = "pvalue", "P.Value")
  res$FDR <- res[[fdrCol]]
  res$PValue <- res[[pvalCol]]

  ## -- gene symbol (ASE table has none) ---------------------------------------
  res <- addGeneName(res, se)
  ## surface gene_name near the front of the written tables
  res <- res[, c("gene_name", setdiff(colnames(res), "gene_name")), drop = FALSE]

  ## -- deliverables: from here we work IN the report dir (delivered to gstore) -
  setwdNew(reportName)
  ## full table
  ezWrite.table(res, "SpliceWiz_ASE_all.txt", row.names = FALSE)
  ## one TSV per event type
  for (et in sort(unique(res$EventType))) {
    ezWrite.table(
      res[res$EventType == et, , drop = FALSE],
      paste0("SpliceWiz_ASE_", et, ".txt"),
      row.names = FALSE
    )
  }
  ## filtered NxtSE for downstream reuse. saveRDS() refuses an NxtSE that still
  ## references out-of-memory (COV) data, so persist with qs2 instead. (Its cov
  ## paths point into the scratch workDir and are gone after the job, so coverage
  ## is not re-plottable from the delivered object -- the ASE/PSI data remain.)
  qs2::qs_save(se, "spliceWiz_filtered_se.qs2")

  ## -- plain bundle for the hosted 'exploreSpliceWiz' Shiny app ---------------
  ## Everything here is plain data.frames/matrix so the explorer needs NO
  ## SpliceWiz in its image. Contract (input file spliceWiz_explore.qs2):
  ##   list(ase, psi[events x samples, 0-1], coldata[sample,condition], meta).
  ## The app is reached at http://fgcz-shiny.uzh.ch/exploreSpliceWiz?data=<gstore
  ## relative report dir> (built below and in SpliceWizApp.rb#next_dataset).
  inc <- SummarizedExperiment::assay(se, "Included")
  exc <- SummarizedExperiment::assay(se, "Excluded")
  psi <- as.data.frame(inc / (inc + exc))
  rownames(psi) <- SummarizedExperiment::rowData(se)$EventName
  bundle <- list(
    ase = res,
    psi = psi,
    coldata = data.frame(
      sample = colnames(se),
      condition = as.character(colData(se)$condition),
      stringsAsFactors = FALSE
    ),
    meta = list(
      comparison = param$comparison, sampleGroup = param$sampleGroup,
      refGroup = param$refGroup, aseMethod = param$aseMethod,
      FDR = param$FDR, deltaPSI = param$deltaPSI
    )
  )
  qs2::qs_save(bundle, "spliceWiz_explore.qs2")

  ## explorer URL from the gstore-relative result dir (matches SpliceWizApp.rb's
  ## next_dataset report_file = File.join(@result_dir, comparison))
  if (ezIsSpecified(param$resultDir)) {
    param$exploreURL <- paste0(
      "http://fgcz-shiny.uzh.ch/exploreSpliceWiz?data=",
      file.path(param$resultDir, param$comparison)
    )
  }

  ## -- FGCZ Quarto report -----------------------------------------------------
  ## quiet = FALSE so a failing chunk surfaces the real quarto/knitr error in the
  ## SUSHI job log (the default quiet = TRUE reduces any failure to an opaque
  ## "Error returned by quarto CLI" -- see the note in makeQuartoReport()).
  makeQuartoReport(
    output = output,
    param = param,
    aseResult = res,
    se = se,
    qmdFile = "SpliceWiz.qmd",
    reportTitle = param$comparison,
    quiet = FALSE
  )
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodSpliceWiz(input=NA, output=NA, param=NA)
##' @description Use this reference class to run SpliceWiz differential
##'   alternative-splicing analysis on STAR BAM files for a two-group comparison.
##' @seealso \code{\link{prepareSpliceWizRef}}
EzAppSpliceWiz <-
  setRefClass(
    "EzAppSpliceWiz",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSpliceWiz
        name <<- "EzAppSpliceWiz"
        appDefaults <<- rbind(
          aseMethod = ezFrame(
            Type = "character",
            DefaultValue = "limma",
            Description = "differential ASE engine: limma (ASE_limma), DESeq (ASE_DESeq) or edgeR (ASE_edgeR)"
          ),
          IRmode = ezFrame(
            Type = "character",
            DefaultValue = "all",
            Description = "intron-retention handling for the ASE test: all, annotated or annotated_binary"
          ),
          strandMode = ezFrame(
            Type = "character",
            DefaultValue = "both",
            Description = "library strandedness; 'both' collates strand-agnostically, otherwise SpliceWiz auto-detects the direction"
          ),
          transcriptTypes = ezFrame(
            Type = "charVector",
            DefaultValue = "protein_coding",
            Description = "restrict the reference to these transcript biotypes before event discovery (as in FeatureCounts)"
          ),
          grouping2 = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "optional secondary co-variate (a Factor/Numeric dataset column) added to the model as ~ grouping2 + condition"
          ),
          FDR = ezFrame(
            Type = "numeric",
            DefaultValue = 0.05,
            Description = "adjusted-p cutoff for calling a significant event"
          ),
          deltaPSI = ezFrame(
            Type = "numeric",
            DefaultValue = 0.1,
            Description = "minimum |delta PSI| for calling a significant event"
          ),
          topN = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "number of top events to draw coverage plots for"
          )
        )
      }
    )
  )
