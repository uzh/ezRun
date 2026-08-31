###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## SpliceWiz: differential alternative splicing (SE / MXE / A5SS / A3SS / RI /
## AFE / ALE) from STAR BAMs, for a two-group comparison.

##' @title Prepare (build once, cache) a SpliceWiz reference
##' @description Builds a SpliceWiz reference from the genome FASTA + GTF that the
##'   \code{EzRef} points at, and caches it next to the GTF keyed by the reference
##'   build -- so it is shared across projects and built only once. Mirrors the
##'   lock/sentinel protocol of \code{\link{getSTARReference}}. If the reference
##'   directory beside the GTF is not writable, the reference is built into the
##'   current (scratch) working directory for this run instead.
##' @param param the SUSHI/ezRun parameter list; \code{param$ezRef} must resolve
##'   \code{refFeatureFile} (GTF) and \code{refFastaFile} (genome FASTA).
##' @return the path to a complete SpliceWiz reference directory.
prepareSpliceWizRef <- function(param) {
  gtfFile <- param$ezRef["refFeatureFile"]
  fastaFile <- param$ezRef["refFastaFile"]
  if (!ezIsSpecified(gtfFile)) {
    stop("refFeatureFile not defined")
  }

  ## cache dir beside the GTF, keyed by the reference build (as getSTARReference does)
  refDir <- sub("\\.gtf$", "_SpliceWizRef", gtfFile)
  ## the reference is complete iff SpliceWiz.ref.gz sits inside it (per buildRef docs)
  sentinel <- function(d) file.exists(file.path(d, "SpliceWiz.ref.gz"))

  ## if the shared location is not writable, build into scratch for this run
  if (!sentinel(refDir) && file.access(dirname(refDir), mode = 2) != 0) {
    refDir <- file.path(getwd(), "SpliceWizRef")
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
    gtf = gtfFile,
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

ezMethodSpliceWiz <- function(input = NA, output = NA, param = NA) {
  library(SpliceWiz)
  library(SummarizedExperiment)

  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add = TRUE)
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

  ## -- stage BAMs (+ .bai) locally: gstore automount can time out mid-job -----
  localBams <- vapply(bamFiles, getBamLocally, character(1))
  toClean <- localBams[localBams != bamFiles]
  if (length(toClean) > 0) {
    on.exit(
      file.remove(c(toClean, paste0(toClean, ".bai"))),
      add = TRUE
    )
  }

  ## -- reference (build once, cache) -----------------------------------------
  refDir <- prepareSpliceWizRef(param)

  nCores <- as.integer(param$cores)

  ## -- per-sample processBAM --------------------------------------------------
  pbDir <- file.path(getwd(), "sw_pb")
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
  nxtDir <- file.path(getwd(), "sw_nxt")
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
  ## keeps them beside the qs2-saved NxtSE so plotCoverage still resolves them in
  ## the separate report render process. Symlinks do NOT work: the setter
  ## normalizePath()s them back to the original sibling dir, re-triggering it.
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
  ## optional batch covariate for a ~ batch + condition model
  if (ezIsSpecified(param$batch)) {
    batchCol <- sub("\\ \\[.*", "", param$batch)
    if (input$hasColumn(batchCol)) {
      colData(se)$batch <- factor(make.names(input$getColumn(batchCol)[colnames(se)]))
      aseArgs$se <- se
      aseArgs$batch1 <- "batch"
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

  ## -- write outputs ----------------------------------------------------------
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
  ## references out-of-memory (COV) data, so persist with qs2 instead.
  qs2::qs_save(se, "spliceWiz_filtered_se.qs2")

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
          batch = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "optional factor column to include as a batch covariate (~ batch + condition)"
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
