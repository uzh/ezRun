###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## rMATS-turbo: replicate MultivariateAnalysis of Transcript Splicing. Detection
## and quantification of differential alternative splicing (SE / MXE / A5SS /
## A3SS / RI) from STAR BAMs, for a two-group comparison. rMATS-turbo is a
## Python/C executable (conda env gi_rmats-turbo); this app shells out to it via
## a staged bash script and parses the *.MATS.JC/JCEC.txt tables back into R for
## the FGCZ Quarto report. Modelled on app-SpliceWiz.R.

## Fixed locations of the rMATS-turbo install. We call rmats.py directly through
## `conda run -n` rather than the run_rmats wrapper: the wrapper's
## setup_environment.sh does `source ~/.bashrc`, whose conda-init block re-activates
## a different env AFTER ours, so `python rmats.py` then runs under the wrong
## interpreter and the python3.14-only `rmatspipeline` C extension fails to import
## (crashed the first real SUSHI run). `conda run -n` activates the target env
## deterministically, independent of the currently-active env and of ~/.bashrc.
RMATS_CONDA <- "/usr/local/ngseq/miniforge3/bin/conda"
RMATS_ENV   <- "gi_rmats-turbo"
RMATS_PY    <- "/usr/local/ngseq/srcm/Tools/rmats-turbo/rmats.py"

##' @title Prepare the GTF fed to rMATS-turbo (optional biotype prefilter)
##' @description rMATS needs only a GTF (no index build). If
##'   \code{param$transcriptTypes} restricts the biotypes, a filtered GTF is
##'   written into \code{workDir} and returned; otherwise the shared GTF that the
##'   \code{EzRef} points at is returned unchanged. Mirrors the transcriptTypes
##'   prefilter of \code{\link{prepareSpliceWizRef}} (as in FeatureCounts/DESeq2).
##' @param param the SUSHI/ezRun parameter list; \code{param$ezRef} must resolve
##'   \code{refFeatureFile} (GTF).
##' @param workDir scratch directory (never delivered to gstore) for the filtered GTF.
##' @return path to the GTF rMATS should use.
prepareRmatsGtf <- function(param, workDir = getwd()) {
  gtfFile <- param$ezRef["refFeatureFile"]
  if (!ezIsSpecified(gtfFile)) {
    stop("refFeatureFile not defined")
  }
  buildGtf <- gtfFile
  if (ezIsSpecified(param$transcriptTypes)) {
    ## SUSHI collapses a multi_selection param to a comma-joined string; split it.
    types <- sort(unique(trimws(unlist(strsplit(as.character(param$transcriptTypes), ",")))))
    annoFile <- param$ezRef["refAnnotationFile"]
    if (ezIsSpecified(annoFile) && file.exists(annoFile)) {
      buildGtf <- file.path(workDir, "genes_filtered.gtf")
      if (!file.exists(buildGtf)) {
        rtracklayer::export(gtfByTxTypes(param, types), buildGtf)
      }
    } else {
      ezLog("rMATS: refAnnotationFile missing -- ignoring transcriptTypes, using full GTF")
    }
  }
  buildGtf
}

##' @title Map an FGCZ strandMode to an rMATS --libType
##' @param strandMode one of both / sense / antisense.
##' @return the rMATS library type string.
rmatsLibType <- function(strandMode) {
  switch(as.character(strandMode),
    both = "fr-unstranded",
    sense = "fr-secondstrand",
    antisense = "fr-firststrand",
    "fr-unstranded"
  )
}

ezMethodRmatsTurbo <- function(input = NA, output = NA, param = NA) {
  cwd <- getwd() # the SUSHI $SCRATCH_DIR; only the report dir is shipped to gstore
  on.exit(setwd(cwd), add = TRUE)
  reportName <- basename(output$getColumn("Report"))
  stopifnot(param$sampleGroup != param$refGroup)

  ## -- two-group resolution (app-SpliceWiz.R / app-JunctionSeq.R pattern) ------
  param$grouping <- sub("\\ \\[.*", "", param$grouping) # strip Ruby [Factor] suffix
  condition <- make.names(input$getColumn(param$grouping))
  names(condition) <- input$getNames()
  param$sampleGroup <- make.names(param$sampleGroup)
  param$refGroup <- make.names(param$refGroup)
  samplesUse <- condition %in% c(param$sampleGroup, param$refGroup)
  input <- input$subset(samplesUse)
  condition <- condition[samplesUse]

  bamFiles <- input$getFullPaths("BAM")

  ## -- scratch working dir: ALL heavy intermediates live here, NOT in the report
  ## dir. SUSHI ships the whole report dir to gstore recursively and rm -rf's the
  ## rest of the scratch dir, so anything here is cleaned and never delivered. --
  workDir <- file.path(cwd, "rmats_work")
  dir.create(workDir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workDir, recursive = TRUE, force = TRUE), add = TRUE)

  ## -- stage BAMs (+ .bai) into workDir (gstore automount can time out mid-job).
  ## getBamLocally copies to getwd() and returns a bare basename, so run it from
  ## workDir and rebuild absolute paths (rMATS needs absolute BAM paths).
  setwd(workDir)
  localBams <- vapply(bamFiles, getBamLocally, character(1))
  setwd(cwd)
  copied <- localBams != bamFiles
  localBams[copied] <- file.path(workDir, basename(localBams[copied]))

  ## -- GTF (optional transcriptTypes prefilter) -------------------------------
  gtf <- prepareRmatsGtf(param, workDir)

  ## -- b1 = sampleGroup, b2 = refGroup so IncLevelDifference = sample - ref -----
  s1 <- names(condition)[condition == param$sampleGroup]
  s2 <- names(condition)[condition == param$refGroup]
  b1Bams <- unname(localBams[s1])
  b2Bams <- unname(localBams[s2])
  b1File <- file.path(workDir, "b1.txt")
  b2File <- file.path(workDir, "b2.txt")
  writeLines(paste(b1Bams, collapse = ","), b1File)
  writeLines(paste(b2Bams, collapse = ","), b2File)

  ## -- read type: from param, else from the dataset 'paired' column ------------
  readType <- as.character(param$readType)
  if (!ezIsSpecified(readType) || readType == "auto") {
    readType <- "paired"
    if (input$hasColumn("paired")) {
      pv <- tolower(as.character(input$getColumn("paired")))
      if (length(pv) > 0 && all(pv %in% c("false", "0", "no", "f"))) {
        readType <- "single"
      }
    }
  }

  libType <- rmatsLibType(param$strandMode)
  nCores <- as.integer(param$cores)
  outDir <- file.path(workDir, "rmats_out")
  tmpDir <- file.path(workDir, "rmats_tmp")

  ## -- optional statistical model / novel splice sites / free options ----------
  extraFlags <- character(0)
  if (ezIsSpecified(param$statModel) && param$statModel == "paired") {
    extraFlags <- c(extraFlags, "--paired-stats")
  }
  if (ezIsSpecified(param$statModel) && param$statModel == "darts") {
    extraFlags <- c(extraFlags, "--darts-model")
  }
  if (ezIsSpecified(param$novelSS) && tolower(as.character(param$novelSS)) == "true") {
    extraFlags <- c(extraFlags, "--novelSS")
  }
  if (ezIsSpecified(param$specialOptions)) {
    extraFlags <- c(extraFlags, param$specialOptions)
  }

  ## -- readLength: use the param if set, else auto-detect in the run script ----
  readLen <- if (ezIsSpecified(param$readLength)) as.character(param$readLength) else ""

  ## -- staged bash script ------------------------------------------------------
  ## A script (not an inline ezSystem string) also sidesteps the ezSystem
  ## pipe+single-quote guard.
  ## The read-length probe uses samtools from the loaded Tools/samtools module
  ## (stays on PATH; the rMATS env has no samtools), outside the conda env.
  ## rMATS itself runs via `conda run -n gi_rmats-turbo python rmats.py` -- NOT the
  ## run_rmats wrapper (see RMATS_CONDA note above); `conda run` propagates the
  ## real exit code so ezSystem still detects failure, and --no-capture-output
  ## streams progress to the job log. 'set -e' (NOT pipefail, so 'samtools | head'
  ## closing the pipe early does not abort) guards the whole script.
  script <- file.path(workDir, "run_rmats.sh")
  scriptLines <- c(
    "#!/bin/bash",
    "set -e",
    sprintf('FIRSTBAM=%s', shQuote(b1Bams[1])),
    sprintf('READLEN=%s', shQuote(readLen)),
    'if [ -z "$READLEN" ]; then',
    '  READLEN=$(samtools view "$FIRSTBAM" | head -n 2000 | awk \'{ if (length($10) > m) m = length($10) } END { print (m > 0 ? m : 100) }\')',
    'fi',
    'echo "rMATS-turbo readLength=$READLEN"',
    paste(
      shQuote(RMATS_CONDA), "run --no-capture-output -n", RMATS_ENV,
      "python", shQuote(RMATS_PY),
      sprintf("--b1 %s --b2 %s", shQuote(b1File), shQuote(b2File)),
      sprintf("--gtf %s", shQuote(gtf)),
      sprintf("-t %s", readType),
      '--readLength "$READLEN" --variable-read-length',
      sprintf("--libType %s", libType),
      sprintf("--nthread %d", nCores),
      sprintf("--cstat %s", as.character(param$cstat)),
      sprintf("--od %s --tmp %s", shQuote(outDir), shQuote(tmpDir)),
      paste(extraFlags, collapse = " ")
    )
  )
  writeLines(scriptLines, script)
  ezSystem(paste("bash", shQuote(script)))

  ## -- parse the rMATS output tables ------------------------------------------
  eventTypes <- c("SE", "MXE", "A5SS", "A3SS", "RI")
  sharedCols <- c("ID", "GeneID", "geneSymbol", "chr", "strand",
                  "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference")
  ## combine the shared columns across event types into one report table; the
  ## full per-event-type files (with the event-specific coordinate columns) are
  ## delivered verbatim below.
  combineShared <- function(count) {
    lst <- lapply(eventTypes, function(et) {
      f <- file.path(outDir, paste0(et, ".MATS.", count, ".txt"))
      if (!file.exists(f)) return(NULL)
      df <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
      if (nrow(df) == 0) return(NULL)
      keep <- df[, intersect(sharedCols, colnames(df)), drop = FALSE]
      keep$EventType <- et
      keep
    })
    lst <- lst[!vapply(lst, is.null, logical(1))]
    if (length(lst) == 0) return(NULL)
    do.call(rbind, lst)
  }

  res <- combineShared("JC")
  if (is.null(res)) {
    stop("rMATS-turbo produced no events -- check the STAR BAMs, read length and strandedness")
  }
  ## rMATS wraps gene ids/symbols in double quotes; drop them.
  res$geneSymbol <- gsub('"', "", res$geneSymbol)
  res$GeneID <- gsub('"', "", res$GeneID)
  res$gene_name <- res$geneSymbol
  res$deltaPSI <- suppressWarnings(as.numeric(res$IncLevelDifference))
  res$abs_deltaPSI <- abs(res$deltaPSI)
  res$PValue <- suppressWarnings(as.numeric(res$PValue))
  res$FDR <- suppressWarnings(as.numeric(res$FDR))
  res$EventName <- paste(res$EventType, res$gene_name, res$ID, sep = ":")
  res$EventRegion <- paste0(res$chr, ":", res$strand)
  ## surface the useful columns near the front
  front <- c("gene_name", "EventType", "EventName", "EventRegion",
             "deltaPSI", "abs_deltaPSI", "PValue", "FDR")
  res <- res[, c(intersect(front, colnames(res)),
                 setdiff(colnames(res), front)), drop = FALSE]

  ## -- explorer bundle: per-replicate PSI matrix (events x samples) -----------
  psiMat <- matrix(
    NA_real_, nrow = nrow(res), ncol = length(c(s1, s2)),
    dimnames = list(res$EventName, c(s1, s2))
  )
  for (i in seq_len(nrow(res))) {
    v1 <- suppressWarnings(as.numeric(strsplit(res$IncLevel1[i], ",")[[1]]))
    v2 <- suppressWarnings(as.numeric(strsplit(res$IncLevel2[i], ",")[[1]]))
    if (length(v1) == length(s1)) psiMat[i, s1] <- v1
    if (length(v2) == length(s2)) psiMat[i, s2] <- v2
  }
  psi <- as.data.frame(psiMat, check.names = FALSE)
  coldata <- data.frame(
    sample = c(s1, s2),
    condition = c(rep(param$sampleGroup, length(s1)), rep(param$refGroup, length(s2))),
    stringsAsFactors = FALSE
  )

  ## -- deliverables: from here we work IN the report dir (delivered to gstore) -
  setwdNew(reportName)
  ezWrite.table(res, "rMATS_JC_all.txt", row.names = FALSE)
  ## full per-event-type files (JC + JCEC, all coordinate/count columns), verbatim
  for (et in eventTypes) {
    for (count in c("JC", "JCEC")) {
      f <- file.path(outDir, paste0(et, ".MATS.", count, ".txt"))
      if (file.exists(f)) {
        file.copy(f, paste0("rMATS_", et, ".MATS.", count, ".txt"), overwrite = TRUE)
      }
    }
  }
  summaryFile <- file.path(outDir, "summary.txt")
  if (file.exists(summaryFile)) file.copy(summaryFile, "rMATS_summary.txt", overwrite = TRUE)
  qs2::qs_save(res, "rmatsTurbo_result.qs2")

  ## plain bundle for the hosted 'exploreRmatsTurbo' Shiny app (no rMATS needed).
  ## Contract (input file rmatsTurbo_explore.qs2):
  ##   list(ase, psi[events x samples, 0-1], coldata[sample,condition], meta).
  bundle <- list(
    ase = res,
    psi = psi,
    coldata = coldata,
    meta = list(
      comparison = param$comparison, sampleGroup = param$sampleGroup,
      refGroup = param$refGroup, statModel = param$statModel,
      FDR = param$FDR, deltaPSI = param$deltaPSI
    )
  )
  qs2::qs_save(bundle, "rmatsTurbo_explore.qs2")

  if (ezIsSpecified(param$resultDir)) {
    param$exploreURL <- paste0(
      "http://fgcz-shiny.uzh.ch/exploreRmatsTurbo?data=",
      file.path(param$resultDir, param$comparison)
    )
  }

  ## -- FGCZ Quarto report -----------------------------------------------------
  ## quiet = FALSE so a failing chunk surfaces the real quarto/knitr error in the
  ## SUSHI job log (default quiet = TRUE reduces it to an opaque CLI error).
  makeQuartoReport(
    output = output,
    param = param,
    aseResult = res,
    coldata = coldata,
    qmdFile = "RmatsTurbo.qmd",
    reportTitle = param$comparison,
    quiet = FALSE
  )
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodRmatsTurbo(input=NA, output=NA, param=NA)
##' @description Use this reference class to run rMATS-turbo differential
##'   alternative-splicing analysis on STAR BAM files for a two-group comparison.
##' @seealso \code{\link{prepareRmatsGtf}}
EzAppRmatsTurbo <-
  setRefClass(
    "EzAppRmatsTurbo",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodRmatsTurbo
        name <<- "EzAppRmatsTurbo"
        appDefaults <<- rbind(
          readLength = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "RNA-seq read length; leave blank to auto-detect from the BAMs (always run with --variable-read-length)"
          ),
          readType = ezFrame(
            Type = "character",
            DefaultValue = "auto",
            Description = "paired / single; 'auto' reads the dataset 'paired' column"
          ),
          strandMode = ezFrame(
            Type = "character",
            DefaultValue = "both",
            Description = "library strandedness -> rMATS --libType (both=fr-unstranded, sense=fr-secondstrand, antisense=fr-firststrand)"
          ),
          statModel = ezFrame(
            Type = "character",
            DefaultValue = "default",
            Description = "rMATS statistical model: default (unpaired) / paired (--paired-stats) / darts (--darts-model)"
          ),
          novelSS = ezFrame(
            Type = "character",
            DefaultValue = "false",
            Description = "detect novel (unannotated) splice sites (--novelSS)"
          ),
          cstat = ezFrame(
            Type = "character",
            DefaultValue = "0.0001",
            Description = "rMATS --cstat: cutoff splicing difference for the null hypothesis test (0 <= c < 1)"
          ),
          transcriptTypes = ezFrame(
            Type = "charVector",
            DefaultValue = "protein_coding",
            Description = "restrict the reference GTF to these transcript biotypes before event discovery (as in FeatureCounts)"
          ),
          grouping2 = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "optional secondary co-variate (reserved; rMATS models a single two-group factor)"
          ),
          FDR = ezFrame(
            Type = "numeric",
            DefaultValue = 0.05,
            Description = "adjusted-p cutoff for calling a significant event (report-side)"
          ),
          deltaPSI = ezFrame(
            Type = "numeric",
            DefaultValue = 0.1,
            Description = "minimum |IncLevelDifference| for calling a significant event (report-side)"
          ),
          topN = ezFrame(
            Type = "numeric",
            DefaultValue = 20,
            Description = "number of top events to draw PSI boxplots for"
          )
        )
      }
    )
  )
