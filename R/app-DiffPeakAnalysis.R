###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodDiffPeakAnalysis <- function(input = NA, output = NA, param = NA) {
  if (identical(as.character(param$sampleGroup), as.character(param$refGroup))) {
    stop(
      "sampleGroup and refGroup must be different (both are '",
      param$sampleGroup, "')."
    )
  }
  if (!param$grouping %in% input$colNames) {
    stop(
      "The grouping column '", param$grouping,
      "' is not present in the input dataset. Available columns: ",
      paste(input$colNames, collapse = ", ")
    )
  }
  if (!param$normMethod %in% c("DESeq2", "readsInPeaks", "TMM")) {
    stop(
      "normMethod must be one of DESeq2, readsInPeaks, TMM; got '",
      param$normMethod, "'."
    )
  }

  groupingAll <- input$getColumn(param$grouping)
  grouping <- groupingAll[groupingAll %in% c(param$refGroup, param$sampleGroup)]

  ## robustness: both comparison groups must have samples
  groupCounts <- table(
    factor(grouping, levels = c(param$refGroup, param$sampleGroup))
  )
  emptyGroups <- names(groupCounts)[groupCounts == 0]
  if (length(emptyGroups) > 0) {
    stop(
      "No samples found for group(s): ", paste(emptyGroups, collapse = ", "),
      " in column '", param$grouping, "'. Observed values: ",
      paste(sort(unique(groupingAll)), collapse = ", ")
    )
  }
  if (sum(groupCounts) < 2) {
    stop(
      "A differential peak comparison needs at least two samples; found ",
      sum(groupCounts), "."
    )
  }
  if (any(groupCounts < 2)) {
    ezLog(paste0(
      "Comparison ", param$sampleGroup, " vs ", param$refGroup,
      " has group(s) without replicates (",
      paste(sprintf("%s=%d", names(groupCounts), as.integer(groupCounts)),
        collapse = ", "),
      "). Continuing with the DESeq2 no-replicate workaround; ",
      "p-values will be conservative."
    ))
  }

  ## samples entering the model: the two comparison groups, or every sample
  ## with a group label when fitAllSamples is set (better dispersion estimates)
  compSamples <- names(grouping)
  if (isTRUE(param$fitAllSamples)) {
    fitGrouping <- groupingAll[!is.na(groupingAll) & groupingAll != ""]
  } else {
    fitGrouping <- grouping
  }

  ## optional second factor for paired/blocked testing (mirrors the DESeq2 app's
  ## 'grouping2'); resolved to a per-sample vector named by sample.
  grouping2 <- NULL
  if (ezIsSpecified(param$grouping2)) {
    if (!param$grouping2 %in% input$colNames) {
      stop(
        "The grouping2 column '", param$grouping2,
        "' is not present in the input dataset. Available columns: ",
        paste(input$colNames, collapse = ", ")
      )
    }
    grouping2 <- input$getColumn(param$grouping2)[names(fitGrouping)]
  }

  commonCols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

  countFiles <- input$getFullPathsList("Count")
  featureCounts <- loadCountFiles(countFiles, fitGrouping, commonCols)

  dds <- generateDESeqDS(featureCounts, commonCols, fitGrouping, grouping2)
  peakAnno <- annotateConsensusPeaks(
    gtfFile = param$ezRef@refFeatureFile,
    fastaFile = param$ezRef@refFastaFile,
    peakFile = countFiles[[1]],
    tool = param$annotationMethod,
    cores = param$cores
  )

  ## BigWig tracks are optional; resolve them before changing directory
  bwFiles <- NULL
  bwCol <- intersect(c("BigWig", "BigWigFile"), input$colNames)
  if (length(bwCol) > 0) {
    bwFiles <- input$getFullPaths(bwCol[1], checkExists = FALSE)[compSamples]
    bwFiles <- bwFiles[!is.na(bwFiles) & file.exists(bwFiles)]
    if (length(bwFiles) == 0) {
      bwFiles <- NULL
    }
  }

  outDir <- file.path(basename(output$getColumn('ResultFolder')))
  cd <- getwd()
  on.exit(setwd(cd), add = TRUE)
  setwdNew(outDir)

  deseqOut <- runDiffPeakDESeq(
    dds, param$sampleGroup, param$refGroup, peakAnno,
    grouping2Name = param$grouping2,
    normMethod = param$normMethod,
    fitAllSamples = isTRUE(param$fitAllSamples),
    lfcShrinkType = param$lfcShrink,
    lfcTest = isTRUE(param$lfcTest),
    lfcThreshold = param$lfcThreshold
  )
  mlDds <- deseqOut$dds
  peakTable <- makeDiffPeakTable(
    deseqOut$res, deseqOut$noReplicates,
    lfcThreshold = param$lfcThreshold, fdrThreshold = param$fdrThreshold
  )
  utils::write.table(
    peakTable |> dplyr::arrange(padj),
    file = "differential_peaks.txt", sep = "\t", quote = FALSE,
    row.names = FALSE
  )

  ## variance-stabilised signal of the two comparison groups (PCA, heatmap)
  vsd <- diffPeakVST(mlDds, deseqOut$noReplicates)
  vsd <- vsd[, colnames(vsd) %in% compSamples]

  deseqInfo <- deseqOut[setdiff(names(deseqOut), c("dds", "res"))]
  normDiag <- diffPeakNormDiagnostic(
    mlDds, param$sampleGroup, param$refGroup, peakTable,
    lfcThreshold = param$lfcThreshold, fdrThreshold = param$fdrThreshold
  )
  peakGc <- tryCatch(
    peakGcContent(peakTable, param$ezRef@refFastaFile),
    error = function(e) {
      ezLog(paste("GC content could not be computed:", conditionMessage(e)))
      NULL
    }
  )

  motifRes <- NULL
  if (isTRUE(param$runMotifs)) {
    motifRes <- tryCatch(
      runHomerKnownMotifs(
        peakTable,
        fastaFile = param$ezRef@refFastaFile,
        cores = param$cores,
        maxPeaks = param$enrichMaxPeaks,
        deNovo = isTRUE(param$motifDeNovo),
        motifSet = resolveHomerMotifSet(param$refBuild),
        noReplicates = deseqOut$noReplicates
      ),
      error = function(e) {
        ezLog(paste("HOMER motif analysis failed:", conditionMessage(e)))
        NULL
      }
    )
  }

  goRes <- NULL
  if (isTRUE(param$runGoOra)) {
    goRes <- tryCatch(
      diffPeakGoOra(
        peakTable, param$refBuild,
        maxPeaks = param$enrichMaxPeaks,
        noReplicates = deseqOut$noReplicates
      ),
      error = function(e) {
        ezLog(paste("GO enrichment failed:", conditionMessage(e)))
        NULL
      }
    )
  }

  profiles <- NULL
  if (!is.null(bwFiles)) {
    profiles <- tryCatch(
      bigwigPeakProfiles(
        bwFiles, grouping[names(bwFiles)], peakTable,
        noReplicates = deseqOut$noReplicates
      ),
      error = function(e) {
        ezLog(paste("BigWig profiles failed:", conditionMessage(e)))
        NULL
      }
    )
  }

  reportTitle <- paste0(
    "Differential Peak Analysis: ",
    param$sampleGroup, " over ", param$refGroup
  )
  makeQuartoReport(
    output = output,
    param = param,
    peakAnno = peakAnno,
    dds = slimDiffPeakDds(mlDds, commonCols),
    vsd = vsd,
    peakTable = peakTable,
    deseqInfo = deseqInfo,
    normDiag = normDiag,
    peakGc = peakGc,
    motifRes = motifRes,
    goRes = goRes,
    profiles = profiles,
    qmdFile = "DiffPeakAnalysis.qmd",
    htmlFile = "00index.html",
    reportTitle = reportTitle,
    buttons = TRUE,
    use.qs2 = TRUE
  )
}


#' @template app-template
##' @templateVar method ezMethodDiffPeakAnalysis(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run a differential peak analysis
##' (ATAC-seq / ChIP-seq) on consensus-peak counts of two groups with DESeq2,
##' including known-motif enrichment (HOMER), GO over-representation, BigWig
##' signal profiles and normalisation / GC diagnostics.
EzAppDiffPeakAnalysis <-
  setRefClass(
    "EzAppDiffPeakAnalysis",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodDiffPeakAnalysis
        name <<- "EzAppDiffPeakAnalysis"
        appDefaults <<- rbind(
          annotationMethod = ezFrame(
            Type = "character",
            DefaultValue = "homer",
            Description = "peak annotation method: homer, chippeakanno or chipseeker"
          ),
          normMethod = ezFrame(
            Type = "character",
            DefaultValue = "DESeq2",
            Description = "size factors: DESeq2 (median of ratios over peaks), readsInPeaks (library size) or TMM"
          ),
          lfcThreshold = ezFrame(
            Type = "numeric",
            DefaultValue = 1,
            Description = "absolute log2 fold-change threshold for candidates"
          ),
          fdrThreshold = ezFrame(
            Type = "numeric",
            DefaultValue = 0.05,
            Description = "adjusted p-value threshold for candidates"
          ),
          lfcShrink = ezFrame(
            Type = "character",
            DefaultValue = "apeglm",
            Description = "log2 fold-change shrinkage reported in an extra column: apeglm, ashr or none"
          ),
          lfcTest = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "test against |log2FC| > lfcThreshold instead of 0 (DESeq2 lfcThreshold)"
          ),
          fitAllSamples = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "fit DESeq2 on all samples of the dataset and extract the contrast"
          ),
          runMotifs = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "known-motif enrichment of up/down candidates with HOMER"
          ),
          motifDeNovo = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "also run HOMER de novo motif discovery (slow)"
          ),
          runGoOra = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "GO BP over-representation of candidate genes (mouse and human)"
          ),
          enrichMaxPeaks = ezFrame(
            Type = "integer",
            DefaultValue = 2000,
            Description = "maximum number of top candidates per direction used for motif and GO enrichment"
          )
        )
      }
    )
  )

##' @description generate file with all counts for sampleGroup and refGroup
##' Count columns are named by the (unique) sample name so per-sample covariates
##' (e.g. a blocking factor) can be attached in \code{generateDESeqDS}.
loadCountFiles <- function(countFiles, grouping, commonCols) {
  countFilesSubset <- countFiles[names(grouping)]

  loadAllTables <- imap(countFilesSubset, function(file_i, listName) {
    data.table::fread(file_i, data.table = FALSE) %>%
      rename(!!listName := matchCounts)
  })

  purrr::reduce(loadAllTables, full_join, by = commonCols)
}

##' @description generate DESeqDataSet from the counts table
##' @param grouping Named (by sample) vector of the primary group per sample.
##' @param grouping2 Optional named (by sample) vector of a second/blocking
##' factor; when supplied the design becomes \code{~ group + grouping2}.
generateDESeqDS <- function(featureCounts, commonCols, grouping, grouping2 = NULL) {
  library(DESeq2)
  countCols <- setdiff(colnames(featureCounts), commonCols)

  countData <- featureCounts %>%
    select(all_of(countCols)) %>%
    as.data.frame()
  countData <- round(countData)

  ## count columns are sample names; build colData from the dataset maps.
  ## The blocking factor is carried as a colData column but the object is always
  ## built with ~group: runDiffPeakDESeq decides whether ~ group + grouping2 is
  ## usable (full rank, enough residual df) before switching the design, so a
  ## confounded second factor never trips DESeq2's construction-time rank check.
  colData <- data.frame(
    sample = countCols,
    group = factor(as.character(grouping[countCols])),
    row.names = countCols,
    stringsAsFactors = FALSE
  )
  if (!is.null(grouping2)) {
    colData$grouping2 <- factor(as.character(grouping2[countCols]))
  }

  dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = ~group
  )
  rowData(dds) <- featureCounts[, commonCols]
  rownames(dds) <- featureCounts$Geneid
  dds$Condition <- dds$group
  dds
}

##' @title Size factors for differential peak testing
##' @description Set DESeq2 size factors with one of three methods.
##' \code{DESeq2}: median of ratios over peaks (assumes most peaks do not
##' change). \code{readsInPeaks}: total reads in peaks per sample (library
##' size; robust to a global shift of accessibility in many peaks).
##' \code{TMM}: edgeR trimmed mean of M-values on the peak counts. The
##' library-size based factors are scaled to a geometric mean of 1; DESeq2
##' factors are returned unchanged.
##' @param dds A \code{DESeqDataSet}.
##' @param normMethod One of "DESeq2", "readsInPeaks", "TMM".
##' @return The \code{DESeqDataSet} with size factors set.
##' @export
setDiffPeakSizeFactors <- function(dds, normMethod = "DESeq2") {
  DESeq2::sizeFactors(dds) <- diffPeakSizeFactors(DESeq2::counts(dds), normMethod)
  dds
}

##' @describeIn setDiffPeakSizeFactors size factors from a count matrix
diffPeakSizeFactors <- function(counts, normMethod = "DESeq2") {
  libSize <- colSums(counts)
  if (normMethod == "DESeq2") {
    return(DESeq2::estimateSizeFactorsForMatrix(counts))
  }
  sizeFactors <- switch(
    normMethod,
    readsInPeaks = libSize,
    TMM = libSize * edgeR::calcNormFactors(counts, lib.size = libSize, method = "TMM"),
    stop("unsupported normMethod: ", normMethod)
  )
  sizeFactors / exp(mean(log(sizeFactors)))
}

##' @title Run DESeq2 for a two-group differential peak comparison
##' @description
##' Fit DESeq2 for one \code{sampleGroup} vs \code{refGroup} comparison and
##' return the annotated result table together with the fitted object. This is
##' robust to comparisons \emph{without biological replicates}: when the design
##' has no residual degrees of freedom, dispersions are estimated by treating
##' the two conditions as a single group (design \code{~1}) and a Wald test is
##' run against the real design. This absorbs the between-condition signal into
##' the dispersion, so p-values are conservative (usually non-significant),
##' while log2 fold changes stay informative for exploratory ranking. See the
##' DESeq2 vignette section on analysis without replicates.
##' @param dds A \code{DESeqDataSet} with a \code{group} column.
##' @param sampleGroup Group of interest (numerator of the contrast).
##' @param refGroup Reference group (denominator of the contrast).
##' @param peakAnno Optional peak annotation data frame joined by \code{peakId}.
##' @param grouping2Name Optional display name of the blocking factor (used only
##' for messaging in the report).
##' @param normMethod Size-factor method, see \code{setDiffPeakSizeFactors}.
##' @param fitAllSamples Keep all samples (all groups) in the fit and extract
##' the contrast; otherwise only the two comparison groups are used.
##' @param lfcShrinkType "apeglm", "ashr" or "none"; the shrunken estimate is
##' added as \code{log2FoldChange_shrunk}, the MLE stays in \code{log2FoldChange}.
##' @param lfcTest Test against |log2FC| > \code{lfcThreshold} instead of 0.
##' @param lfcThreshold Threshold used when \code{lfcTest} is TRUE.
##' @return A list with \code{dds} (fitted), \code{res} (annotated tibble),
##' \code{noReplicates} (logical), \code{blockingUsed} (logical),
##' \code{grouping2Name}, \code{sampleGroup}, \code{refGroup},
##' \code{normMethod}, \code{fitAllSamples}, \code{lfcShrinkType},
##' \code{lfcTest} and \code{nFitSamples}.
##' @export
runDiffPeakDESeq <- function(dds, sampleGroup, refGroup, peakAnno = NULL,
                             grouping2Name = NULL, normMethod = "DESeq2",
                             fitAllSamples = FALSE, lfcShrinkType = "none",
                             lfcTest = FALSE, lfcThreshold = 1) {
  library(DESeq2)

  ## keep only the two comparison groups (unless all samples are fitted) and
  ## drop unused factor levels
  if (!isTRUE(fitAllSamples)) {
    keep <- as.character(dds$group) %in% c(sampleGroup, refGroup)
    dds <- dds[, keep]
  }
  dds$group <- droplevels(factor(as.character(dds$group)))
  dds$group <- relevel(dds$group, ref = refGroup)
  dds$Condition <- dds$group

  ## Decide whether the optional blocking factor (grouping2) can be used for a
  ## paired/blocked design. Drop it (with a warning) when it is constant,
  ## confounded with the groups (rank-deficient), or leaves no residual df.
  blockingUsed <- FALSE
  hasBlock <- "grouping2" %in% colnames(SummarizedExperiment::colData(dds))
  if (hasBlock) {
    dds$grouping2 <- droplevels(factor(dds$grouping2))
    cd <- as.data.frame(SummarizedExperiment::colData(dds))
    usable <- tryCatch({
      if (nlevels(dds$grouping2) < 2) {
        FALSE
      } else {
        mm <- stats::model.matrix(~ group + grouping2, data = cd)
        qr(mm)$rank == ncol(mm) && (nrow(mm) - qr(mm)$rank) >= 1
      }
    }, error = function(e) FALSE)
    if (usable) {
      DESeq2::design(dds) <- ~ group + grouping2
      blockingUsed <- TRUE
    } else {
      DESeq2::design(dds) <- ~group
      ezLog(paste0(
        "Second factor '",
        if (is.null(grouping2Name)) "grouping2" else grouping2Name,
        "' is not usable for ", sampleGroup, " vs ", refGroup,
        " (single level, confounded with the groups, or too few residual ",
        "degrees of freedom); it was dropped and ~group is used instead."
      ))
    }
  } else {
    DESeq2::design(dds) <- ~group
  }

  dds <- setDiffPeakSizeFactors(dds, normMethod)

  ## residual degrees of freedom on the CHOSEN design decide replicate handling
  colDataDf <- as.data.frame(SummarizedExperiment::colData(dds))
  modelMat <- stats::model.matrix(DESeq2::design(dds), data = colDataDf)
  noReplicates <- (nrow(modelMat) - ncol(modelMat)) < 1

  if (!noReplicates) {
    ## DESeq() keeps the pre-set size factors
    mlDds <- DESeq2::DESeq(dds)
  } else {
    ## a block term cannot be estimated without residual df: fall back to ~group
    DESeq2::design(dds) <- ~group
    blockingUsed <- FALSE
    ezLog(paste0(
      "No biological replicates for ", sampleGroup, " vs ", refGroup,
      "; estimating dispersion with a blind (~1) design. p-values are ",
      "conservative and fold changes are for exploratory ranking only."
    ))
    mlDds <- dds
    blind <- mlDds
    DESeq2::design(blind) <- ~1
    blind <- tryCatch(
      DESeq2::estimateDispersions(blind, fitType = "parametric"),
      error = function(e) DESeq2::estimateDispersions(blind, fitType = "local")
    )
    DESeq2::dispersions(mlDds) <- DESeq2::dispersions(blind)
    DESeq2::dispersionFunction(mlDds) <- DESeq2::dispersionFunction(blind)
    mlDds <- DESeq2::nbinomWaldTest(mlDds)
  }

  contrast <- c("group", sampleGroup, refGroup)
  useLfcTest <- isTRUE(lfcTest) && !noReplicates
  if (useLfcTest) {
    rawRes <- DESeq2::results(
      mlDds, contrast = contrast,
      lfcThreshold = lfcThreshold, altHypothesis = "greaterAbs"
    )
  } else {
    rawRes <- DESeq2::results(mlDds, contrast = contrast)
  }
  res <- rawRes |>
    as.data.frame() |>
    tibble::rownames_to_column("peakId") |>
    tibble::as_tibble()

  shrunk <- NULL
  if (!noReplicates && lfcShrinkType %in% c("apeglm", "ashr")) {
    shrunk <- shrinkDiffPeakLfc(mlDds, sampleGroup, refGroup, lfcShrinkType)
  }
  res$log2FoldChange_shrunk <- if (is.null(shrunk)) {
    NA_real_
  } else {
    unname(shrunk[res$peakId])
  }

  if (!is.null(peakAnno)) {
    res <- dplyr::left_join(res, peakAnno, by = "peakId")
  }

  list(
    dds = mlDds,
    res = res,
    noReplicates = noReplicates,
    blockingUsed = blockingUsed,
    grouping2Name = if (blockingUsed) grouping2Name else NULL,
    sampleGroup = sampleGroup,
    refGroup = refGroup,
    normMethod = normMethod,
    fitAllSamples = isTRUE(fitAllSamples),
    lfcShrinkType = if (is.null(shrunk)) "none" else lfcShrinkType,
    lfcTest = useLfcTest,
    nFitSamples = ncol(mlDds)
  )
}

##' @title Reduce a fitted DESeqDataSet to what downstream users need
##' @description Keeps counts, colData (incl. size factors) and the peak
##' coordinates; drops the fitted assays (mu, H, Cook's distances) and the
##' per-peak fit results, which are already in \code{differential_peaks.txt}.
##' The design formula is re-homed to \code{baseenv()}: it was created inside
##' \code{runDiffPeakDESeq}, so its environment holds every intermediate object
##' of that call and would be serialised with it (160 MB instead of 5 MB).
slimDiffPeakDds <- function(mlDds, commonCols) {
  SummarizedExperiment::assays(mlDds) <- SummarizedExperiment::assays(mlDds)["counts"]
  keepCols <- intersect(commonCols, colnames(S4Vectors::mcols(mlDds)))
  S4Vectors::mcols(mlDds) <- S4Vectors::mcols(mlDds)[, keepCols, drop = FALSE]
  designFormula <- DESeq2::design(mlDds)
  environment(designFormula) <- baseenv()
  DESeq2::design(mlDds) <- designFormula
  mlDds
}

##' @title Shrunken log2 fold changes for the comparison
##' @description apeglm needs the contrast as a model coefficient; the
##' reference level is \code{refGroup}, so the coefficient is
##' \code{group_<sampleGroup>_vs_<refGroup>}. If it is not available (or apeglm
##' fails) ashr on the contrast is used. Returns NULL when shrinkage fails.
##' @return Named numeric vector (by peak) or NULL.
shrinkDiffPeakLfc <- function(mlDds, sampleGroup, refGroup, lfcShrinkType) {
  coefName <- paste0("group_", make.names(sampleGroup), "_vs_", make.names(refGroup))
  contrast <- c("group", sampleGroup, refGroup)
  shrunk <- NULL
  if (lfcShrinkType == "apeglm" && coefName %in% DESeq2::resultsNames(mlDds)) {
    shrunk <- tryCatch(
      DESeq2::lfcShrink(mlDds, coef = coefName, type = "apeglm", quiet = TRUE),
      error = function(e) {
        ezLog(paste("apeglm shrinkage failed, trying ashr:", conditionMessage(e)))
        NULL
      }
    )
  }
  if (is.null(shrunk)) {
    shrunk <- tryCatch(
      DESeq2::lfcShrink(mlDds, contrast = contrast, type = "ashr", quiet = TRUE),
      error = function(e) {
        ezLog(paste("ashr shrinkage failed:", conditionMessage(e)))
        NULL
      }
    )
  }
  if (is.null(shrunk)) {
    return(NULL)
  }
  stats::setNames(shrunk$log2FoldChange, rownames(shrunk))
}

##' @title Build the candidate differential-peak table
##' @description
##' Add \code{direction} and \code{candidate} flags to a DESeq2 result table.
##' When replicates are available a peak is a candidate when it passes both the
##' fold-change and FDR thresholds; without replicates (unreliable p-values)
##' the fold-change threshold alone is used.
##' @param res Annotated result tibble from \code{runDiffPeakDESeq}.
##' @param noReplicates Logical, whether the comparison lacked replicates.
##' @param lfcThreshold Absolute log2 fold-change threshold (default 1).
##' @param fdrThreshold Adjusted p-value threshold (default 0.05).
##' @return A tibble with \code{direction} and \code{candidate} columns added.
##' @export
makeDiffPeakTable <- function(res, noReplicates, lfcThreshold = 1, fdrThreshold = 0.05) {
  tbl <- res |>
    dplyr::filter(!is.na(log2FoldChange)) |>
    dplyr::mutate(direction = dplyr::if_else(log2FoldChange > 0, "up", "down"))
  if (isTRUE(noReplicates)) {
    tbl <- tbl |>
      dplyr::mutate(candidate = abs(log2FoldChange) >= lfcThreshold)
  } else {
    tbl <- tbl |>
      dplyr::mutate(
        candidate = !is.na(padj) &
          abs(log2FoldChange) >= lfcThreshold &
          padj < fdrThreshold
      )
  }
  tbl
}

##' @title Variance-stabilizing transformation robust to missing replicates
##' @description Wrapper around \code{varianceStabilizingTransformation} that
##' uses a blind (design \code{~1}) transformation when the comparison has no
##' replicates, which is the only variant that can be computed in that case.
##' @param mlDds A fitted \code{DESeqDataSet}.
##' @param noReplicates Logical, whether the comparison lacked replicates.
##' @return A \code{DESeqTransform} object.
##' @export
diffPeakVST <- function(mlDds, noReplicates) {
  DESeq2::varianceStabilizingTransformation(mlDds, blind = isTRUE(noReplicates))
}

##' @title Rank candidate peaks for downstream enrichment
##' @description Top \code{maxPeaks} candidates of one direction, by adjusted
##' p-value (or absolute fold change without replicates).
##' @return The subset of \code{peakTable}.
topDiffPeaks <- function(peakTable, direction, maxPeaks, noReplicates = FALSE) {
  cand <- peakTable[peakTable$candidate %in% TRUE & peakTable$direction == direction, ,
    drop = FALSE]
  if (nrow(cand) == 0) {
    return(cand)
  }
  if (isTRUE(noReplicates)) {
    ord <- order(-abs(cand$log2FoldChange))
  } else {
    ord <- order(cand$padj, -abs(cand$log2FoldChange), na.last = TRUE)
  }
  cand[head(ord, maxPeaks), , drop = FALSE]
}

##' @title Normalisation sensitivity of a differential peak comparison
##' @description Compares the size factors actually used with the three
##' alternatives and approximates the results under each alternative without
##' refitting: switching size factors shifts every log2 fold change by the same
##' offset (difference of the mean log2 size factors of the two groups), so the
##' Wald statistic is recomputed as (LFC + offset) / lfcSE with the fitted
##' standard errors and adjusted with BH over the tested peaks. This is an
##' approximation intended to reveal global shifts (e.g. many peaks gaining
##' accessibility), not a replacement for a refit.
##' @return A list with \code{sizeFactors} (per sample and method),
##' \code{balance} (up/down candidates per method) and \code{used}.
##' @export
diffPeakNormDiagnostic <- function(mlDds, sampleGroup, refGroup, peakTable,
                                   lfcThreshold = 1, fdrThreshold = 0.05) {
  counts <- DESeq2::counts(mlDds)
  groups <- as.character(mlDds$group)
  methods <- c("DESeq2", "readsInPeaks", "TMM")
  sfUsed <- DESeq2::sizeFactors(mlDds)
  sfAll <- sapply(methods, function(m) {
    tryCatch(diffPeakSizeFactors(counts, m), error = function(e) {
      rep(NA_real_, ncol(counts))
    })
  })
  rownames(sfAll) <- colnames(counts)
  groupOffset <- function(sf) {
    mean(log2(sf[groups == sampleGroup])) - mean(log2(sf[groups == refGroup]))
  }
  offUsed <- groupOffset(sfUsed)
  tested <- peakTable[!is.na(peakTable$padj) & !is.na(peakTable$lfcSE), , drop = FALSE]
  balance <- lapply(methods, function(m) {
    off <- offUsed - groupOffset(sfAll[, m])
    lfc <- tested$log2FoldChange + off
    p <- 2 * stats::pnorm(-abs(lfc / tested$lfcSE))
    padj <- stats::p.adjust(p, method = "BH")
    cand <- padj < fdrThreshold & abs(lfc) >= lfcThreshold
    data.frame(
      method = m,
      lfcOffset = off,
      medianLfc = stats::median(lfc),
      up = sum(cand & lfc > 0, na.rm = TRUE),
      down = sum(cand & lfc < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  sizeFactors <- data.frame(
    sample = colnames(counts),
    group = groups,
    totalReadsInPeaks = colSums(counts),
    used = sfUsed,
    sfAll,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  list(sizeFactors = sizeFactors, balance = do.call(rbind, balance))
}

##' @title GC content and width of peaks
##' @param peakTable Table with \code{peakId}, \code{Chr}, \code{Start},
##' \code{End} (1-based, inclusive).
##' @param fastaFile Indexed genome FASTA.
##' @return data.frame with \code{peakId}, \code{gc} (fraction G+C among
##' A/C/G/T) and \code{width}.
##' @export
peakGcContent <- function(peakTable, fastaFile) {
  fa <- Rsamtools::FaFile(fastaFile)
  seqLengths <- GenomeInfoDb::seqlengths(Rsamtools::seqinfo(fa))
  ok <- peakTable$Chr %in% names(seqLengths) &
    peakTable$End <= seqLengths[as.character(peakTable$Chr)] &
    peakTable$Start >= 1
  ok[is.na(ok)] <- FALSE
  gc <- rep(NA_real_, nrow(peakTable))
  if (any(ok)) {
    gr <- GenomicRanges::GRanges(
      peakTable$Chr[ok],
      IRanges::IRanges(peakTable$Start[ok], peakTable$End[ok])
    )
    seqs <- Biostrings::getSeq(fa, gr)
    nGc <- Biostrings::letterFrequency(seqs, "GC")[, 1]
    nAcgt <- Biostrings::letterFrequency(seqs, "ACGT")[, 1]
    gc[ok] <- ifelse(nAcgt > 0, nGc / nAcgt, NA_real_)
  }
  data.frame(
    peakId = peakTable$peakId,
    gc = gc,
    width = peakTable$End - peakTable$Start + 1,
    stringsAsFactors = FALSE
  )
}

##' @title HOMER known-motif set for a reference build
##' @description "vertebrates" for vertebrate genomes, "all" otherwise.
resolveHomerMotifSet <- function(refBuild) {
  species <- sub("/.*", "", refBuild)
  vertebrates <- c(
    "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus", "Danio_rerio",
    "Gallus_gallus", "Sus_scrofa", "Bos_taurus", "Macaca_mulatta",
    "Macaca_fascicularis", "Canis_familiaris", "Canis_lupus_familiaris",
    "Ovis_aries", "Equus_caballus", "Oryctolagus_cuniculus", "Xenopus_tropicalis",
    "Xenopus_laevis", "Felis_catus", "Callithrix_jacchus", "Oryzias_latipes"
  )
  if (species %in% vertebrates) "vertebrates" else "all"
}

##' @title Write peaks as a BED file for HOMER
##' @description BED is 0-based half-open; peak tables are 1-based inclusive.
writeHomerBed <- function(peaks, file) {
  bed <- data.frame(
    peaks$Chr, as.integer(peaks$Start) - 1L, as.integer(peaks$End),
    peaks$peakId, 0, "+"
  )
  utils::write.table(bed, file, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE)
  invisible(file)
}

##' @title Known-motif enrichment of differential peaks with HOMER
##' @description For up and down candidates separately (top \code{maxPeaks}),
##' runs \code{findMotifsGenome.pl} on the genome FASTA (so any reference
##' works) with 200 bp windows centred on the peaks. The background is the
##' set of tested, non-candidate consensus peaks (sampled to at most
##' \code{maxBackground}), which matches the accessibility and GC profile of
##' the targets better than random genomic sequence. Up and down run in
##' parallel, each with half of \code{cores}.
##' @return A list with \code{up}/\code{down} data.frames (NULL when there are
##' too few candidates), \code{nBackground}, \code{motifSet} and the HOMER
##' result directories.
##' @export
runHomerKnownMotifs <- function(peakTable, fastaFile, cores = 4, maxPeaks = 2000,
                                deNovo = FALSE, motifSet = "vertebrates",
                                noReplicates = FALSE, minPeaks = 20,
                                maxBackground = 50000, size = 200) {
  bgPeaks <- peakTable[!(peakTable$candidate %in% TRUE) & !is.na(peakTable$log2FoldChange), ,
    drop = FALSE]
  if (nrow(bgPeaks) > maxBackground) {
    set.seed(42)
    bgPeaks <- bgPeaks[sort(sample(nrow(bgPeaks), maxBackground)), , drop = FALSE]
  }
  bgFile <- writeHomerBed(bgPeaks, "homer_background.bed")
  prepDir <- file.path(tempdir(), "homer_preparsed")
  dir.create(prepDir, showWarnings = FALSE, recursive = TRUE)
  coresEach <- max(1L, floor(as.integer(cores) / 2L))

  runOne <- function(direction) {
    targets <- topDiffPeaks(peakTable, direction, maxPeaks, noReplicates)
    if (nrow(targets) < minPeaks) {
      ezLog(paste0("Only ", nrow(targets), " ", direction,
        " candidates; motif analysis skipped."))
      return(NULL)
    }
    targetFile <- writeHomerBed(targets, paste0("homer_", direction, ".bed"))
    outDir <- paste0("homerMotifs_", direction)
    cmd <- paste(
      "findMotifsGenome.pl", targetFile, fastaFile, outDir,
      "-size", size, "-bg", bgFile, "-p", coresEach,
      "-mset", motifSet, "-preparsedDir", prepDir,
      if (deNovo) "" else "-nomotif",
      "> ", paste0(outDir, ".log"), "2>&1"
    )
    ezSystem(cmd)
    knownFile <- file.path(outDir, "knownResults.txt")
    if (!file.exists(knownFile)) {
      stop("HOMER did not write ", knownFile)
    }
    known <- utils::read.delim(knownFile, check.names = FALSE, quote = "",
      stringsAsFactors = FALSE)
    parseHomerKnown(known, nTargets = nrow(targets))
  }
  results <- parallel::mclapply(c("up", "down"), runOne, mc.cores = 2)
  names(results) <- c("up", "down")
  failed <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop(paste(unlist(results[failed]), collapse = "; "))
  }
  list(
    up = results$up,
    down = results$down,
    nBackground = nrow(bgPeaks),
    motifSet = motifSet,
    deNovo = deNovo,
    dirs = c(up = "homerMotifs_up", down = "homerMotifs_down")
  )
}

##' @title Tidy a HOMER knownResults.txt table
##' @return data.frame with motif, family, source, consensus, log10 p-value,
##' q-value, target / background percentages and their ratio.
parseHomerKnown <- function(known, nTargets) {
  pct <- function(x) as.numeric(sub("%", "", x, fixed = TRUE))
  nameCol <- "Motif Name"
  targetCol <- grep("^% of Target", colnames(known), value = TRUE)[1]
  bgCol <- grep("^% of Background", colnames(known), value = TRUE)[1]
  motifName <- known[[nameCol]]
  data.frame(
    rank = seq_len(nrow(known)),
    motif = sub("/.*$", "", motifName),
    tf = sub("\\(.*$", "", motifName),
    family = ifelse(grepl("\\(", motifName),
      sub("^[^(]*\\(([^)]*)\\).*$", "\\1", motifName), NA_character_),
    source = sub("^[^/]*/", "", motifName),
    consensus = known[["Consensus"]],
    log10P = -known[["Log P-value"]] / log(10),
    qValue = as.numeric(known[["q-value (Benjamini)"]]),
    pctTarget = pct(known[[targetCol]]),
    pctBackground = pct(known[[bgCol]]),
    enrichment = pct(known[[targetCol]]) / pmax(pct(known[[bgCol]]), 0.01),
    nTargets = nTargets,
    stringsAsFactors = FALSE
  )
}

##' @title Gene identifiers of peaks for over-representation analysis
##' @description Uses Ensembl gene IDs (version stripped) when an annotation
##' column holds them, else gene symbols from \code{geneName}.
##' @return list(ids = character per peak, keyType = "ENSEMBL" or "SYMBOL")
diffPeakGeneIds <- function(peakTable) {
  ensPattern <- "^ENS[A-Z]*G[0-9]+"
  for (col in intersect(c("Entrez ID", "Nearest Ensembl", "geneName", "geneId"),
    colnames(peakTable))) {
    x <- as.character(peakTable[[col]])
    if (mean(grepl(ensPattern, x[!is.na(x) & x != ""])) > 0.5) {
      ids <- sub("\\.[0-9]+$", "", x)
      ids[!grepl(ensPattern, ids)] <- NA_character_
      return(list(ids = ids, keyType = "ENSEMBL"))
    }
  }
  if ("geneName" %in% colnames(peakTable)) {
    ids <- as.character(peakTable$geneName)
    ids[ids == ""] <- NA_character_
    return(list(ids = ids, keyType = "SYMBOL"))
  }
  NULL
}

##' @title GO biological-process over-representation of candidate genes
##' @description Nearest genes of the top candidates per direction are tested
##' against the nearest genes of all tested peaks as universe (the correct
##' background for peak-based gene lists). Only mouse and human references are
##' supported (org.Mm.eg.db / org.Hs.eg.db); other species return NULL.
##' @return list(up, down) of \code{enrichResult} objects (or NULL), plus
##' \code{keyType}, \code{orgDb} and \code{nUniverse}.
##' @export
diffPeakGoOra <- function(peakTable, refBuild, maxPeaks = 2000, noReplicates = FALSE,
                          ont = "BP") {
  species <- sub("/.*", "", refBuild)
  orgDb <- switch(species,
    Mus_musculus = "org.Mm.eg.db",
    Homo_sapiens = "org.Hs.eg.db",
    NULL
  )
  if (is.null(orgDb) || !requireNamespace(orgDb, quietly = TRUE)) {
    ezLog(paste("GO enrichment is only available for mouse and human; skipped for", species))
    return(NULL)
  }
  geneIds <- diffPeakGeneIds(peakTable)
  if (is.null(geneIds)) {
    ezLog("No gene annotation for GO enrichment; skipped.")
    return(NULL)
  }
  library(clusterProfiler)
  peakTable$goId <- geneIds$ids
  tested <- !is.na(peakTable$padj) | isTRUE(noReplicates)
  universe <- unique(stats::na.omit(peakTable$goId[tested]))
  runOne <- function(direction) {
    top <- topDiffPeaks(peakTable, direction, maxPeaks, noReplicates)
    genes <- unique(stats::na.omit(top$goId))
    if (length(genes) < 10) {
      return(NULL)
    }
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene = genes, universe = universe,
        OrgDb = getExportedValue(orgDb, orgDb), keyType = geneIds$keyType,
        ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
        minGSSize = 10, maxGSSize = 500,
        readable = geneIds$keyType != "SYMBOL"
      ),
      error = function(e) {
        ezLog(paste("enrichGO failed for", direction, ":", conditionMessage(e)))
        NULL
      }
    )
    ego
  }
  list(
    up = runOne("up"),
    down = runOne("down"),
    keyType = geneIds$keyType,
    orgDb = orgDb,
    ont = ont,
    nUniverse = length(universe)
  )
}

##' @title BigWig signal profiles around differential peaks
##' @description For the top \code{maxPeaks} up and down candidates, extracts
##' BigWig signal in \code{binSize} bins +-\code{extend} bp around the peak
##' centre (EnrichedHeatmap::normalizeToMatrix) per sample, then averages per
##' group. Rows keep the candidate ranking.
##' @param bwFiles Named (by sample) BigWig paths.
##' @param groups Named (by sample) group labels.
##' @return list(up, down), each a list with per-group mean matrices
##' (\code{normalizedMatrix}) and a long per-sample metaprofile data.frame.
##' @export
bigwigPeakProfiles <- function(bwFiles, groups, peakTable, maxPeaks = 1000,
                               extend = 2000, binSize = 50, noReplicates = FALSE) {
  bwInfo <- GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(bwFiles[[1]]))
  seqLevelsBw <- Reduce(intersect, lapply(bwFiles, function(f) {
    GenomeInfoDb::seqlevels(rtracklayer::BigWigFile(f))
  }))
  runOne <- function(direction) {
    top <- topDiffPeaks(peakTable, direction, maxPeaks, noReplicates)
    top <- top[top$Chr %in% seqLevelsBw, , drop = FALSE]
    if (nrow(top) < 10) {
      return(NULL)
    }
    mid <- round((top$Start + top$End) / 2)
    centres <- GenomicRanges::GRanges(top$Chr, IRanges::IRanges(mid, mid))
    names(centres) <- top$peakId
    windows <- GenomicRanges::resize(centres, width = 2 * extend + 1, fix = "center")
    GenomeInfoDb::seqlevels(windows) <- GenomeInfoDb::seqlevels(bwInfo)[
      GenomeInfoDb::seqlevels(bwInfo) %in% GenomeInfoDb::seqlevels(windows)]
    GenomeInfoDb::seqinfo(windows) <- bwInfo[GenomeInfoDb::seqlevels(windows)]
    windows <- GenomicRanges::trim(windows)
    mats <- lapply(bwFiles, function(f) {
      sig <- rtracklayer::import(f, which = windows)
      EnrichedHeatmap::normalizeToMatrix(
        sig, centres, value_column = "score", extend = extend, w = binSize,
        mean_mode = "w0", background = 0
      )
    })
    groupMats <- lapply(split(names(bwFiles), groups[names(bwFiles)]), function(s) {
      m <- Reduce(`+`, lapply(mats[s], unclass)) / length(s)
      attributes(m) <- attributes(mats[[s[1]]])
      m
    })
    pos <- seq(-extend + binSize / 2, extend - binSize / 2, by = binSize)
    meta <- do.call(rbind, lapply(names(mats), function(s) {
      data.frame(
        sample = s, group = groups[[s]], position = pos,
        signal = colMeans(unclass(mats[[s]]), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    list(groupMats = groupMats, meta = meta, peakIds = top$peakId)
  }
  list(up = runOne("up"), down = runOne("down"), extend = extend, binSize = binSize)
}


##' @description generate peaks annotation file using the method specified with tool
annotateConsensusPeaks <- function(gtfFile, peakFile, fastaFile, tool, cores) {
  switch(
    as.character(tool),
    chippeakanno = {
      library(ChIPpeakAnno)
      library(GenomicRanges)
      library(rtracklayer)
      gtf <- rtracklayer::import(gtfFile)
      if ('gene' %in% unique(gtf$type)) {
        idx = gtf$type == 'gene'
      } else if ('transcript' %in% unique(gtf$type)) {
        idx = gtf$type == 'transcript'
      } else if ('start_codon' %in% unique(gtf$type)) {
        idx = gtf$type == 'start_codon'
      } else {
        ezLog('gtf is incompatabible. Peak annotation skipped!')
        return(NULL)
      }
      gtf = gtf[idx]
      if (grepl('gtf$', gtfFile)) {
        names_gtf = make.unique(gtf$'gene_id')
      } else {
        names_gtf = make.unique(gtf$'ID')
      }
      names(gtf) = names_gtf
      myPeaks = ezRead.table(peakFile)
      peaksRD = makeGRangesFromDataFrame(
        myPeaks,
        keep.extra.columns = TRUE,
        start.field = "Start",
        end.field = "End",
        seqnames.field = "Chr"
      )
      annotChIPpeak <- annotatePeakInBatch(
        peaksRD,
        AnnotationData = gtf,
        output = 'nearestStart',
        multiple = FALSE,
        FeatureLocForDistance = 'TSS'
      )
      annotChIPpeak <- as.data.frame(annotChIPpeak)
      annotChIPpeak <- annotChIPpeak %>%
        rename(
          "feature_start" = "start_position",
          "feature_end" = "end_position"
        )
      annotChIPpeak <- annotChIPpeak %>%
        rename("peakId" = "peak", "geneName" = "feature")
      return(annotChIPpeak)
    },
    chipseeker = {
      library(ChIPseeker)
      library(GenomicRanges)
      library(rtracklayer)
      myTxDB <- txdbmaker::makeTxDbFromGFF(file = gtfFile, format = 'gtf')
      myPeaks = ezRead.table(peakFile)
      myPeaks$peakId = rownames(myPeaks)
      peaksRD = makeGRangesFromDataFrame(
        myPeaks,
        keep.extra.columns = TRUE,
        start.field = "Start",
        end.field = "End",
        seqnames.field = "Chr"
      )
      annotChIPseeker <- annotatePeak(
        peaksRD,
        TxDb = myTxDB,
        tssRegion = c(-1000, 1000),
        verbose = FALSE
      )
      annotChIPseeker <- data.frame(annotChIPseeker@anno)
      keepColChIPSeeker <- c(
        "seqnames",
        "peakId",
        "annotation",
        "geneId",
        "transcriptId",
        "distanceToTSS",
        "geneChr",
        "geneStart",
        "geneEnd",
        "geneLength",
        "geneStrand"
      )
      annotChIPseeker <- annotChIPseeker[, keepColChIPSeeker]
      annotChIPseeker <- annotChIPseeker %>% rename("geneName" = "geneId")
      return(annotChIPseeker)
    },
    homer = {
      myPeaks <- ezRead.table(peakFile)
      bedFileCols <- c("Chr", "Start", "End")
      bedFile <- myPeaks[, bedFileCols]
      bedFile$Names <- rownames(myPeaks)
      bedFile$Score <- 0
      bedFile$Strand <- myPeaks[["Strand"]]
      bedFileName <- "peaks.bed"
      write.table(
        bedFile,
        bedFileName,
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      annoFile <- "annotatedPeaks.txt"
      cmd = paste(
        "annotatePeaks.pl",
        bedFileName,
        fastaFile,
        "-gtf",
        gtfFile,
        "-cpu",
        cores,
        ">",
        annoFile
      )
      ezSystem(cmd)
      if (!file.exists(annoFile) || file.size(annoFile) == 0) {
        stop("HOMER annotatePeaks.pl did not produce ", annoFile)
      }
      annotHomer <- ezRead.table(annoFile, row.names = NULL)
      colnames(annotHomer)[1] <- 'peakId'
      annotHomer <- annotHomer %>% rename("geneName" = "Gene Name")
      return(annotHomer)
    }
  )
}
