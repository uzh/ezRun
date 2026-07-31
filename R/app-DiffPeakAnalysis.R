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

  grouping <- input$getColumn(param$grouping)
  grouping <- grouping[grouping %in% c(param$refGroup, param$sampleGroup)]

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
    grouping2 <- input$getColumn(param$grouping2)[names(grouping)]
  }

  ## robustness: both comparison groups must have samples
  groupCounts <- table(
    factor(grouping, levels = c(param$refGroup, param$sampleGroup))
  )
  emptyGroups <- names(groupCounts)[groupCounts == 0]
  if (length(emptyGroups) > 0) {
    stop(
      "No samples found for group(s): ", paste(emptyGroups, collapse = ", "),
      " in column '", param$grouping, "'. Observed values: ",
      paste(sort(unique(input$getColumn(param$grouping))), collapse = ", ")
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

  commonCols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

  countFiles <- input$getFullPathsList("Count")
  featureCounts <- loadCountFiles(countFiles, grouping, commonCols)

  dds <- generateDESeqDS(featureCounts, commonCols, grouping, grouping2)
  peakAnno <- annotateConsensusPeaks(
    gtfFile = param$ezRef@refFeatureFile,
    fastaFile = param$ezRef@refFastaFile,
    peakFile = countFiles[[1]],
    tool = param$annotationMethod,
    cores = param$cores
  )
  outDir <- file.path(basename(output$getColumn('ResultFolder')))
  cd <- getwd()
  on.exit(setwd(cd), add = TRUE)
  setwdNew(outDir)
  reportTitle <- paste0(
    "Differential Peak Analysis: ",
    param$sampleGroup, " over ", param$refGroup
  )
  makeQuartoReport(
    output = output,
    param = param,
    peakAnno = peakAnno,
    dds = dds,
    qmdFile = "DiffPeakAnalysis.qmd",
    htmlFile = "00index.html",
    reportTitle = reportTitle,
    buttons = TRUE,
    use.qs2 = TRUE
  )
}


#' @template app-template
##' @templateVar method ezMethodDiffPeakAnalysis(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run a differential expression analysis with the application edgeR on two groups.
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
            DefaultValue = "chippeakanno",
            Description = "peak annotation method"
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
##' @return A list with \code{dds} (fitted), \code{res} (annotated tibble),
##' \code{noReplicates} (logical), \code{blockingUsed} (logical),
##' \code{grouping2Name}, \code{sampleGroup} and \code{refGroup}.
##' @export
runDiffPeakDESeq <- function(dds, sampleGroup, refGroup, peakAnno = NULL,
                             grouping2Name = NULL) {
  library(DESeq2)

  ## keep only the two comparison groups and drop unused factor levels
  keep <- as.character(dds$group) %in% c(sampleGroup, refGroup)
  dds <- dds[, keep]
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

  ## residual degrees of freedom on the CHOSEN design decide replicate handling
  colDataDf <- as.data.frame(SummarizedExperiment::colData(dds))
  modelMat <- stats::model.matrix(DESeq2::design(dds), data = colDataDf)
  noReplicates <- (nrow(modelMat) - ncol(modelMat)) < 1

  if (!noReplicates) {
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
    mlDds <- DESeq2::estimateSizeFactors(dds)
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

  rawRes <- DESeq2::results(mlDds, contrast = c("group", sampleGroup, refGroup))
  res <- rawRes |>
    as.data.frame() |>
    tibble::rownames_to_column("peakId") |>
    tibble::as_tibble()
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
    refGroup = refGroup
  )
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


##' @description generate peaks annotation file using the method specified with tool
annotateConsensusPeaks <- function(gtfFile, peakFile, fastaFile, tool, cores) {
  switch(
    as.character(tool),
    chippeakanno = {
      require(ChIPpeakAnno)
      require(GenomicRanges)
      require(rtracklayer)
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
      require(ChIPseeker)
      require(GenomicRanges)
      require(rtracklayer)
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
      cmd = paste(
        "annotatePeaks.pl",
        bedFileName,
        fastaFile,
        "-gtf",
        gtfFile,
        "-cpu",
        cores,
        "> annotatedPeaks.txt"
      )
      system(cmd)
      annotHomer <- ezRead.table('annotatedPeaks.txt', row.names = NULL)
      colnames(annotHomer)[1] <- 'peakId'
      annotHomer <- annotHomer %>% rename("geneName" = "Gene Name")
      return(annotHomer)
    }
  )
}
