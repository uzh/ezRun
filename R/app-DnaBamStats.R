###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Run a dataset-level DNA BAM stats report
##' @description
##' Maintained dataset-level DNA BAM stats workflow.
##' This app collects per-sample DNA BAM QC metrics into one shared
##' `resultList`, renders a single cohort-level report, and keeps the report
##' stable even when only a subset of metrics is available.
##' @template input-template
##' @template output-template
##' @template param-template
##' @template htmlFile-template
##' @return Character scalar indicating success.
##' @template roxygen-template
ezMethodDnaBamStats <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  library(ATACseqQC)
  library(GenomicAlignments)
  library(S4Vectors)

  param <- normalize_dna_bamstats_param(param)

  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)

  setwdNew(basename(output$getColumn("Report")))

  dataset <- input$meta
  samples <- input$getNames()
  bamFiles <- input$getFullPaths("BAM")

  resultList <- list()
  for (sm in samples) {
    ezLog(sm)
    nReads <- if ("Read Count" %in% colnames(dataset)) {
      dataset[sm, "Read Count"]
    } else {
      NA_real_
    }
    resultList[[sm]] <- tryCatch(
      {
        get_dna_bamstats_skeleton(
          bamFile = bamFiles[sm],
          sampleName = sm,
          param = param,
          nReads = nReads
        )
      },
      error = function(e) {
        list(
          error = paste(
            "Failed while collecting DNA BAM stats for sample",
            sm,
            ":",
            conditionMessage(e)
          )
        )
      }
    )
    if ("error" %in% names(resultList[[sm]])) {
      writeErrorReport(
        htmlFile,
        param = param,
        error = resultList[[sm]][["error"]]
      )
      return("Error")
    }
  }

  reportData <- build_dna_bamstats_report_data(
    dataset = dataset,
    resultList = resultList
  )
  write_dna_bamstats_support_files(
    reportData = reportData
  )
  if (isTRUE(param$runQualimap)) {
    reportData$qualimapMulti <- get_dna_qualimap_multi_sample_summary(
      samples = samples,
      resultList = resultList
    )
  }

  makeRmdReport(
    dataset = dataset,
    param = param,
    resultList = resultList,
    reportData = reportData,
    rmdFile = "DnaBamStats.Rmd",
    reportTitle = "DNA BAM Stats",
    selfContained = TRUE
  )

  rm(resultList)
  gc()
  return("Success")
}

##' @title Collect minimal per-sample DNA BAM stats for the report skeleton
##' @description
##' Gather the DNA BAM QC values needed to render the maintained cohort-level
##' report while degrading cleanly when optional external-tool steps fail.
##' @param bamFile Full path to the BAM file.
##' @param sampleName Sample name used in the report.
##' @param param ezRun parameter list.
##' @param nReads Total read count from the input dataset metadata.
##' @return Named list with per-sample report values.
##' @template roxygen-template
get_dna_bamstats_skeleton <- function(
  bamFile,
  sampleName,
  param,
  nReads = NA
) {
  result <- list(
    sampleName = sampleName,
    bamFile = bamFile,
    nReads = nReads,
    multiMatchInFileTable = NULL,
    mappingRate = NA_real_,
    dupRate = NA_real_,
    avgCoverage = NA_real_,
    errorRate = NA_real_,
    insertRate = NA_real_,
    delRate = NA_real_,
    allDuplicates = NA_real_,
    optDuplicates = NA_real_,
    readPairs = NA_real_,
    opticalDupRate = NA_real_,
    pctCov1x = NA_real_,
    pctCov10x = NA_real_,
    pctCov20x = NA_real_,
    meanInsertSize = NA_real_,
    medianInsertSize = NA_real_,
    meanMapQuality = NA_real_,
    gcContent = NA_real_,
    fragmentSizePlot = NULL,
    libComplexityPlot = NULL,
    notes = "Dataset-level DNA BAM QC metrics collected for cohort reporting.",
    qualimapDir = NULL,
    picardMetricsFile = NULL
  )

  if (file.exists(bamFile) && !is.na(nReads)) {
    result$multiMatchInFileTable <- tryCatch(
      {
        getBamMultiMatching(
          param = param,
          bamFile = bamFile,
          nReads = nReads
        )
      },
      error = function(e) {
        ezLog(
          "Failed to collect multi-matching counts for ",
          sampleName,
          ": ",
          conditionMessage(e),
          level = "warn"
        )
        NULL
      }
    )
  }

  if (isTRUE(param$paired) && file.exists(bamFile)) {
    pairedPlotStats <- get_dna_paired_end_plots(
      bamFile = bamFile,
      sampleName = sampleName
    )
    result[names(pairedPlotStats)] <- pairedPlotStats
    if ("notes" %in% names(pairedPlotStats)) {
      result$notes <- c(result$notes, pairedPlotStats$notes)
    }
  }

  if (isTRUE(param$runQualimap) && file.exists(bamFile)) {
    qualimapStats <- tryCatch(
      {
        get_dna_qualimap_stats(
          bamFile = bamFile,
          sampleName = sampleName,
          param = param,
          nReads = nReads
        )
      },
      error = function(e) {
        ezLog(
          "Failed to collect Qualimap metrics for ",
          sampleName,
          ": ",
          conditionMessage(e),
          level = "warn"
        )
        list(
          notes = paste(
            "Qualimap metrics unavailable for",
            sampleName,
            "-",
            conditionMessage(e)
          )
        )
      }
    )
    if ("notes" %in% names(qualimapStats)) {
      result$notes <- c(result$notes, qualimapStats$notes)
      qualimapStats$notes <- NULL
    }
    result[names(qualimapStats)] <- qualimapStats
  }

  if (isTRUE(param$runPicard) && file.exists(bamFile)) {
    picardStats <- tryCatch(
      {
        get_dna_picard_dup_stats(
          bamFile = bamFile,
          sampleName = sampleName,
          param = param
        )
      },
      error = function(e) {
        ezLog(
          "Failed to collect Picard duplicate metrics for ",
          sampleName,
          ": ",
          conditionMessage(e),
          level = "warn"
        )
        list(
          notes = paste(
            "Picard duplicate metrics unavailable for",
            sampleName,
            "-",
            conditionMessage(e)
          )
        )
      }
    )
    if ("notes" %in% names(picardStats)) {
      result$notes <- c(result$notes, picardStats$notes)
      picardStats$notes <- NULL
    }
    result[names(picardStats)] <- picardStats
  }

  result$notes <- unique(result$notes)

  return(result)
}

##' @title Build report-side helper objects for DNA BAM stats
##' @description Collect the cohort tables, QC thresholds, flags, and stable
##' sample ordering used by the report layout.
##' @param dataset Input dataset metadata.
##' @param resultList Named list of per-sample results.
##' @return Named list with report helper objects.
##' @template roxygen-template
build_dna_bamstats_report_data <- function(dataset, resultList) {
  samples <- rownames(dataset)
  overviewTable <- get_dna_bamstats_overview_table(
    samples = samples,
    resultList = resultList
  )
  qcThresholds <- get_dna_bamstats_qc_thresholds()
  fullSummaryTable <- get_dna_bamstats_full_summary_table(
    overviewTable = overviewTable,
    qcThresholds = qcThresholds
  )
  reviewOrder <- get_dna_bamstats_review_order(fullSummaryTable)
  fullSummaryTable <- fullSummaryTable[
    match(reviewOrder, fullSummaryTable$Sample),
    ,
    drop = FALSE
  ]

  reportData <- list(
    overviewTable = overviewTable,
    qcThresholds = qcThresholds,
    fullSummaryTable = fullSummaryTable,
    quickViewTable = get_dna_bamstats_quick_view_table(
      fullSummaryTable = fullSummaryTable
    ),
    inputDatasetTable = get_dna_bamstats_input_dataset_table(dataset),
    resourceTable = get_dna_bamstats_resource_table(
      samples = samples,
      resultList = resultList
    ),
    flagSummaryTable = get_dna_bamstats_flag_summary_table(
      fullSummaryTable = fullSummaryTable,
      qcThresholds = qcThresholds
    ),
    reviewOrder = reviewOrder,
    qualimapMulti = NULL
  )

  return(reportData)
}

##' @title Build cohort-level DNA BAM stats overview table
##' @description Assemble the main per-sample numeric metrics into one cohort
##' table that can be reused by the report and written as a sidecar file.
##' @param samples Character vector of sample names.
##' @param resultList Named list of per-sample results.
##' @return Data frame with one row per sample and one column per metric.
##' @template roxygen-template
get_dna_bamstats_overview_table <- function(samples, resultList) {
  ## Display column -> per-sample result field. Extend this map to add metrics.
  metricFields <- c(
    MappingRate = "mappingRate",
    DuplicationRate = "dupRate",
    OpticalDupRate = "opticalDupRate",
    AverageCoverage = "avgCoverage",
    PctCov1x = "pctCov1x",
    PctCov10x = "pctCov10x",
    PctCov20x = "pctCov20x",
    MeanInsertSize = "meanInsertSize",
    MedianInsertSize = "medianInsertSize",
    MeanMapQuality = "meanMapQuality",
    GcContent = "gcContent",
    ErrorRate = "errorRate",
    InsertionRate = "insertRate",
    DeletionRate = "delRate",
    AllDuplicates = "allDuplicates",
    OpticalDuplicates = "optDuplicates",
    ReadPairs = "readPairs"
  )

  overviewTable <- data.frame(Sample = samples, stringsAsFactors = FALSE)
  for (colName in names(metricFields)) {
    field <- metricFields[[colName]]
    overviewTable[[colName]] <- vapply(
      samples,
      function(sm) {
        value <- resultList[[sm]][[field]]
        if (is.null(value) || length(value) != 1) {
          return(NA_real_)
        }
        as.numeric(value)
      },
      numeric(1)
    )
  }

  return(overviewTable)
}

##' @title Return default QC thresholds for DNA BAM stats
##' @description Provide the conservative absolute review thresholds agreed for
##' the first version of the report.
##' @return Data frame describing thresholded QC flags.
##' @template roxygen-template
get_dna_bamstats_qc_thresholds <- function() {
  qcThresholds <- data.frame(
    Flag = c(
      "LowMappingRate",
      "LowAverageCoverage",
      "HighDuplicationRate",
      "HighOpticalDuplicationRate"
    ),
    Metric = c(
      "MappingRate",
      "AverageCoverage",
      "DuplicationRate",
      "OpticalDupRate"
    ),
    Direction = c("below", "below", "above", "above"),
    Threshold = c(70, 10, 50, 20),
    Unit = c("%", "x", "%", "%"),
    stringsAsFactors = FALSE
  )

  return(qcThresholds)
}

##' @title Build full DNA BAM stats summary table
##' @description Combine numeric cohort metrics with per-flag boolean columns
##' and compact sample-level flag summaries.
##' @param overviewTable Cohort metric table.
##' @param qcThresholds Threshold definition table.
##' @return Data frame with one raw summary row per sample.
##' @template roxygen-template
get_dna_bamstats_full_summary_table <- function(overviewTable, qcThresholds) {
  fullSummaryTable <- overviewTable

  for (i in seq_len(nrow(qcThresholds))) {
    flagName <- qcThresholds$Flag[i]
    metricName <- qcThresholds$Metric[i]
    thresholdValue <- qcThresholds$Threshold[i]
    direction <- qcThresholds$Direction[i]
    metricValues <- fullSummaryTable[[metricName]]

    fullSummaryTable[[flagName]] <- if (direction == "below") {
      !is.na(metricValues) & metricValues < thresholdValue
    } else {
      !is.na(metricValues) & metricValues > thresholdValue
    }
  }

  fullSummaryTable$FlagCount <- vapply(
    seq_len(nrow(fullSummaryTable)),
    function(i) {
      sum(unlist(fullSummaryTable[i, qcThresholds$Flag, drop = FALSE]))
    },
    numeric(1)
  )
  fullSummaryTable$TriggeredFlags <- vapply(
    seq_len(nrow(fullSummaryTable)),
    function(i) {
      get_dna_bamstats_triggered_flags_text(
        sampleRow = fullSummaryTable[i, , drop = FALSE],
        qcThresholds = qcThresholds,
        formatted = FALSE
      )
    },
    character(1)
  )

  fullSummaryTable <- fullSummaryTable[, c(
    "Sample",
    "FlagCount",
    "TriggeredFlags",
    "MappingRate",
    "DuplicationRate",
    "OpticalDupRate",
    "AverageCoverage",
    "PctCov1x",
    "PctCov10x",
    "PctCov20x",
    "MeanInsertSize",
    "MedianInsertSize",
    "MeanMapQuality",
    "GcContent",
    "ErrorRate",
    "InsertionRate",
    "DeletionRate",
    qcThresholds$Flag,
    "AllDuplicates",
    "OpticalDuplicates",
    "ReadPairs"
  )]

  return(fullSummaryTable)
}

##' @title Build quick-view DNA BAM stats table
##' @description Format the compact all-sample triage table shown in the first
##' report tab.
##' @param fullSummaryTable Raw full summary table.
##' @return Data frame with formatted quick-view values.
##' @template roxygen-template
get_dna_bamstats_quick_view_table <- function(fullSummaryTable) {
  quickViewTable <- data.frame(
    Sample = fullSummaryTable$Sample,
    FlagCount = fullSummaryTable$FlagCount,
    MappingRate = fullSummaryTable$MappingRate,
    DuplicationRate = fullSummaryTable$DuplicationRate,
    OpticalDupRate = fullSummaryTable$OpticalDupRate,
    AverageCoverage = fullSummaryTable$AverageCoverage,
    stringsAsFactors = FALSE
  )

  return(quickViewTable)
}

##' @title Build a reduced input dataset table for the report
##' @description Keep only the provenance columns that matter for this report
##' and shorten file paths for readability.
##' @param dataset Input dataset metadata.
##' @return Data frame with reduced plain-text provenance columns.
##' @template roxygen-template
get_dna_bamstats_input_dataset_table <- function(dataset) {
  inputDatasetTable <- data.frame(
    Sample = rownames(dataset),
    stringsAsFactors = FALSE
  )

  candidateCols <- c("Species", "refBuild", "Read Count", "BAM", "BAI")
  keepCols <- intersect(candidateCols, colnames(dataset))
  for (colName in keepCols) {
    values <- dataset[, colName]
    if (colName %in% c("BAM", "BAI")) {
      values <- basename(values)
    }
    cleanColName <- gsub(" ", "", colName, fixed = TRUE)
    inputDatasetTable[[cleanColName]] <- values
  }

  return(inputDatasetTable)
}

##' @title Build the per-sample resource table for the report
##' @description Collect the drilldown artifacts linked from the Resources tab.
##' @param samples Character vector of sample names.
##' @param resultList Named list of per-sample results.
##' @return Data frame with one row per sample and one drilldown column per
##' resource type.
##' @template roxygen-template
get_dna_bamstats_resource_table <- function(samples, resultList) {
  resourceTable <- data.frame(
    Sample = samples,
    QualimapReport = vapply(
      samples,
      function(sm) {
        qualimapDir <- resultList[[sm]]$qualimapDir
        if (is.null(qualimapDir) || is.na(qualimapDir) || qualimapDir == "") {
          return(NA_character_)
        }
        file.path(qualimapDir, "qualimapReport.html")
      },
      character(1)
    ),
    PicardMetrics = vapply(
      samples,
      function(sm) {
        if (is.null(resultList[[sm]]$picardMetricsFile)) {
          return(NA_character_)
        }
        resultList[[sm]]$picardMetricsFile
      },
      character(1)
    ),
    FragmentSizePlot = vapply(
      samples,
      function(sm) {
        if (is.null(resultList[[sm]]$fragmentSizePlot)) {
          return(NA_character_)
        }
        resultList[[sm]]$fragmentSizePlot
      },
      character(1)
    ),
    LibraryComplexityPlot = vapply(
      samples,
      function(sm) {
        if (is.null(resultList[[sm]]$libComplexityPlot)) {
          return(NA_character_)
        }
        resultList[[sm]]$libComplexityPlot
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  return(resourceTable)
}

##' @title Build the per-flag summary table for quick-view plotting
##' @description Count how many distinct samples trigger each conservative
##' absolute QC rule.
##' @param fullSummaryTable Raw full summary table.
##' @param qcThresholds Threshold definition table.
##' @return Data frame with one row per flag and a sample count.
##' @template roxygen-template
get_dna_bamstats_flag_summary_table <- function(fullSummaryTable, qcThresholds) {
  flagSummaryTable <- data.frame(
    Flag = qcThresholds$Flag,
    SamplesAffected = vapply(
      qcThresholds$Flag,
      function(flagName) {
        sum(fullSummaryTable[[flagName]], na.rm = TRUE)
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  return(flagSummaryTable)
}

##' @title Compute stable review-priority sample ordering
##' @description Sort samples by flag burden first and then by a small set of
##' severity proxies that stay stable across all metric plots.
##' @param fullSummaryTable Raw full summary table.
##' @return Character vector of sample names in display order.
##' @template roxygen-template
get_dna_bamstats_review_order <- function(fullSummaryTable) {
  ordering <- order(
    -fullSummaryTable$FlagCount,
    fullSummaryTable$MappingRate,
    -fullSummaryTable$DuplicationRate,
    fullSummaryTable$AverageCoverage,
    fullSummaryTable$Sample,
    na.last = TRUE
  )

  return(fullSummaryTable$Sample[ordering])
}

##' @title Build a triggered-flag text summary for one sample
##' @description Assemble one compact flag string using the canonical flag
##' ordering and optionally including formatted threshold comparisons.
##' @param sampleRow One-row data frame from the full summary table.
##' @param qcThresholds Threshold definition table.
##' @param formatted Whether to format values with report units.
##' @return Character scalar containing the triggered flag summary.
##' @template roxygen-template
get_dna_bamstats_triggered_flags_text <- function(
  sampleRow,
  qcThresholds,
  formatted = FALSE
) {
  triggeredFlags <- character()

  for (i in seq_len(nrow(qcThresholds))) {
    flagName <- qcThresholds$Flag[i]
    metricName <- qcThresholds$Metric[i]
    thresholdValue <- qcThresholds$Threshold[i]
    unit <- qcThresholds$Unit[i]
    direction <- qcThresholds$Direction[i]

    if (isTRUE(sampleRow[[flagName]][1])) {
      observedValue <- sampleRow[[metricName]][1]
      if (formatted) {
        observedLabel <- format_dna_bamstats_metric(observedValue, unit = unit)
        thresholdLabel <- format_dna_bamstats_metric(thresholdValue, unit = unit)
      } else {
        observedLabel <- as.character(observedValue)
        thresholdLabel <- as.character(thresholdValue)
      }
      comparator <- if (direction == "below") "<" else ">"
      triggeredFlags <- c(
        triggeredFlags,
        paste0(
          flagName,
          " (",
          observedLabel,
          " ",
          comparator,
          " ",
          thresholdLabel,
          ")"
        )
      )
    }
  }

  if (length(triggeredFlags) == 0) {
    return("")
  }

  return(paste(triggeredFlags, collapse = "; "))
}

##' @title Format a DNA BAM stats metric for report display
##' @description Apply one-decimal formatting and optional unit suffixes.
##' @param values Numeric vector.
##' @param unit Optional unit suffix such as `"%"` or `"x"`.
##' @return Character vector with formatted metric values.
##' @template roxygen-template
format_dna_bamstats_metric <- function(values, unit = "") {
  labels <- ifelse(
    is.na(values),
    NA_character_,
    sprintf("%.1f", values)
  )

  if (unit != "") {
    labels <- ifelse(
      is.na(labels),
      NA_character_,
      paste0(labels, unit)
    )
  }

  return(labels)
}

##' @title Write DNA BAM stats support files
##' @description Write cohort-level helper files that are useful outside the
##' HTML report, including the summary metric table used during comparison.
##' @param reportData Named list of report-side helper objects.
##' @return Invisible `NULL`.
##' @template roxygen-template
write_dna_bamstats_support_files <- function(reportData) {
  if (!is.null(reportData$fullSummaryTable)) {
    ezWrite.table(
      reportData$fullSummaryTable,
      file = "metrics_summary.tsv",
      sep = "\t",
      row.names = FALSE
    )
  }

  return(invisible(NULL))
}

##' @title Normalize DnaBamStats parameter types
##' @description Convert Sushi-style string parameters to the types expected by
##' the ezRun implementation.
##' @param param ezRun parameter list.
##' @return Modified parameter list with normalized logical and numeric values.
##' @template roxygen-template
normalize_dna_bamstats_param <- function(param) {
  logicalKeys <- c("paired", "runQualimap", "runPicard")
  for (key in logicalKeys) {
    if (!is.null(param[[key]])) {
      value <- param[[key]]
      if (is.character(value)) {
        param[[key]] <- tolower(value) %in% c("true", "t", "1", "yes")
      } else {
        param[[key]] <- isTRUE(value)
      }
    }
  }

  numericKeys <- c("cores", "ram", "scratch", "pixelDist")
  for (key in numericKeys) {
    if (!is.null(param[[key]]) && is.character(param[[key]])) {
      suppressWarnings({
        numValue <- as.numeric(param[[key]])
      })
      if (!is.na(numValue)) {
        param[[key]] <- numValue
      }
    }
  }

  return(param)
}

##' @title Create paired-end DNA BAM QC plots
##' @description Generate fragment-size and library-complexity plots for one BAM
##' file when paired-end sequencing is used.
##' @param bamFile Full path to the BAM file.
##' @param sampleName Sample name used for output files.
##' @return Named list with plot paths and note text.
##' @template roxygen-template
get_dna_paired_end_plots <- function(bamFile, sampleName) {
  result <- list(
    fragmentSizePlot = NULL,
    libComplexityPlot = NULL,
    notes = character()
  )

  fragmentSizePlot <- paste0("fragmentSize_", sampleName, ".png")
  tryCatch(
    {
      png(fragmentSizePlot, width = 600, height = 500, res = 90)
      fragSizeDist(bamFile, sampleName)
      dev.off()
      result$fragmentSizePlot <- fragmentSizePlot
      result$notes <- c(
        result$notes,
        paste("Fragment-size plot created:", basename(fragmentSizePlot))
      )
    },
    error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      ezLog(
        "Failed to create fragment-size plot for ",
        sampleName,
        ": ",
        conditionMessage(e),
        level = "warn"
      )
    }
  )

  libComplexityPlot <- paste0("libComplexity_", sampleName, ".png")
  tryCatch(
    {
      png(libComplexityPlot, width = 600, height = 500, res = 90)
      ez_dna_lib_complexity(
        readsDupFreq(bamFile),
        main = paste(sampleName, "Library Complexity")
      )
      dev.off()
      result$libComplexityPlot <- libComplexityPlot
      result$notes <- c(
        result$notes,
        paste("Library-complexity plot created:", basename(libComplexityPlot))
      )
    },
    error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      ezLog(
        "Failed to create library-complexity plot for ",
        sampleName,
        ": ",
        conditionMessage(e),
        level = "warn"
      )
    }
  )

  return(result)
}

##' @title Collect per-sample DNA BAM metrics from Qualimap
##' @description
##' Run `qualimap bamqc` for one BAM file and parse a stable subset of summary
##' values from `genome_results.txt` for the first refactor iteration.
##' @param bamFile Full path to the BAM file.
##' @param sampleName Sample name used for output staging.
##' @param param ezRun parameter list.
##' @param nReads Total read count from the input dataset metadata.
##' @return Named list with parsed Qualimap metrics.
##' @template roxygen-template
get_dna_qualimap_stats <- function(
  bamFile,
  sampleName,
  param,
  nReads = NA
) {
  qualimapDir <- paste0(sampleName, "_qualimap")
  ## Qualimap's launcher injects the Java-8-only `-XX:MaxPermSize` flag unless
  ## JAVA_OPTS is set. That flag is fatal on the jdk21 that the Tools/Picard
  ## module puts on PATH, so pass JAVA_OPTS (which also carries the heap size)
  ## to bypass it and let Qualimap run under jdk21.
  cmd <- paste0(
    "unset DISPLAY; ",
    "JAVA_OPTS=\"-Xms32m -Xmx", param[["ram"]], "G\" ",
    "qualimap bamqc",
    " -bam ", bamFile,
    " -c -nt ", param[["cores"]],
    " -outdir ", qualimapDir
  )
  ezSystem(cmd)

  qmFile <- file.path(qualimapDir, "genome_results.txt")
  stopifnot(file.exists(qmFile))

  allData <- readLines(qmFile)
  parsed <- list(
    dupRate = get_qualimap_percent(
      allData,
      "^.*duplication rate = "
    ),
    errorRate = get_qualimap_numeric(
      allData,
      "^.*general error rate = "
    ),
    insertRate = get_qualimap_percent(
      allData,
      "mapped reads with insertion percentage = "
    ),
    delRate = get_qualimap_percent(
      allData,
      "mapped reads with deletion percentage = "
    ),
    avgCoverage = get_qualimap_coverage(
      allData,
      "mean coverageData = "
    ),
    pctCov1x = get_dna_qualimap_coverage_breadth(allData, 1),
    pctCov10x = get_dna_qualimap_coverage_breadth(allData, 10),
    pctCov20x = get_dna_qualimap_coverage_breadth(allData, 20),
    meanInsertSize = get_qualimap_value(
      allData,
      "mean insert size = ",
      strip = ","
    ),
    medianInsertSize = get_qualimap_value(
      allData,
      "median insert size = ",
      strip = ","
    ),
    meanMapQuality = get_qualimap_value(
      allData,
      "mean mapping quality = "
    ),
    gcContent = get_qualimap_value(
      allData,
      "GC percentage = ",
      strip = "%"
    ),
    qualimapDir = qualimapDir
  )

  qualimapReads <- get_qualimap_read_count(
    allData,
    "^.*number of reads = "
  )
  if (!is.na(qualimapReads) && isTRUE(param$paired)) {
    qualimapReads <- qualimapReads / 2
  }

  if (!is.na(qualimapReads) && !is.na(nReads) && nReads > 0) {
    parsed$mappingRate <- 100 * qualimapReads / nReads
  } else {
    parsed$mappingRate <- NA_real_
  }

  if (is.na(parsed$dupRate)) {
    dupReads <- get_qualimap_read_count(
      allData,
      "^.*number of duplicated reads \\(flagged\\) = "
    )
    if (!is.na(dupReads) && isTRUE(param$paired)) {
      dupReads <- dupReads / 2
    }
    if (!is.na(dupReads) && !is.na(qualimapReads) && qualimapReads > 0) {
      parsed$dupRate <- 100 * dupReads / qualimapReads
    }
  }

  parsed$notes <- paste(
    "Qualimap metrics collected from",
    basename(qmFile)
  )

  return(parsed)
}

##' @title Collect dataset-level Qualimap multi-sample report path
##' @description Run `qualimap multi-bamqc` across all per-sample Qualimap
##' output directories and return the generated report path when available.
##' @param samples Character vector of sample names.
##' @param resultList Named list of per-sample DNA BAM stats results.
##' @return Named list with the multi-sample report path.
##' @template roxygen-template
get_dna_qualimap_multi_sample_summary <- function(samples, resultList) {
  sampleKey <- "sampleKey.txt"
  qualimapInfo <- data.frame(
    SampleName = samples,
    Folder = vapply(
      samples,
      function(sm) {
        resultList[[sm]]$qualimapDir
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  qualimapInfo <- qualimapInfo[
    !is.na(qualimapInfo$Folder) &
      qualimapInfo$Folder != "" &
      vapply(qualimapInfo$Folder, dir.exists, logical(1)),
    ,
    drop = FALSE
  ]
  if (nrow(qualimapInfo) == 0) {
    ezLog(
      "Skipping Qualimap multi-sample summary because no per-sample Qualimap directories are available.",
      level = "warn"
    )
    return(list(reportFile = NULL))
  }

  keyTable <- data.frame(
    SampleName = qualimapInfo$SampleName,
    Folder = qualimapInfo$Folder,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    keyTable,
    file = sampleKey,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  multiBamqcOk <- tryCatch(
    {
      ## JAVA_OPTS bypasses Qualimap's Java-8-only `-XX:MaxPermSize` flag so
      ## multi-bamqc also runs under the jdk21 on PATH (see get_dna_qualimap_stats).
      ezSystem(
        paste0(
          "unset DISPLAY; ",
          "JAVA_OPTS=\"-Xms32m -Xmx4G\" ",
          "qualimap multi-bamqc -d sampleKey.txt -outdir qualimap_MultiSample"
        )
      )
      TRUE
    },
    error = function(e) {
      ezLog(
        "Failed to build Qualimap multi-sample summary: ",
        conditionMessage(e),
        level = "warn"
      )
      FALSE
    }
  )
  if (!isTRUE(multiBamqcOk)) {
    return(list(reportFile = NULL))
  }

  reportFile <- file.path(
    "qualimap_MultiSample",
    "multisampleBamQcReport.html"
  )
  return(list(
    reportFile = if (file.exists(reportFile)) reportFile else NULL
  ))
}

##' @title Plot DNA library complexity
##' @description Plot the observed and extrapolated library complexity from an
##' ATACseqQC duplicate-frequency histogram.
##' @param histFile Duplicate-frequency histogram from `readsDupFreq()`.
##' @param times Number of bootstrap repetitions.
##' @param interpolate.sample.sizes Relative sample sizes used for interpolation.
##' @param extrapolate.sample.sizes Relative sample sizes used for extrapolation.
##' @param main Plot title.
##' @return Invisible data frame with estimated complexity values.
##' @template roxygen-template
ez_dna_lib_complexity <- function(
  histFile,
  times = 100,
  interpolate.sample.sizes = seq(0.1, 1, by = 0.1),
  extrapolate.sample.sizes = seq(5, 20, by = 5),
  main = c()
) {
  total <- histFile[, 1] %*% histFile[, 2]
  suppressWarnings({
    set.seed(42)
    result <- ds.rSAC.bootstrap(histFile, r = 1, times = times)
  })
  sequences <- c(interpolate.sample.sizes, extrapolate.sample.sizes)
  estimates <- data.frame(
    relative.size = sequences,
    values = rep(NA, length(sequences))
  )
  for (i in seq_along(sequences)) {
    estimates$values[i] <- result$f(sequences[i])
  }
  suppressWarnings({
    estimates$reads <- estimates$relative.size * total
  })
  plot(
    x = estimates$reads / 10^6,
    y = estimates$values / 10^6,
    type = "o",
    xlab = expression(
      Putative ~ sequenced ~ fragments ~ x ~ 10^6
    ),
    ylab = expression(
      Distinct ~ fragments ~ x ~ 10^6
    ),
    main = main
  )
  return(invisible(estimates))
}
environment(ez_dna_lib_complexity) <- asNamespace("ATACseqQC")

##' @title Collect per-sample DNA duplicate metrics from Picard
##' @description
##' Run Picard `MarkDuplicates` for one BAM file and parse the duplicate metrics
##' table needed for the cohort overview.
##' @param bamFile Full path to the BAM file.
##' @param sampleName Sample name used for output staging.
##' @param param ezRun parameter list.
##' @return Named list with parsed duplicate metrics.
##' @template roxygen-template
get_dna_picard_dup_stats <- function(
  bamFile,
  sampleName,
  param
) {
  metricFn <- paste0(sampleName, "_picardDupReport.txt")
  tempBam <- paste0(sampleName, "_picard_tmp.bam")
  cmd <- paste(
    prepareJavaTools("picard"),
    "MarkDuplicates",
    paste0("I=", bamFile),
    paste0("O=", tempBam),
    paste0("M=", metricFn),
    "REMOVE_DUPLICATES=false",
    paste0("OPTICAL_DUPLICATE_PIXEL_DISTANCE=", param$pixelDist),
    "> /dev/null"
  )
  ezSystem(cmd)
  ezSystem(paste("rm", tempBam))

  stopifnot(file.exists(metricFn))
  duplicateStats <- read.table(
    metricFn,
    skip = 6,
    nrows = 1,
    header = TRUE,
    sep = "\t"
  )
  allDuplicates <- if ("READ_PAIR_DUPLICATES" %in% colnames(duplicateStats)) {
    duplicateStats$READ_PAIR_DUPLICATES
  } else {
    NA_real_
  }
  optDuplicates <- if ("READ_PAIR_OPTICAL_DUPLICATES" %in% colnames(duplicateStats)) {
    duplicateStats$READ_PAIR_OPTICAL_DUPLICATES
  } else {
    NA_real_
  }
  readPairs <- if ("READ_PAIRS_EXAMINED" %in% colnames(duplicateStats)) {
    duplicateStats$READ_PAIRS_EXAMINED
  } else {
    NA_real_
  }

  result <- list(
    allDuplicates = allDuplicates,
    optDuplicates = optDuplicates,
    readPairs = readPairs,
    opticalDupRate = NA_real_,
    picardMetricsFile = metricFn
  )

  if (!is.na(result$optDuplicates) &&
    !is.na(result$allDuplicates) &&
    result$allDuplicates > 0) {
    result$opticalDupRate <- 100 * result$optDuplicates / result$allDuplicates
  }

  result$notes <- paste(
    "Picard duplicate metrics collected from",
    basename(metricFn)
  )

  return(result)
}

##' @title Parse a numeric value from Qualimap text output
##' @description Extract a numeric scalar from the first matching line.
##' @param lines Character vector read from `genome_results.txt`.
##' @param prefix Regular-expression prefix to strip from the matching line.
##' @param strip Character vector of literal tokens to remove after stripping
##' the prefix.
##' @return Numeric scalar or `NA_real_`.
##' @template roxygen-template
get_qualimap_value <- function(lines, prefix, strip = character()) {
  idx <- grep(prefix, lines)
  if (length(idx) == 0) {
    return(NA_real_)
  }
  value <- sub(prefix, "", lines[idx[1]])
  for (token in strip) {
    value <- gsub(token, "", value, fixed = TRUE)
  }
  suppressWarnings(as.numeric(value))
}

##' @title Parse a numeric value from Qualimap text output
##' @description Extract a numeric scalar from the first matching line.
##' @param lines Character vector read from `genome_results.txt`.
##' @param prefix Regular-expression prefix to strip from the matching line.
##' @return Numeric scalar or `NA_real_`.
##' @template roxygen-template
get_qualimap_numeric <- function(lines, prefix) {
  get_qualimap_value(lines = lines, prefix = prefix)
}

##' @title Parse a percent value from Qualimap text output
##' @description Extract a percent scalar from the first matching line.
##' @param lines Character vector read from `genome_results.txt`.
##' @param prefix Regular-expression prefix to strip from the matching line.
##' @return Numeric scalar or `NA_real_`.
##' @template roxygen-template
get_qualimap_percent <- function(lines, prefix) {
  get_qualimap_value(
    lines = lines,
    prefix = prefix,
    strip = "%"
  )
}

##' @title Parse a coverage value from Qualimap text output
##' @description Extract a coverage scalar from the first matching line.
##' @param lines Character vector read from `genome_results.txt`.
##' @param prefix Regular-expression prefix to strip from the matching line.
##' @return Numeric scalar or `NA_real_`.
##' @template roxygen-template
get_qualimap_coverage <- function(lines, prefix) {
  get_qualimap_value(
    lines = lines,
    prefix = prefix,
    strip = "X"
  )
}

##' @title Parse a read-count value from Qualimap text output
##' @description Extract a count scalar from the first matching line.
##' @param lines Character vector read from `genome_results.txt`.
##' @param prefix Regular-expression prefix to strip from the matching line.
##' @return Numeric scalar or `NA_real_`.
##' @template roxygen-template
get_qualimap_read_count <- function(lines, prefix) {
  get_qualimap_value(
    lines = lines,
    prefix = prefix,
    strip = c(",", "%")
  )
}

##' @title Parse coverage breadth from Qualimap text output
##' @description Extract the percent of reference covered at >= x-fold depth
##' from a `There is a N% of reference with a coverageData >= xX` line.
##' @param lines Character vector read from `genome_results.txt`.
##' @param xThreshold Integer coverage-depth threshold.
##' @return Numeric percent or `NA_real_`.
##' @template roxygen-template
get_dna_qualimap_coverage_breadth <- function(lines, xThreshold) {
  pattern <- paste0("coverageData >= ", xThreshold, "X")
  idx <- grep(pattern, lines, fixed = TRUE)
  if (length(idx) == 0) {
    return(NA_real_)
  }
  matched <- regmatches(
    lines[idx[1]],
    regexpr("[0-9.]+% of reference", lines[idx[1]])
  )
  if (length(matched) == 0) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(sub("% of reference", "", matched)))
}

##' @title Coerce a scalar legacy metric to numeric
##' @description Return a numeric scalar when possible, otherwise `NA_real_`.
##' @param x Value to coerce.
##' @return Numeric scalar or `NA_real_`.
##' @template roxygen-template
dna_bamstats_scalar_or_na <- function(x) {
  if (length(x) != 1 || is.null(x)) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(x))
}

##' @title Extract the first numeric token from a legacy text line
##' @description Parse the first numeric token, allowing comma separators.
##' @param x Character scalar containing a legacy metric line.
##' @return Numeric scalar or `NA_real_`.
##' @template roxygen-template
dna_bamstats_extract_first_number <- function(x) {
  if (length(x) != 1 || is.null(x) || is.na(x)) {
    return(NA_real_)
  }
  match <- regexpr("[0-9][0-9,]*\\.?[0-9]*", x)
  if (match[1] < 0) {
    return(NA_real_)
  }
  value <- regmatches(x, match)[1]
  value <- gsub(",", "", value, fixed = TRUE)
  suppressWarnings(as.numeric(value))
}

##' @title Extract a percentage from parentheses in legacy text
##' @description Parse a percentage stored as `(100%)` from a legacy line.
##' @param x Character scalar containing a legacy metric line.
##' @return Numeric scalar or `NA_real_`.
##' @template roxygen-template
dna_bamstats_extract_percent_in_parens <- function(x) {
  if (length(x) != 1 || is.null(x) || is.na(x)) {
    return(NA_real_)
  }
  match <- regexpr("\\(([0-9]+\\.?[0-9]*)%\\)", x, perl = TRUE)
  if (match[1] < 0) {
    return(NA_real_)
  }
  value <- regmatches(x, match)[1]
  value <- sub("^\\(", "", value)
  value <- sub("%\\)$", "", value)
  suppressWarnings(as.numeric(value))
}

##' @title Get one matching legacy Qualimap line
##' @description Return the first matching line or `NA_character_`.
##' @param lines Character vector read from `genome_results.txt`.
##' @param prefix Regular-expression prefix used to find the line.
##' @return Character scalar or `NA_character_`.
##' @template roxygen-template
get_dna_bamstats_legacy_qualimap_line <- function(lines, prefix) {
  idx <- grep(prefix, lines)
  if (length(idx) == 0) {
    return(NA_character_)
  }
  lines[idx[1]]
}

##' @title Read legacy DNA BAM stats parameters
##' @description Read a historical `parameters.tsv` and normalize the values to
##' the current DNA BAM stats parameter contract.
##' @param path Path to a legacy `parameters.tsv`.
##' @return Named parameter list compatible with `normalize_dna_bamstats_param`.
##' @template roxygen-template
read_dna_bamstats_legacy_parameters <- function(path) {
  paramTable <- ezRead.table(
    path,
    header = FALSE,
    row.names = NULL
  )
  param <- as.list(paramTable[[2]])
  names(param) <- paramTable[[1]]
  normalize_dna_bamstats_param(param)
}

##' @title Read a legacy DNA BAM stats input dataset
##' @description Convert a historical `input_dataset.tsv` into the minimal
##' dataset shape used by the maintained DNA BAM stats report.
##' @param path Path to a legacy `input_dataset.tsv`.
##' @return Data frame with one row per sample.
##' @template roxygen-template
read_dna_bamstats_legacy_input_dataset <- function(path) {
  inputDataset <- ezRead.table(
    path,
    row.names = NULL
  )
  data.frame(
    BAM = file.path("/srv/gstore/projects", inputDataset[["BAM [File]"]]),
    BAI = file.path("/srv/gstore/projects", inputDataset[["BAI [File]"]]),
    Species = inputDataset[["Species"]],
    refBuild = inputDataset[["refBuild"]],
    paired = tolower(inputDataset[["paired"]]) %in% c("true", "t", "1", "yes"),
    "Read Count" = suppressWarnings(as.numeric(inputDataset[["Read Count"]])),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = inputDataset[["Name"]]
  )
}

##' @title Parse legacy Picard duplicate metrics
##' @description Parse a historical Picard duplicate-metrics text file when it
##' exists, otherwise return `NA` values with the optional drilldown path.
##' @param metricFile Path to a legacy Picard metrics text file.
##' @param drilldownFile Optional fallback file path to keep as a drilldown link.
##' @return Named list of duplicate metrics.
##' @template roxygen-template
parse_dna_bamstats_legacy_picard_metrics <- function(
  metricFile = NULL,
  drilldownFile = NULL
) {
  if (is.null(metricFile) || !file.exists(metricFile)) {
    return(list(
      allDuplicates = NA_real_,
      optDuplicates = NA_real_,
      readPairs = NA_real_,
      opticalDupRate = NA_real_,
      picardMetricsFile = drilldownFile
    ))
  }
  duplicateStats <- read.table(
    metricFile,
    skip = 6,
    nrows = 1,
    header = TRUE,
    sep = "\t"
  )
  allDuplicates <- if ("READ_PAIR_DUPLICATES" %in% colnames(duplicateStats)) {
    duplicateStats$READ_PAIR_DUPLICATES
  } else {
    NA_real_
  }
  optDuplicates <- if ("READ_PAIR_OPTICAL_DUPLICATES" %in% colnames(duplicateStats)) {
    duplicateStats$READ_PAIR_OPTICAL_DUPLICATES
  } else {
    NA_real_
  }
  readPairs <- if ("READ_PAIRS_EXAMINED" %in% colnames(duplicateStats)) {
    duplicateStats$READ_PAIRS_EXAMINED
  } else {
    NA_real_
  }
  opticalDupRate <- NA_real_
  if (!is.na(optDuplicates) &&
    !is.na(allDuplicates) &&
    allDuplicates > 0) {
    opticalDupRate <- 100 *
      optDuplicates /
      allDuplicates
  }
  list(
    allDuplicates = allDuplicates,
    optDuplicates = optDuplicates,
    readPairs = readPairs,
    opticalDupRate = opticalDupRate,
    picardMetricsFile = metricFile
  )
}

##' @title Parse Qualimap metrics from legacy artifacts
##' @description Recover the current core DNA BAM stats metrics from a
##' historical `genome_results.txt`, including older file shapes that only
##' expose duplicate counts or mapping percentages indirectly.
##' @param qmFile Path to a legacy `genome_results.txt`.
##' @param nReads Total read count from the historical dataset metadata.
##' @return Named list with recovered Qualimap metrics.
##' @template roxygen-template
parse_dna_bamstats_legacy_qualimap_metrics <- function(
  qmFile,
  nReads = NA_real_
) {
  if (is.null(qmFile) || !file.exists(qmFile)) {
    return(list(
      mappingRate = NA_real_,
      dupRate = NA_real_,
      avgCoverage = NA_real_,
      errorRate = NA_real_,
      insertRate = NA_real_,
      delRate = NA_real_
    ))
  }
  allData <- readLines(qmFile)
  totalReadsLine <- get_dna_bamstats_legacy_qualimap_line(
    allData,
    "^.*number of reads = "
  )
  mappedReadsLine <- get_dna_bamstats_legacy_qualimap_line(
    allData,
    "^.*number of mapped reads = "
  )
  dupReadsLine <- get_dna_bamstats_legacy_qualimap_line(
    allData,
    "^.*number of duplicated reads \\((flagged|estimated)\\) = "
  )
  qualimapReads <- dna_bamstats_extract_first_number(totalReadsLine)
  mappedReads <- dna_bamstats_extract_first_number(mappedReadsLine)
  nReads <- dna_bamstats_scalar_or_na(nReads)
  mappingRate <- NA_real_
  if (!is.na(mappedReads) && !is.na(nReads) && nReads > 0) {
    mappingRate <- 100 * mappedReads / nReads
  } else {
    mappingRate <- dna_bamstats_extract_percent_in_parens(mappedReadsLine)
  }
  dupRate <- dna_bamstats_scalar_or_na(get_qualimap_percent(
    allData,
    "^.*duplication rate = "
  ))
  if (is.na(dupRate) && !is.na(qualimapReads) && qualimapReads > 0) {
    dupReads <- dna_bamstats_extract_first_number(dupReadsLine)
    if (!is.na(dupReads)) {
      dupRate <- 100 * dupReads / qualimapReads
    }
  }
  list(
    mappingRate = mappingRate,
    dupRate = dupRate,
    avgCoverage = dna_bamstats_scalar_or_na(get_qualimap_coverage(
      allData,
      "mean coverageData = "
    )),
    errorRate = dna_bamstats_scalar_or_na(get_qualimap_numeric(
      allData,
      "^.*general error rate = "
    )),
    insertRate = dna_bamstats_scalar_or_na(get_qualimap_percent(
      allData,
      "mapped reads with insertion percentage = "
    )),
    delRate = dna_bamstats_scalar_or_na(get_qualimap_percent(
      allData,
      "mapped reads with deletion percentage = "
    ))
  )
}

##' @title Copy one legacy directory tree
##' @description Copy a directory into a target location, replacing any
##' existing content at the target path.
##' @param fromDir Source directory.
##' @param toDir Target directory path.
##' @return Invisible target path.
##' @template roxygen-template
copy_dna_bamstats_legacy_directory <- function(fromDir, toDir) {
  dir.create(dirname(toDir), recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(toDir)) {
    unlink(toDir, recursive = TRUE, force = TRUE)
  }
  ok <- file.copy(fromDir, dirname(toDir), recursive = TRUE)
  stopifnot(ok)
  invisible(toDir)
}

##' @title Stage legacy DNABamStats artifacts
##' @description Copy the per-sample files needed to rerender an older
##' DNABamStats run with the maintained dataset-level report.
##' @param sourceRoot Root directory of the historical DNABamStats run.
##' @param targetLegacyRoot Target staging directory.
##' @param samples Character vector of sample names.
##' @return Invisible target directory.
##' @template roxygen-template
copy_dna_bamstats_legacy_sample_artifacts <- function(
  sourceRoot,
  targetLegacyRoot,
  samples
) {
  dir.create(targetLegacyRoot, recursive = TRUE, showWarnings = FALSE)
  for (sampleName in samples) {
    sampleDir <- file.path(sourceRoot, sampleName)
    if (dir.exists(sampleDir)) {
      ok <- file.copy(sampleDir, targetLegacyRoot, recursive = TRUE)
      stopifnot(ok)
    }
    sampleSamstat <- file.path(sourceRoot, paste0(sampleName, ".samstat.html"))
    if (file.exists(sampleSamstat)) {
      ok <- file.copy(sampleSamstat, targetLegacyRoot, overwrite = TRUE)
      stopifnot(ok)
    }
    samplePicardPdf <- file.path(sourceRoot, paste0(sampleName, ".picard.pdf"))
    if (file.exists(samplePicardPdf)) {
      ok <- file.copy(samplePicardPdf, targetLegacyRoot, overwrite = TRUE)
      stopifnot(ok)
    }
  }
  for (metaFile in c("dataset.tsv", "input_dataset.tsv", "parameters.tsv")) {
    metaPath <- file.path(sourceRoot, metaFile)
    if (file.exists(metaPath)) {
      ok <- file.copy(metaPath, targetLegacyRoot, overwrite = TRUE)
      stopifnot(ok)
    }
  }
  invisible(targetLegacyRoot)
}

##' @title Build one per-sample result from legacy artifacts
##' @description Assemble one current-style DNA BAM stats sample result from
##' historical Qualimap, Picard, and image artifacts.
##' @param sampleName Sample name.
##' @param dataset Legacy dataset converted to the current shape.
##' @param qualimapDir Directory containing `genome_results.txt`.
##' @param picardMetricsFile Optional Picard metrics text file.
##' @param picardDrilldownFile Optional Picard drilldown file.
##' @param fragmentSizePlot Optional fragment-size image file.
##' @param libComplexityPlot Optional library-complexity image file.
##' @param multiMatchRow Optional named vector of multi-match counts.
##' @param notes Character note describing the legacy import source.
##' @return Named list compatible with `build_dna_bamstats_report_data()`.
##' @template roxygen-template
build_dna_bamstats_result_from_legacy_sample <- function(
  sampleName,
  dataset,
  qualimapDir = NULL,
  picardMetricsFile = NULL,
  picardDrilldownFile = NULL,
  fragmentSizePlot = NULL,
  libComplexityPlot = NULL,
  multiMatchRow = NULL,
  notes = "Imported from legacy DNA BAM stats artifacts."
) {
  nReads <- dataset[sampleName, "Read Count"]
  qmFile <- if (!is.null(qualimapDir)) {
    file.path(qualimapDir, "genome_results.txt")
  } else {
    NULL
  }
  qualimapStats <- parse_dna_bamstats_legacy_qualimap_metrics(
    qmFile = qmFile,
    nReads = nReads
  )
  picardStats <- parse_dna_bamstats_legacy_picard_metrics(
    metricFile = picardMetricsFile,
    drilldownFile = picardDrilldownFile
  )

  multiMatchInFileTable <- NULL
  if (!is.null(multiMatchRow)) {
    multiMatchRow <- multiMatchRow[!is.na(multiMatchRow)]
    multiMatchInFileTable <- as.numeric(multiMatchRow)
    names(multiMatchInFileTable) <- names(multiMatchRow)
  }

  c(
    list(
      sampleName = sampleName,
      bamFile = dataset[sampleName, "BAM"],
      nReads = nReads,
      multiMatchInFileTable = multiMatchInFileTable,
      fragmentSizePlot = fragmentSizePlot,
      libComplexityPlot = libComplexityPlot,
      notes = notes,
      qualimapDir = qualimapDir
    ),
    qualimapStats,
    picardStats
  )
}

##' @title Build current DNA BAM stats report objects from legacy runs
##' @description Import a historical DNAQC or DNABamStats directory and convert
##' it into `dataset`, `param`, `resultList`, and `reportData` objects that the
##' maintained dataset-level report can render directly.
##' @param sourceRoot Root directory of the historical run or report directory.
##' @param inputDatasetPath Path to the historical `input_dataset.tsv`.
##' @param parametersPath Path to the historical `parameters.tsv`.
##' @param importMode Either `"dnaqc"` or `"dnabamstats"`.
##' @return Named list with `dataset`, `param`, `resultList`, `reportData`, and
##' `artifactRoot`.
##' @template roxygen-template
build_dna_bamstats_legacy_report_objects <- function(
  sourceRoot,
  inputDatasetPath,
  parametersPath,
  importMode = c("dnaqc", "dnabamstats")
) {
  importMode <- match.arg(importMode)
  dataset <- read_dna_bamstats_legacy_input_dataset(inputDatasetPath)
  param <- read_dna_bamstats_legacy_parameters(parametersPath)
  param$name <- if (importMode == "dnaqc") {
    "DNA_QC_Compatibility_ReRender"
  } else {
    "DNABamStats_Compatibility_ReRender"
  }

  resultList <- list()
  qualimapMulti <- list(reportFile = NULL)

  if (importMode == "dnaqc") {
    multiMatchTable <- ezRead.table(
      file.path(sourceRoot, "read-alignment-statistics.txt"),
      row.names = NULL
    )
    rownames(multiMatchTable) <- multiMatchTable$Sample
    multiMatchTable$Sample <- NULL

    for (sampleName in rownames(dataset)) {
      resultList[[sampleName]] <- build_dna_bamstats_result_from_legacy_sample(
        sampleName = sampleName,
        dataset = dataset,
        qualimapDir = file.path(sourceRoot, sampleName),
        picardMetricsFile = file.path(
          sourceRoot,
          paste0(sampleName, "_picardDupReport.txt")
        ),
        fragmentSizePlot = file.path(
          sourceRoot,
          paste0("fragmentSize_", sampleName, ".png")
        ),
        libComplexityPlot = file.path(
          sourceRoot,
          paste0("libComplexity_", sampleName, ".png")
        ),
        multiMatchRow = multiMatchTable[sampleName, , drop = TRUE],
        notes = "Imported from legacy DNAQC artifacts."
      )
    }
    qualimapMulti <- list(
      reportFile = file.path(
        sourceRoot,
        "qualimap_MultiSample",
        "multisampleBamQcReport.html"
      )
    )
  } else {
    for (sampleName in rownames(dataset)) {
      resultList[[sampleName]] <- build_dna_bamstats_result_from_legacy_sample(
        sampleName = sampleName,
        dataset = dataset,
        qualimapDir = file.path(sourceRoot, sampleName),
        picardDrilldownFile = file.path(
          sourceRoot,
          paste0(sampleName, ".picard.pdf")
        ),
        notes = "Imported from legacy DNABamStats artifacts."
      )
    }
  }

  reportData <- build_dna_bamstats_report_data(
    dataset = dataset,
    resultList = resultList
  )
  reportData$qualimapMulti <- qualimapMulti
  reportData$legacyAvailabilityNotes <- get_dna_bamstats_legacy_availability_notes(
    reportData = reportData,
    resultList = resultList,
    importMode = importMode
  )

  list(
    dataset = dataset,
    param = param,
    resultList = resultList,
    reportData = reportData,
    artifactRoot = sourceRoot
  )
}

##' @title Summarize legacy artifact limitations for the report
##' @description Create section-level explanatory notes when historical artifact
##' sets cannot support all metrics or drilldowns in the maintained report.
##' @param reportData Named report-data list created for a legacy import.
##' @param importMode Either `"dnaqc"` or `"dnabamstats"`.
##' @return Named list of character vectors keyed by report section.
##' @template roxygen-template
get_dna_bamstats_legacy_availability_notes <- function(
  reportData,
  resultList,
  importMode = c("dnaqc", "dnabamstats")
) {
  importMode <- match.arg(importMode)
  notes <- list(
    alignment = character(),
    duplicates = character(),
    resources = character()
  )

  fullSummaryTable <- reportData$fullSummaryTable
  resourceTable <- reportData$resourceTable

  if (importMode == "dnabamstats") {
    if (all(is.na(fullSummaryTable$OpticalDupRate))) {
      notes$duplicates <- c(
        notes$duplicates,
        paste(
          "OpticalDupRate and related Picard duplicate counts are unavailable",
          "for this legacy DNABamStats import because the historical run",
          "retained Picard PDF drilldowns but not parseable duplicate-metrics",
          "text files."
        )
      )
    }
    if (all(is.na(fullSummaryTable$MappingRate))) {
      notes$alignment <- c(
        notes$alignment,
        "MappingRate could not be recovered from the retained legacy Qualimap artifacts."
      )
    }
  }

  if (is.null(reportData$qualimapMulti$reportFile) ||
    is.na(reportData$qualimapMulti$reportFile) ||
    !file.exists(reportData$qualimapMulti$reportFile)) {
    notes$alignment <- c(
      notes$alignment,
      paste(
        "The multi-sample Qualimap summary is unavailable for this legacy",
        "import because the historical artifact set did not retain a usable",
        "multi-bamqc output."
      )
    )
  }

  if (all(vapply(resultList, function(x) {
    is.null(x$multiMatchInFileTable)
  }, logical(1)))) {
    notes$alignment <- c(
      notes$alignment,
      paste(
        "Multi-matching counts are unavailable for this legacy import because",
        "the historical artifact set did not retain the dataset-level",
        "alignment summary used by the maintained report."
      )
    )
  }

  missingFrag <- sum(is.na(resourceTable$FragmentSizePlot) | resourceTable$FragmentSizePlot == "")
  if (missingFrag > 0) {
    notes$resources <- c(
      notes$resources,
      paste(
        "Fragment-size plots are unavailable for",
        missingFrag,
        "sample(s) because the corresponding legacy image artifacts were not retained."
      )
    )
  }

  missingLib <- sum(is.na(resourceTable$LibraryComplexityPlot) | resourceTable$LibraryComplexityPlot == "")
  if (missingLib > 0) {
    notes$resources <- c(
      notes$resources,
      paste(
        "Library-complexity plots are unavailable for",
        missingLib,
        "sample(s) because the corresponding legacy image artifacts were not retained."
      )
    )
  }

  notes
}

##' @title Render the maintained DNA BAM stats report from legacy artifacts
##' @description Stage a historical DNAQC or DNABamStats run into a writable
##' directory and render the current dataset-level report from the recovered
##' objects.
##' @param sourceRoot Root directory of the historical run or report directory.
##' @param inputDatasetPath Path to the historical `input_dataset.tsv`.
##' @param parametersPath Path to the historical `parameters.tsv`.
##' @param targetDir Output directory for the rerendered report.
##' @param importMode Either `"dnaqc"` or `"dnabamstats"`.
##' @param htmlFile Name of the rendered HTML file.
##' @param templateFile Optional template path overriding the installed package
##' template.
##' @return Invisible output directory path.
##' @template roxygen-template
render_dna_bamstats_legacy_report <- function(
  sourceRoot,
  inputDatasetPath,
  parametersPath,
  targetDir,
  importMode = c("dnaqc", "dnabamstats"),
  htmlFile = "00index.html",
  templateFile = NULL
) {
  importMode <- match.arg(importMode)

  dir.create(dirname(targetDir), recursive = TRUE, showWarnings = FALSE)
  unlink(targetDir, recursive = TRUE, force = TRUE)
  dir.create(targetDir, recursive = TRUE, showWarnings = FALSE)

  stagedSourceRoot <- if (importMode == "dnaqc") {
    file.path(targetDir, basename(sourceRoot))
  } else {
    file.path(targetDir, "legacy_artifacts")
  }
  if (importMode == "dnaqc") {
    copy_dna_bamstats_legacy_directory(sourceRoot, stagedSourceRoot)
  } else {
    dataset <- read_dna_bamstats_legacy_input_dataset(inputDatasetPath)
    copy_dna_bamstats_legacy_sample_artifacts(
      sourceRoot = sourceRoot,
      targetLegacyRoot = stagedSourceRoot,
      samples = rownames(dataset)
    )
  }

  legacyObjects <- build_dna_bamstats_legacy_report_objects(
    sourceRoot = stagedSourceRoot,
    inputDatasetPath = inputDatasetPath,
    parametersPath = parametersPath,
    importMode = importMode
  )

  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)
  setwd(targetDir)

  write_dna_bamstats_support_files(
    reportData = legacyObjects$reportData
  )

  if (is.null(templateFile)) {
    makeRmdReport(
      dataset = legacyObjects$dataset,
      param = legacyObjects$param,
      resultList = legacyObjects$resultList,
      reportData = legacyObjects$reportData,
      htmlFile = htmlFile,
      rmdFile = "DnaBamStats.Rmd",
      reportTitle = "DNA BAM Stats",
      selfContained = TRUE
    )
  } else {
    file.copy(templateFile, targetDir, overwrite = TRUE)
    renderEnv <- new.env(parent = globalenv())
    renderEnv$dataset <- legacyObjects$dataset
    renderEnv$param <- legacyObjects$param
    renderEnv$resultList <- legacyObjects$resultList
    renderEnv$reportData <- legacyObjects$reportData
    renderEnv$alignmentCountBarPlot <- alignmentCountBarPlot
    rmarkdown::render(
      input = file.path(targetDir, basename(templateFile)),
      envir = renderEnv,
      output_dir = targetDir,
      output_file = htmlFile,
      quiet = TRUE
    )
  }

  invisible(targetDir)
}

##' @template app-template
##' @templateVar method ezMethodDnaBamStats(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
EzAppDnaBamStats <-
  setRefClass(
    "EzAppDnaBamStats",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodDnaBamStats
        name <<- "EzAppDnaBamStats"
        appDefaults <<- rbind(
          pixelDist = ezFrame(
            Type = "integer",
            DefaultValue = 2500,
            Description = "optical duplicate pixel distance for Picard-based duplicate metrics"
          ),
          runQualimap = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "run Qualimap bamqc to derive mapping, coverage, duplication, and error-rate metrics"
          ),
          runPicard = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "run Picard MarkDuplicates to derive duplicate and optical-duplicate metrics"
          )
        )
      }
    )
  )
