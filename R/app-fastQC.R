###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodFastQC <- function(input = NA, output = NA, param = NA,
                           htmlFile = "00index.html") {
  setwdNew(basename(output$getColumn("Report")))

  # Preprocessing
  isUBam <- input$readType() == "bam"
  if (isTRUE(isUBam)) {
    if (isTRUE(param$perLibrary)) {
      fastqInput <- ezMethodBam2Fastq(
        input = input, param = param,
        OUTPUT_PER_RG = FALSE
      )
    } else {
      ## We only support one uBam when it's per cell mode
      stopifnot(input$getLength() == 1L)
      fastqInput <- ezMethodBam2Fastq(
        input = input, param = param,
        OUTPUT_PER_RG = TRUE
      )
    }
    input <- fastqInput$copy()
  }
  dataset <- input$meta
  samples <- rownames(dataset)

  ans4Report <- list() # a list of results for rmarkdown report
  ans4Report[["dataset"]] <- dataset

  if (!is.null(dataset$"Read Count")) {
    readCount <- signif(dataset$"Read Count" / 1e6, digits = 3)
    names(readCount) <- rownames(dataset)
  } else {
    stop("'Read Count' has to be specified in the dataset.")
  }
  ans4Report[["Read Counts"]] <- readCount

  if (sum(dataset$`Read Count`) > 1e9) {
    doSubsample <- TRUE # subsample to 1Mio reads
    input <- ezMethodSubsampleFastq(input = input, param = param)
    dataset <- input$meta
  } else {
    doSubsample <- FALSE
  }

  files <- c()
  for (sm in samples) {
    files[paste0(sm, "_R1")] <- input$getFullPaths("Read1")[sm]
    if (isTRUE(param$paired)) {
      files[paste0(sm, "_R2")] <- input$getFullPaths("Read2")[sm]
    }
  }
  nFiles <- length(files)

  ## guess the names of the report directories that will be creatd by fastqc
  reportDirs <- sub("\\.(fastq|fq|bam)(\\.gz)*$", "_fastqc", basename(files))
  stopifnot(!any(duplicated(reportDirs)))
  filesUse <- files[!file.exists(reportDirs)]
  if (length(filesUse) > 0) {
    cmd <- paste(
      "fastqc", "--extract -o . -t", min(param$cores, 8),
      "-a", FASTQC_ADAPTERS, "--kmers 7",
      param$cmdOptions, paste(filesUse, collapse = " "),
      "> fastqc.out", "2> fastqc.err"
    )
    memory.size(max = F)
    gc()
    memory.size(max = F)
    result <- ezSystem(cmd)
  }
  statusToPng <- c(PASS = "tick.png", WARN = "warning.png", FAIL = "error.png")

  ## Copy the style files and templates
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c(
      "fgcz.css", "FastQC.Rmd", "FastQC_overview.Rmd",
      "fgcz_header.html", "banner.png"
    )
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)

  ## collect the overview table
  plots <- c(
    "Per base sequence quality" = "per_base_quality.png",
    "Per sequence quality scores" = "per_sequence_quality.png",
    "Per tile sequence quality" = "per_tile_quality.png",
    "Per base sequence content" = "per_base_sequence_content.png",
    "Per sequence GC content" = "per_sequence_gc_content.png",
    "Per base N content" = "per_base_n_content.png",
    "Sequence Length Distribution" = "sequence_length_distribution.png",
    "Sequence Duplication Levels" = "duplication_levels.png",
    "Adapter Content" = "adapter_content.png",
    "Kmer Content" = "kmer_profiles.png"
  )

  ## make for each plot type an html report with all samples
  plotPages <- sub(".png", ".html", plots)
  for (i in 1:length(plots)) {
    plotPage <- plotPages[i]
    pngs <- file.path(reportDirs, "Images", plots[i])
    rmarkdown::render(
      input = "FastQC_overview.Rmd", envir = new.env(),
      output_dir = ".", output_file = plotPage, quiet = TRUE
    )
  }

  ## Each sample can have different number of reports.
  ## Especially per tile sequence quality
  nrReports <- sapply(
    reportDirs,
    function(x) {
      smy <- ezRead.table(file.path(x, "summary.txt"),
        row.names = NULL, header = FALSE
      )
      nrow(smy)
    }
  )
  i <- which.max(nrReports)
  smy <- ezRead.table(file.path(reportDirs[i], "summary.txt"),
    row.names = NULL, header = FALSE
  )
  rowNames <- paste0(
    "<a href=", reportDirs, "/fastqc_report.html>",
    names(files), "</a>"
  )

  tbl <- ezMatrix("", rows = rowNames, cols = smy[[2]])
  for (i in 1:nFiles) {
    smy <- ezRead.table(file.path(reportDirs[i], "summary.txt"),
      row.names = NULL, header = FALSE
    )

    href <- paste0(
      reportDirs[i], "/fastqc_report.html#M",
      0:(ncol(tbl) - 1)
    )[colnames(tbl) %in% smy[[2]]]
    img <- paste0(reportDirs[1], "/Icons/", statusToPng[smy[[1]]])
    tbl[i, colnames(tbl) %in% smy[[2]]] <- paste0(
      "<a href=", href,
      "><img src=", img, "></a>"
    )
  }
  colnames(tbl) <- ifelse(colnames(tbl) %in% names(plotPages),
    paste0(
      "<a href=", plotPages[colnames(tbl)],
      ">", colnames(tbl),
      "</a>"
    ),
    colnames(tbl)
  )

  ans4Report[["Fastqc quality measures"]] <- tbl

  qualMatrixList <- ezMclapply(files, getQualityMatrix, mc.cores = param$cores)
  ans4Report[["Per Base Read Quality"]] <- qualMatrixList

  ## debug
  ## save(ans4Report, file="ans4Report.rda")

  ## generate the main reports
  rmarkdown::render(
    input = "FastQC.Rmd", envir = new.env(),
    output_dir = ".", output_file = htmlFile, quiet = TRUE
  )

  prepareRmdLib()

  ## generate multiQC report
  ezSystem("multiqc .")

  ## Cleaning
  if (isTRUE(isUBam) || isTRUE(doSubsample)) {
    file.remove(files)
  }
  unlink(paste0(reportDirs, ".zip"), recursive = TRUE)

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFastQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
##' @section Functions:
##' \itemize{
##'   \item{\code{plotReadCountToLibConc(dataset, colname): }}
##'   {Plots \code{colname} from \code{dataset} against read counts in millions.}
##'   \item{\code{getQualityMatrix(inputFile): }}
##'   {Gets a quality count matrix from a fastq or gziped fastq.gz file with dimensions read quality and read length.}
##'   \item{\code{plotQualityMatrixAsHeatmap(qualMatrixList, isR2=FALSE, xScale=1, yScale=1): }}
##'   {Returns a png table of quality matrices interpreted as heatmaps.}
##'   \item{\code{plotQualityHeatmap(result, name=NULL, colorRange=c(0,sqrt(40)), colors=gray((1:256)/256), main=NULL, pngFileName=NULL, xScale=1, yScale=1): }}
##'   {Creates and returns the images used by \code{plotQualityMatrixAsHeatmap()}.}
##' }
EzAppFastqc <-
  setRefClass("EzAppFastqc",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodFastQC
        name <<- "EzAppFastqc"
        appDefaults <<- rbind(perLibrary = ezFrame(Type = "logical", DefaultValue = TRUE, Description = "Run FastQC per library or per cell for single cell experiment"))
      }
    )
  )

plotReadCountToLibConc <- function(dataset, colname) {
  if (colname %in% colnames(dataset) && nrow(dataset) > 1) {
    if (!all(dataset[[colname]] == 0) && !any(is.na(dataset[[colname]]))) {
      ## LibCon column can all be 0 or NA. Then don't plot.
      dataset <- dataset[order(dataset$"Read Count", decreasing = T), ]
      dataset$"Read Count" <- dataset$"Read Count" / 10^6
      corResult <- cor.test(dataset$"Read Count", dataset[[colname]],
        method = "spearman"
      )
      regressionResult <- lm(dataset[[colname]] ~ dataset$"Read Count")
      label <- sub(" \\[.*", "", colname)

      ## plotly
      require(plotly)
      # a function to calculate your abline
      xmin <- min(dataset$"Read Count") - 5
      xmax <- max(dataset$"Read Count") + 5
      intercept <- regressionResult$coefficients[1]
      slope <- regressionResult$coefficients[2]
      p_abline <- function(x, a, b) {
        y <- a * x + b
        return(y)
      }

      p <- plot_ly(
        x = dataset$"Read Count", y = dataset[[colname]],
        text = rownames(dataset)
      ) %>%
        add_markers() %>%
        add_text(textposition = "top right") %>%
        plotly::layout(showlegend = FALSE)
      a <- list(
        x = max(dataset$"Read Count"),
        y = max(dataset[[colname]]),
        text = paste0("r=", round(corResult$estimate, 2)),
        xref = "x",
        yref = "y",
        showarrow = FALSE
      )
      p <- p %>% plotly::layout(
        shapes = list(
          type = "line", line = list(dash = "dash"),
          x0 = xmin, x1 = xmax,
          y0 = p_abline(xmin, slope, intercept),
          y1 = p_abline(xmax, slope, intercept)
        ),
        annotations = a,
        title = label, yaxis = list(title = label),
        xaxis = list(title = "Counts [Mio]")
      )
      return(p)
    }
  }
}

plotQualityMatrixAsHeatmapGG2 <- function(qualMatrixList, isR2 = FALSE,
                                          xScale = 1, yScale = 1) {
  colorsGray <- gray((30:256) / 256)
  minPercent <- 0
  maxPercent <- sqrt(40)
  minDiff <- -5
  maxDiff <- 5
  plotList <- list() # to store all the ggplot objects

  ## test if R2 exists
  index <- list("R1" = which(!isR2))

  if (any(isR2)) {
    stopifnot(sum(isR2) == sum(!isR2))
    index[["R2"]] <- which(isR2)
  }
  for (nm in names(index)) {
    plotList[[nm]] <- list()
    idx <- index[[nm]]
    ## Plot the color key for the average quality heatmap R1_1
    # colorKeyFile = paste0("averageReadsQuality-Key_", nm, ".png")
    by.label <- 1
    at <- seq(from = minPercent, to = maxPercent, by = by.label)
    p <- ezColorLegendGG2(
      colorRange = c(minPercent, maxPercent),
      colors = colorsGray, vertical = FALSE, by.label = by.label,
      at = at, labels = as.character(at^2)
    )
    plotList[[nm]][["Avg Qual Colors"]] <- p

    # result = ezMatrix(0, dim=dim(qualMatrixList[[idx[1]]]))
    result <- ezMatrix(0, dim = apply(sapply(qualMatrixList[idx], dim), 1, max))
    resultCount <- result
    for (i in idx) {
      qm <- qualMatrixList[[i]]
      result[1:nrow(qm), 1:ncol(qm)] <- result[1:nrow(qm), 1:ncol(qm)] + qm
      resultCount[1:nrow(qm), 1:ncol(qm)] <- resultCount[1:nrow(qm), 1:ncol(qm)] + 1
    }
    result <- result / resultCount
    ## The hard way to deal with NaN in result
    result <- sweep(result,
      MARGIN = 2,
      STATS = colSums(result, na.rm = TRUE), FUN = "/"
    )
    avgQual <- signif(result * 100, digits = 3)
    p <- plotQualityHeatmapGG2(
      result = sqrt(avgQual),
      colorRange = c(minPercent, maxPercent),
      colors = colorsGray,
      main = paste("averageReadsQuality", nm, sep = "_"),
      xScale = xScale, yScale = yScale
    )
    plotList[[nm]][["Average"]] <- p

    ## plot the difference quality heatmap for R1_1
    at <- seq(from = minDiff, to = maxDiff, by = by.label)
    p <- ezColorLegendGG2(
      colorRange = c(minDiff, maxDiff),
      colors = getBlueRedScale(), vertical = FALSE,
      by.label = by.label,
      at = at, labels = as.character(at)
    )
    plotList[[nm]][["Diff Qual Colors"]] <- p

    for (sampleName in names(qualMatrixList[idx])) {
      qm <- qualMatrixList[[sampleName]]
      qm <- sweep(qm,
        MARGIN = 2,
        STATS = colSums(qm, na.rm = TRUE), FUN = "/"
      )
      diffResult <- signif(qm * 100, digits = 3) - avgQual[1:nrow(qm), 1:ncol(qm)]
      p <- plotQualityHeatmapGG2(diffResult,
        colorRange = c(minDiff, maxDiff),
        colors = getBlueRedScale(),
        main = paste("diffReadsQuality", sampleName,
          sep = "_"
        ),
        xScale = xScale, yScale = yScale
      )
      plotList[[nm]][[sampleName]] <- p
    }
  }
  return(plotList)
}

plotQualityHeatmapGG2 <- function(result, name = NULL, colorRange = c(0, sqrt(40)),
                                  colors = gray((1:256) / 256), main = NULL,
                                  xScale = 1, yScale = 1) {
  require(reshape2)
  ## some ugly controls of labels
  labCol <- seq(0, ncol(result), by = 10)
  labCol[1] <- 1
  labRow <- seq(0, nrow(result) - 1, by = 5)

  result[result > colorRange[2]] <- colorRange[2]
  result[result < colorRange[1]] <- colorRange[1]
  toPlot <- melt(result)
  p <- ggplot(data = toPlot, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_gradientn(colours = colors, limits = colorRange) +
    theme_bw() +
    scale_y_continuous(
      breaks = seq(1, nrow(result), by = 5), labels = labRow,
      expand = c(0, 0)
    ) +
    scale_x_continuous(breaks = labCol, labels = labCol, expand = c(0, 0)) +
    theme(
      panel.border = element_blank(), panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5)
    ) +
    xlab("Read Position") +
    ylab("Read Quality") +
    ggtitle(main)
  return(p)
}

### Sample up to 300k reads from a fastq file or bam file.
### Calculate the quality matrix
getQualityMatrix <- function(fn) {
  ## This implementation is faster than the FastqStreamer.
  require(ShortRead)
  nReads <- 3e5

  if (grepl("\\.bam$", fn)) {
    cmd <- paste("samtools flagstat", fn)
    cmdOutput <- ezSystem(cmd, intern = TRUE)
    nrTotal <- eval(parse(text = sub(" in total.*", "", cmdOutput[1])))

    tempBamFn <- paste(Sys.getpid(), "temp.bam", sep = "-")
    on.exit(file.remove(tempBamFn), add = TRUE)
    cmd <- paste(
      "samtools view -s", nReads / nrTotal, "-b", fn,
      ">", tempBamFn
    )
    ezSystem(cmd)
    tempFastqFn <- paste(Sys.getpid(), "temp.fastq", sep = "-")
    on.exit(file.remove(tempFastqFn), add = TRUE)
    bam2fastq(
      bamFn = tempBamFn, OUTPUT_PER_RG = FALSE,
      fastqFns = tempFastqFn, paired = FALSE
    )
    reads <- readFastq(tempFastqFn)
  } else {
    f <- FastqSampler(fn, nReads) ## we sample no more than 300k reads.
    reads <- yield(f)
  }
  qual <- quality(reads)

  maxReadLength <- max(width(qual))
  qualMatrix <- as(qual, "matrix")
  maxQuality <- max(qualMatrix, na.rm = TRUE)

  qualCountMatrix <- ezMatrix(0, rows = 0:maxQuality, cols = 1:maxReadLength)
  for (basePos in 1:ncol(qualMatrix)) {
    qualCountByPos <- table(qualMatrix[, basePos])
    qualCountMatrix[names(qualCountByPos), basePos] <- qualCountByPos
  }
  return(qualCountMatrix)
}

plateStatistics <- function(dataset,
                            colname = c(
                              "Read Count",
                              "LibConc_qPCR [Characteristic]",
                              "LibConc_100_800bp [Characteristic]"
                            )) {
  colsExist <- colname %in% colnames(dataset)
  if (any(!colsExist)) {
    warning("No column ", colname[!colsExist], " in dataset!")
    colname <- colname[colsExist]
  }
  colsNumeric <- sapply(dataset[, colname, drop = FALSE], is, "numeric")
  if (any(!colsNumeric)) {
    warning("The column ", colname[!colsNumeric], " is non-numeric.")
    colname <- colname[colsNumeric]
  }
  colsNA <- is.na(colSums(dataset[, colname, drop = FALSE]))
  if (any(colsNA)) {
    message("The column ", colname[colsNA], " has NA!")
    colname <- colname[!colsNA]
  }
  colsZero <- colSums(dataset[, colname, drop = FALSE]) == 0
  if (any(colsZero)) {
    message("The column ", colname[colsZero], " is empty!")
    colname <- colname[!colsZero]
  }
  if (length(colname) == 0L) {
    warning("No suitable columns left in dataset!")
    return(NA)
  }

  if (!is.null(dataset$`PlatePosition [Characteristic]`) &&
    !any(is.na(dataset$`PlatePosition [Characteristic]`))) {
    # `PlatePosition [Characteristic]` may not exist or is empty

    plateChar <- dataset$`PlatePosition [Characteristic]`
    ## Plate position should be in the format of *_C4;
    ## * is the plate number
    ## C is the column; 4 is the row
    require(stringr, quietly = TRUE)
    plateNumber <- sub("_[[:alpha:]]\\d+$", "", plateChar)
    datasetByPlate <- split(dataset, plateNumber)
    ans <- list()
    for (i in seq_len(length(datasetByPlate))) {
      platePos <- str_extract(
        datasetByPlate[[i]]$`PlatePosition [Characteristic]`,
        "_[[:alpha:]]\\d+$"
      )
      if (any(is.na(platePos))) {
        warning("The PlatePosition format is not supported!")
        return(NA)
      }
      plateRow <- str_extract(platePos, "[[:alpha:]]")
      plateCol <- as.numeric(str_extract(platePos, "\\d+"))
      ans[[names(datasetByPlate)[i]]] <- list()
      for (oneCol in colname) {
        ## always the entire plate should be shown which is either 8x12 or 16x24 ....
        counts <- datasetByPlate[[i]][[oneCol]]
        countMatrix <- ezMatrix(NA,
          rows = LETTERS[1:ifelse(max(plateRow) > "I", 16, 8)],
          cols = seq_len(ifelse(max(as.integer(plateCol)) > 12, 24, 12))
        )
        for (j in seq_len(length(counts))) {
          countMatrix[plateRow[j], plateCol[j]] <- counts[j]
        }
        ans[[names(datasetByPlate)[i]]][[oneCol]] <- countMatrix
      }
    }
    return(ans)
  } else {
    warning("PlatePosition [Characteristic] information is not available!")
    return(NA)
  }
}

heatmapPlate <- function(x, title = "unnamed", center = TRUE, log10 = TRUE, ...) {
  require(plotly)
  ## do not plot if there are only NA or zeros
  if (all(x %in% c(NA, 0))) {
    return(NULL)
  }
  if (isTRUE(log10)) {
    ## shift zeros a bit
    isZero <- x == 0
    isZero[is.na(isZero)] <- FALSE
    x[isZero] <- min(0.25 * x[x > 0], na.rm = TRUE)

    x <- log10(x)
    if (isTRUE(center)) {
      medianX <- median(x, na.rm = TRUE)
      p <- plot_ly(
        z = x, x = colnames(x), y = rownames(x), type = "heatmap",
        zmin = medianX - log10(2), zmax = medianX + log10(2),
        hoverinfo = "text",
        text = matrix(paste0(
          "10^", format(x, digits = 3), "=",
          10^x
        ), ncol = ncol(x)),
        ...
      )
    } else {
      p <- plot_ly(
        z = x, x = colnames(x), y = rownames(x), type = "heatmap",
        hoverinfo = "text",
        text = matrix(paste0(
          "10^", format(x, digits = 3), "=",
          10^x
        ), ncol = ncol(x)),
        ...
      )
    }
  } else {
    if (isTRUE(center)) {
      medianX <- median(x, na.rm = TRUE)
      p <- plot_ly(
        z = x, x = colnames(x), y = rownames(x), type = "heatmap",
        zmin = 0, zmax = medianX * 2, ...
      )
    } else {
      p <- plot_ly(z = x, x = colnames(x), y = rownames(x), type = "heatmap", ...)
    }
  }

  p <- p %>% plotly::layout(
    xaxis = list(autotick = FALSE, dtick = 1),
    yaxis = list(autorange = "reversed"),
    title = title
  )
  p
}
