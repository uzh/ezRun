###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @describeIn plotBamStat Gets the table containing the counts of each type.
getTypeCountTable = function(resultList, name) {
  tbl = data.frame(row.names = rownames(resultList[[1]][[name]]))
  for (sm in names(resultList)) {
    counts = resultList[[sm]][[name]]
    tbl[sm] = signif(counts[rownames(tbl), "count"] / 1e6, digits = 4)
  }
  return(tbl)
}

##' @describeIn plotBamStat Gets the table containing the coverage of each type.
getTypeCoverageTable = function(resultList, name) {
  tbl = data.frame(row.names = rownames(resultList[[1]][[name]]))
  for (sm in names(resultList)) {
    counts = resultList[[sm]][[name]]
    tbl[sm] = counts[rownames(tbl), "count"] / counts[rownames(tbl), "width"]
  }
  return(tbl)
}

##' @describeIn plotBamStat Assigns the rows of a type count or coverage table to the
##' sections feature types (1), mRNA regions (2), unannotated (3) and chromosomes (4).
get_type_table_sections <- function(rowNames, seqNames) {
  ## the mRNA rows are created in getTargetTypeCounts()
  mrnaRows <- c(
    "mRNA Exons",
    "mRNA Introns",
    "mRNA Promoter 2kb",
    "mRNA Downstream 2kb"
  )
  sections <- rep(1L, length(rowNames))
  sections[rowNames %in% mrnaRows] <- 2L
  sections[rowNames == "unannotated"] <- 3L
  sections[rowNames %in% seqNames] <- 4L
  sections
}

##' @describeIn plotBamStat Renders a type count or coverage table as html, with the rows
##' ordered into sections (feature types, mRNA regions, unannotated, chromosomes) and
##' a thick line between the sections.
kable_type_table <- function(tbl, seqNames, caption, footnoteText) {
  sections <- get_type_table_sections(rownames(tbl), seqNames)
  rowOrder <- order(sections)
  tbl <- tbl[rowOrder, , drop = FALSE]
  sectionEnds <- which(diff(sections[rowOrder]) != 0)
  k <- knitr::kable(tbl, row.names = TRUE, format = "html", caption = caption) |>
    kableExtra::kable_styling(
      bootstrap_options = "striped",
      full_width = FALSE,
      position = "left"
    )
  if (length(sectionEnds) > 0) {
    k <- kableExtra::row_spec(
      k,
      sectionEnds,
      extra_css = "border-bottom: 3px solid #333;"
    )
  }
  kableExtra::footnote(k, general = footnoteText)
}

##' @describeIn plotBamStat Prints the table of the rows \code{rowsUse} and, collapsed
##' below it, the table of the low-abundance rows \code{rowsLow}.
print_type_tables <- function(
  tbl,
  rowsUse,
  rowsLow,
  seqNames,
  digits,
  caption,
  footnoteText,
  lowSummary
) {
  print(kable_type_table(
    signif(tbl[rowsUse, , drop = FALSE], digits = digits),
    seqNames,
    caption = caption,
    footnoteText = footnoteText
  ))
  if (length(rowsLow) > 0) {
    cat("\n<details>\n<summary>", lowSummary, "</summary>\n\n", sep = "")
    print(kable_type_table(
      signif(tbl[rowsLow, , drop = FALSE], digits = digits),
      seqNames,
      caption = paste(caption, "of low-abundance features"),
      footnoteText = footnoteText
    ))
    cat("\n</details>\n\n")
  }
  invisible(NULL)
}

alignmentCountBarPlot <- function(
  mmCounts,
  relative = FALSE,
  file = NULL,
  multiCountColors = NULL
) {
  library(plotly)
  title <- ifelse(relative, "alignment proportions", "total alignments")
  multiCount = as.integer(colnames(mmCounts))
  isSmall = multiCount <= 3
  if (any(!isSmall)) {
    mmCounts = cbind(
      mmCounts[, isSmall, drop = FALSE],
      ">3" = rowSums(mmCounts[, !isSmall, drop = FALSE])
    )
  }

  if (!is.null(file)) {
    ezWrite.table(mmCounts, file = file)
  }

  colnames(mmCounts) = paste(colnames(mmCounts), "hit(s)")
  if (is.null(multiCountColors)) {
    multiCountColors = c(
      "0 hit(s)" = "gray",
      "1 hit(s)" = "blue",
      "2 hit(s)" = "cyan",
      "3 hit(s)" = "green",
      ">3 hit(s)" = "orange"
    )
  }
  stopifnot(colnames(mmCounts) %in% names(multiCountColors))

  x = mmCounts[, rev(colnames(mmCounts)), drop = F]
  if (isTRUE(relative)) {
    x <- sweep(x, MARGIN = 1, STATS = rowSums(x), FUN = "/")
  }
  x <- as.data.frame(x)
  x$"sample" <- rownames(x)
  m <- list(
    l = 80,
    r = 80,
    b = 200,
    t = 100,
    pad = 0
  )

  p <- plot_ly(
    x,
    x = ~sample,
    y = as.formula(paste0("~`", colnames(x)[1], "`")),
    type = "bar",
    name = colnames(x)[1],
    marker = list(color = multiCountColors[colnames(x)[1]])
  )
  for (i in 2:(ncol(x) - 1)) {
    p <- p %>%
      add_trace(
        y = as.formula(paste0("~`", colnames(x)[i], "`")),
        name = colnames(x)[i],
        marker = list(color = multiCountColors[colnames(x)[i]])
      )
  }
  p <- p %>%
    plotly::layout(
      yaxis = list(title = ifelse(relative, "Percentage of reads", "Count")),
      barmode = 'stack',
      title = title,
      margin = m
    )
  p
}

##' @describeIn plotBamStat Plots the position specific error rates.
plotPosSpecificErrorRate = function(
  errorRate,
  main = "Per base mismatch rate"
) {
  par(mfrow = c(1, 2))
  trimmedClippedRateMatrix = rbind(
    errorRate[["trimmedRate"]],
    errorRate[["clippedRate"]]
  )
  errorRateMatrix = errorRate[["errorRate"]]
  colnames(trimmedClippedRateMatrix) = 1:ncol(trimmedClippedRateMatrix)
  barplot(
    trimmedClippedRateMatrix,
    col = c("darkorange", "blue"),
    space = 0,
    border = NA,
    xlab = "Base Position in Read",
    ylab = "Trimmed-Clipped Rate",
    main = paste(main, "trimmed-clipped rate"),
    ylim = c(0, 1), #max(0.1, 1.2*max(colSums(trimmedClippedRateMatrix)))),
    xpd = FALSE
  )
  legend(
    "top",
    c("trimmed", "clipped"),
    col = c("darkorange", "blue"),
    pch = 20
  )
  yMax = max(0.1, max(errorRateMatrix, na.rm = TRUE))
  barCol = ifelse(yMax > 0.1, "red", "gray")
  barplot(
    errorRateMatrix,
    space = 0,
    border = NA,
    xlab = "Base Position in Read",
    ylab = "Mismatch Rate",
    main = paste(main, "mismatch rate"),
    ylim = c(0, yMax),
    col = barCol,
    xpd = FALSE
  )
  if (yMax > 0.1) {
    warning("The error rate range is larger than 0.1")
  }
}
