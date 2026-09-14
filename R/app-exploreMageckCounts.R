###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodMageckCountQC <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  setwdNew(basename(output$getColumn("Report")))
  rawData <- loadMageckCountDataset(input, param)
  if (isError(rawData)) {
    writeErrorReport(htmlFile, param = param, error = rawData$error)
    return("Error")
  }
  metadata(rawData)$output <- output
  makeQuartoReport(
    rawData = rawData,
    qmdFile = "ExploreMageckCounts.qmd",
    reportTitle = "ExploreMageckCounts",
    colour = isTRUE(param$colour),
    number = isTRUE(param$number)
  )
  return("Success")
}

##' @title Load a MAGeCK sgRNA count dataset into a SummarizedExperiment
##' @description Merges the per-sample MAGeCK \code{.count.txt} files (columns
##'   \code{sgRNA}, \code{Gene}, <count>) into one sgRNA x sample count matrix.
##'   The count column is generically named by mageck (e.g. \code{sample1}), so it
##'   is renamed to the SUSHI sample name per file before merging. Control sgRNAs
##'   are flagged from the library \code{*_MAGeCK_Ctrl.csv} and from control gene
##'   labels. No genome annotation is used -- these are guide counts, not genes.
loadMageckCountDataset <- function(input, param) {
  require(SummarizedExperiment)
  require(data.table)

  files <- input$getFullPaths("Count")
  if (length(files) == 0L) {
    return(list(error = "No 'Count' files found in the input dataset."))
  }
  sampleNames <- names(files)
  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    return(list(error = paste0(
      "Count file(s) not found:<br>",
      paste(missing, collapse = "<br>")
    )))
  }

  ## Read each sample, keep sgRNA + Gene + the single count column, and rename
  ## that count column to the sample name (mageck names it generically).
  perSample <- lapply(seq_along(files), function(i) {
    dt <- data.table::fread(files[[i]])
    countCol <- setdiff(colnames(dt), c("sgRNA", "Gene"))
    if (!all(c("sgRNA", "Gene") %in% colnames(dt)) || length(countCol) < 1L) {
      stop(
        "unexpected MAGeCK count format in ",
        files[[i]],
        " (need sgRNA, Gene, <count>)"
      )
    }
    out <- dt[, c("sgRNA", "Gene", countCol[1]), with = FALSE]
    data.table::setnames(out, countCol[1], sampleNames[i])
    out
  })

  merged <- perSample[[1]]
  for (i in seq_along(perSample)[-1]) {
    merged <- merge(
      merged,
      perSample[[i]][, c("sgRNA", sampleNames[i]), with = FALSE],
      by = "sgRNA",
      all = TRUE
    )
  }
  merged <- as.data.frame(merged)

  counts <- as.matrix(merged[, sampleNames, drop = FALSE])
  rownames(counts) <- merged[["sgRNA"]]
  counts[is.na(counts)] <- 0
  storage.mode(counts) <- "double"
  gene <- merged[["Gene"]]

  ## Control sgRNAs: from the library control-guide file, plus common control
  ## gene labels (Non-Targeting Control, Control, Safe Harbor).
  isControl <- rep(FALSE, nrow(counts))
  libName <- unique(input$getColumn("libName"))
  libName <- libName[nzchar(libName)][1]
  if (!is.na(libName) && dir.exists(libName)) {
    ctrlFile <- list.files(
      libName,
      pattern = "_MAGeCK_Ctrl\\.csv$",
      full.names = TRUE
    )
    if (length(ctrlFile) == 1L) {
      ctrlIds <- readLines(ctrlFile)
      ctrlIds <- trimws(ctrlIds[nzchar(ctrlIds)])
      isControl <- rownames(counts) %in% ctrlIds
    }
  }
  isControl <- isControl |
    grepl("non[- ]?targeting|^control$|safe.?harbor", gene, ignore.case = TRUE)

  ## Sample metadata, ordered to match the count columns.
  dataset <- input$meta
  dataset <- dataset[sampleNames, , drop = FALSE]

  design <- tryCatch(ezDesignFromDataset(dataset, param), error = function(e) NULL)
  if (is.null(design) || ncol(design) == 0L) {
    conds <- rep("all", ncol(counts))
    design <- data.frame(row.names = sampleNames)
  } else {
    conds <- ezConditionsFromDesign(design, maxFactors = 2)
  }

  rawData <- SummarizedExperiment(
    assays = SimpleList(counts = counts, presentFlag = counts > 0),
    rowData = DataFrame(Gene = gene, isControl = isControl),
    colData = DataFrame(dataset, conds = conds),
    metadata = list(
      param = param,
      design = design,
      featureLevel = "sgRNA",
      type = "Counts",
      countName = "count"
    )
  )
  return(rawData)
}

##' @template app-template
##' @templateVar method ezMethodMageckCountQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Exploratory sgRNA/gene count QC for pooled CRISPR screens
##'   (library representation, sample clustering, essential-gene separation) from
##'   MAGeCK count output.
EzAppExploreMageckCounts <-
  setRefClass(
    "EzAppExploreMageckCounts",
    contains = "EzApp",
    methods = list(
      ## MAGeCK count QC. Normalisation via DESeq2 size factors / edgeR (the same
      ## median-of-ratios idea MAGeCK uses); essential/non-essential separation
      ## uses the bundled Hart CEGv2 / NEGv1 gene sets.
      citation = function() {
        c(
          "Li, W. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology 15, 554 (2014). https://doi.org/10.1186/s13059-014-0554-4",
          "Rehrauer, H. et al. ezRun: An R meta-package for the analysis of Next Generation Sequencing Data. https://github.com/uzh/ezRun",
          "Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8",
          "Robinson, M.D., McCarthy, D.J. & Smyth, G.K. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26(1), 139-140 (2010). https://doi.org/10.1093/bioinformatics/btp616",
          "Hart, T. et al. Evaluation and Design of Genome-Wide CRISPR/SpCas9 Knockout Screens. G3 7(8), 2719-2727 (2017). https://doi.org/10.1534/g3.117.041277"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodMageckCountQC
        name <<- "EzAppExploreMageckCounts"
        appDefaults <<- rbind(
          normMethod = ezFrame(
            Type = "character",
            DefaultValue = "deseq2",
            Description = "count normalisation: deseq2 (size factors), tmm, cpm, or logMean"
          ),
          refGroup = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "reference/plasmid/T0 condition used as baseline for essential-gene depletion + ROC; empty = skip that analysis"
          ),
          backgroundExpression = ezFrame(
            Type = "numeric",
            DefaultValue = 5,
            Description = "pseudo-count added before log2 transform"
          ),
          topGeneSize = ezFrame(
            Type = "numeric",
            DefaultValue = 100,
            Description = "number of most-variable sgRNAs used for top-feature plots"
          ),
          nSampleClusters = ezFrame(
            Type = "numeric",
            DefaultValue = 6,
            Description = "number of sample clusters (cutree) in the high-variance heatmap"
          )
        )
      }
    )
  )
