###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodDiffPeakAnalysis <- function(input = NA, output = NA, param = NA){
    grouping <- input$getColumn(param$grouping)
    stopifnot(param$sampleGroup != param$refGroup)
    grouping <- grouping[grouping %in% c(param$refGroup, param$sampleGroup)]
    commonCols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

    countFiles <- input$getFullPathsList("Count")
    featureCounts <- loadCountFiles(countFiles, grouping, commonCols)

    dds <- generateDESeqDS(featureCounts, commonCols, grouping)
    outDir <- file.path(basename(output$getColumn('ResultFolder')), 'diffpeak_analysis')
    cd = getwd()
    setwdNew(outDir)
    makeRmdReport(
      output = output, param = param, dds=dds, selfContained = TRUE,
      rmdFile = "DiffPeakAnalysis.Rmd", htmlFile = "DifferentialPeakAnalysisReport.html",
      reportTitle = 'Differential Peak Analysis', use.qs2 = TRUE
    )
    setwd(cd)
}


#' @template app-template
##' @templateVar method ezMethodDiffPeakAnalysis(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run a differential expression analysis with the application edgeR on two groups.
EzAppDiffPeakAnalysis <-
  setRefClass("EzAppDiffPeakAnalysis",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodDiffPeakAnalysis
        name <<- "EzAppDiffPeakAnalysis"
      }
      )
    )

##' @description generate file with all counts for sampleGroup and refGroup
loadCountFiles <- function(countFiles, grouping, commonCols){
  countFilesSubset <- countFiles[names(grouping)]

  loadAllTables <- imap(countFilesSubset, function(file_i, listName) {
    group <- grouping[[listName]]
    rep_id <- match(listName, names(grouping)[grouping == group])
    new_col <- paste0(group, "_REP", rep_id)
    data.table::fread(file_i, data.table=FALSE) %>%
      rename(!!new_col := matchCounts)
  })

  reduce(loadAllTables, full_join, by = commonCols)
}

##' @description generate DESeqDataSet from the counts table
generateDESeqDS <- function(featureCounts, commonCols, grouping){
  library(DESeq2)
  countCols <- setdiff(colnames(featureCounts), commonCols)

  countData <- featureCounts %>%
    select(all_of(countCols)) %>%
    as.data.frame()
  countData <- round(countData)

  colData <- ezFrame(
    sample = countCols,
    group  = sub("_REP.*", "", countCols),
    replicate = sub(".*_REP", "", countCols),
    row.names = countCols
  )

  dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData   = colData,
    design    = ~ group
  )
  rowData(dds) <- featureCounts[, commonCols]
  dds$Condition <- dds$group
  dds
}