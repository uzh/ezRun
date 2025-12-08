###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodDiffAnalysisPeak <- function(input = NA, output = NA, param = NA){
    grouping <- input$getColumn(param$grouping)
    stopifnot(param$sampleGroup != param$refGroup)

    countFiles <- input$getFullPathsList("Count")
    featureCounts <- loadCountFiles(countFiles, param, grouping)
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
loadCountFiles <- function(countFiles, param, grouping){
  countFilesSubset <- countFiles[names(grouping)[grouping %in% c(param$refGroup, param$sampleGroup)]]

  loadAllTables <- imap(countFilesSubset, function(file_i, listName) {
    group <- grouping[[listName]]
    rep_id <- match(listName, names(grouping)[grouping == group])
    new_col <- paste0(group, "_REP", rep_id)
    data.table::fread(file_i, data.table=FALSE) %>%
      rename(!!new_col := matchCounts)
  })

  reduce(loadAllTables, full_join, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
}

