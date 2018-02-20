###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

CollectAlignmentSummaryMetrics <- function(inBam, fastaFn,
                                           metricLevel=c("ALL_READS", "SAMPLE",
                                                         "LIBRARY", "READ_GROUP")){
  metricLevel <- match.arg(metricLevel)
  setEnvironments("picard")
  outputFn <- tempfile(pattern="CollectAlignmentSummaryMetrics",
                       fileext=".txt")
  cmd <- paste("java -jar", Sys.getenv("Picard_jar"),
               "CollectAlignmentSummaryMetrics",
               paste0("R=", fastaFn),
               paste0("I=", inBam),
               paste0("O=", outputFn),
               paste0("METRIC_ACCUMULATION_LEVEL=", metricLevel))
  ezSystem(cmd)
  
  metrics <- ezRead.table(outputFn, comment.char="#", row.names=NULL)
  nameColumns <- c("SAMPLE", "LIBRARY", "READ_GROUP")
  metrics <- metrics[ , c(nameColumns, setdiff(colnames(metrics), nameColumns))]
  return(metrics)
}
