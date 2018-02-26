###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

CollectAlignmentSummaryMetrics <- function(inBam, fastaFn,
                                           metricLevel=c("ALL_READS", "SAMPLE",
                                                         "LIBRARY", "READ_GROUP")){
  require(matrixStats)
  metricLevel <- match.arg(metricLevel)
  setEnvironments("picard")
  outputFn <- tempfile(pattern="CollectAlignmentSummaryMetrics",
                       fileext=".txt")
  cmd <- paste("java -Xmx50g -jar", Sys.getenv("Picard_jar"),
               "CollectAlignmentSummaryMetrics",
               paste0("R=", fastaFn),
               paste0("I=", inBam),
               paste0("O=", outputFn),
               "METRIC_ACCUMULATION_LEVEL=null", ## clear the default "ALL_READS"
               paste0("METRIC_ACCUMULATION_LEVEL=", metricLevel))
  ezSystem(cmd)
  
  metrics <- ezRead.table(outputFn, comment.char="#", row.names=NULL)
  metrics <- metrics[ ,!colAlls(is.na(metrics))] ## Remove the NA columns of nameColumns
  nameColumns <- c("SAMPLE", "LIBRARY", "READ_GROUP")
  indexName <- intersect(colnames(metrics), nameColumns)
  rownames(metrics) <- metrics[ ,indexName]
  metrics[[indexName]] <- NULL
  return(metrics)
}

CollectRnaSeqMetrics <- function(inBam, gtfFn, featAnnoFn,
                                 strandMode=c("both", "sense", "antisense"),
                                 metricLevel=c("ALL_READS", "SAMPLE",
                                               "LIBRARY", "READ_GROUP")
                                 ){
  setEnvironments("UCSC")
  setEnvironments("samtools")
  setEnvironments("picard")
  
  strandMode <- match.arg(strandMode)
  metricLevel <- match.arg(metricLevel)
  
  strandMode <- switch(strandMode,
                       "both"="NONE",
                       "sense"="FIRST_READ_TRANSCRIPTION_STRAND",
                       "antisense"="SECOND_READ_TRANSCRIPTION_STRAND")
  ## GTF to refFlat
  refFlatFn <- tempfile(pattern="refFlat", fileext=".txt")
  cmd <- paste("gtfToGenePred -genePredExt -geneNameAsName2",
               "-ignoreGroupsWithoutExons", gtfFn, refFlatFn)
  ezSystem(cmd)
  refFlat <- read.table(refFlatFn, sep="\t", header=FALSE,
                        stringsAsFactors = FALSE)
  refFlat <- refFlat[ , c(12, 1:10)]
  write.table(refFlat, file=refFlatFn, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)
  
  ## RIBOSOMAL_INTERVALS
  riboFn <- tempfile(pattern="ribosomal", fileext=".interval")
  cmd <- paste("samtools view -H", inBam, ">", riboFn)
  ezSystem(cmd)
  
  gtf <- ezReadGff(gtfFn)
  featAnno <- ezFeatureAnnotation(featAnnoFn, dataFeatureType="transcript")
  featAnno <- featAnno[featAnno$type == "rRNA", ]
  featAnno <- featAnno[ , c("seqid", "start", "end", "strand", "transcript_id")]
  write.table(featAnno, file=riboFn, quote=FALSE, sep="\t",
              append=TRUE, row.names=FALSE, col.names = FALSE)
  
  ## CollectRnaSeqMetrics
  outputFn <- tempfile(pattern="CollectRnaSeqMetrics",
                       fileext=".txt")
  cmd <- paste("java -Xmx50g -jar", Sys.getenv("Picard_jar"), ## it needs big RAM
               "CollectRnaSeqMetrics",
               paste0("I=", inBam),
               paste0("O=", outputFn),
               paste0("REF_FLAT=", refFlatFn),
               paste0("STRAND=", strandMode),
               paste0("RIBOSOMAL_INTERVALS=", riboFn),
               "METRIC_ACCUMULATION_LEVEL=null", ## clear the default "ALL_READS"
               paste0("METRIC_ACCUMULATION_LEVEL=", metricLevel)
               )
  ezSystem(cmd)
  metricsAll <- readLines(outputFn)
  outputFn1 <- tempfile(pattern="CollectRnaSeqMetrics_1_",
                        fileext=".txt")
  breakLine <- grep("^## HISTOGRAM", metricsAll)
  writeLines(head(metricsAll, breakLine-2), con=outputFn1)
  
  metrics <- ezRead.table(outputFn1, comment.char="#", row.names=NULL)
  metrics <- metrics[ ,!colAlls(is.na(metrics))]
  nameColumns <- c("SAMPLE", "LIBRARY", "READ_GROUP")
  indexName <- intersect(colnames(metrics), nameColumns)
  rownames(metrics) <- metrics[ ,indexName]
  metrics[[indexName]] <- NULL
  
  outputFn2 <- tempfile(pattern="CollectRnaSeqMetrics_2_",
                        fileext=".txt")
  writeLines(metricsAll[breakLine:(102+breakLine)], con=outputFn2)
  histogram <- ezRead.table(outputFn2, comment.char="#", row.names=NULL)
  
  ans <- list(metrics, histogram)
  return(ans)
}
