###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

CollectAlignmentSummaryMetrics <- function(inBams, fastaFn, paired=FALSE,
                                           metricLevel=c("ALL_READS", "SAMPLE",
                                                         "LIBRARY", "READ_GROUP"),
                                           mc.cores=ezThreads()){
  require(matrixStats)
  metricLevel <- match.arg(metricLevel)
  setEnvironments("picard")
  
  outputFns <- tempfile(pattern=inBams,
                       fileext=".CollectAlignmentSummaryMetrics")
  ## TODO: do check the ADAPTER_SEQUENCES option
  cmd <- paste("java -Xmx3G", 
               "-jar", Sys.getenv("Picard_jar"),
               "CollectAlignmentSummaryMetrics",
               paste0("R=", fastaFn),
               paste0("I=", inBams),
               paste0("O=", outputFns),
               "METRIC_ACCUMULATION_LEVEL=null", ## clear the default "ALL_READS"
               paste0("METRIC_ACCUMULATION_LEVEL=", metricLevel),
               "> /dev/null")
  ezMclapply(cmd, ezSystem, mc.preschedule=FALSE, mc.cores=mc.cores)
  
  metrics <- lapply(outputFns, ezRead.table, comment.char="#", row.names=NULL)
  metrics <- do.call(rbind, metrics)
  metrics <- metrics[ ,!colAlls(is.na(metrics))] ## Remove the NA columns of nameColumns
  if (paired){
    metrics = metrics[ metrics$CATEGORY == "PAIR", ]
  }
  rownames(metrics) = names(inBams)
  # nameColumns <- c("SAMPLE", "LIBRARY", "READ_GROUP")
  # indexName <- intersect(colnames(metrics), nameColumns)
  # rownames(metrics) <- metrics[ ,indexName]
  # metrics[[indexName]] <- NULL
  return(metrics)
}

CollectRnaSeqMetrics <- function(inBams, gtfFn, featAnnoFn,
                                 strandMode=c("both", "sense", "antisense"),
                                 metricLevel=c("ALL_READS", "SAMPLE",
                                               "LIBRARY", "READ_GROUP"),
                                 mc.cores=ezThreads()
                                 ){
  require(matrixStats)
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
  cmd <- paste("samtools view -H", inBams[1], ">", riboFn)
  ezSystem(cmd)
  
  gtf <- ezReadGff(gtfFn)
  featAnno <- ezFeatureAnnotation(featAnnoFn, dataFeatureType="transcript")
  featAnno <- featAnno[featAnno$type == "rRNA", ]
  featAnno <- featAnno[ , c("seqid", "start", "end", "strand", "transcript_id")]
  write.table(featAnno, file=riboFn, quote=FALSE, sep="\t",
              append=TRUE, row.names=FALSE, col.names = FALSE)
  
  ## CollectRnaSeqMetrics
  outputFns <- tempfile(pattern=inBams,
                        fileext=".CollectRnaSeqMetrics")
  cmd <- paste("java -Xmx3G", "-jar", Sys.getenv("Picard_jar"),
               "CollectRnaSeqMetrics",
               paste0("I=", inBams),
               paste0("O=", outputFns),
               paste0("REF_FLAT=", refFlatFn),
               paste0("STRAND=", strandMode),
               paste0("RIBOSOMAL_INTERVALS=", riboFn),
               "METRIC_ACCUMULATION_LEVEL=null", ## clear the default "ALL_READS"
               paste0("METRIC_ACCUMULATION_LEVEL=", metricLevel),
               "> /dev/null"
               )
  ezMclapply(cmd, ezSystem, mc.preschedule=FALSE, mc.cores=mc.cores/2)
  ## To save memory
  
  for(outputFn in outputFns){
    ## Remove the stuff after HISTOGRAM
    metricsAll <- readLines(outputFn)
    breakLine <- grep("^## HISTOGRAM", metricsAll)
    if(length(breakLine) >= 1){
      writeLines(head(metricsAll, breakLine-2), con=outputFn)
    }
  }
  metrics <- lapply(outputFns, ezRead.table, comment.char="#", row.names=NULL)
  metrics <- do.call(rbind, metrics)
  metrics <- metrics[ ,!colAlls(is.na(metrics))] ## Remove the NA columns of nameColumns
  rownames(metrics) = names(inBams)
  # nameColumns <- c("SAMPLE", "LIBRARY", "READ_GROUP")
  # indexName <- intersect(colnames(metrics), nameColumns)
  # rownames(metrics) <- metrics[ ,indexName]
  # metrics[[indexName]] <- NULL
  
  # 
  # writeLines(head(metricsAll, breakLine-2), con=outputFn1)
  # 
  # metrics <- ezRead.table(outputFn1, comment.char="#", row.names=NULL)
  # metrics <- metrics[ ,!colAlls(is.na(metrics))]
  # nameColumns <- c("SAMPLE", "LIBRARY", "READ_GROUP")
  # indexName <- intersect(colnames(metrics), nameColumns)
  # rownames(metrics) <- metrics[ ,indexName]
  # metrics[[indexName]] <- NULL
  # 
  # outputFn2 <- tempfile(pattern="CollectRnaSeqMetrics_2_",
  #                       fileext=".txt")
  # writeLines(metricsAll[breakLine:(102+breakLine)], con=outputFn2)
  # histogram <- ezRead.table(outputFn2, comment.char="#", row.names=NULL)
  # 
  # ans <- list(metrics, histogram)
  return(metrics)
}

DuplicationMetrics <- function(inBams, mc.cores=ezThreads()){
  setEnvironments("picard")
  
  outputBams <- paste0(sub(".bam$", "", basename(inBams)),
                       "-", seq_along(inBams), "-markedDup.bam")
  outputFns <- sub('bam$', "metrics", outputBams)
  
  cmd <- paste("java -XX:ParallelGCThreads=4 -Xmx3G -jar",
               Sys.getenv("Picard_jar"), "MarkDuplicates",
               paste0("I=", inBams),
               paste0("O=", outputBams),
               paste0("M=", outputFns),
               "> /dev/null"
               )
  ezMclapply(cmd, ezSystem, mc.preschedule=FALSE, mc.cores=mc.cores)
  
  metrics <- lapply(outputFns, ezRead.table, comment.char="#", row.names=NULL, blank.lines.skip=TRUE, nrows=1)
  metrics <- do.call(rbind, metrics)
  rownames(metrics) = names(inBams)
  # metrics <- metrics[ ,!colAlls(is.na(metrics))] ## Remove the NA columns of nameColumns
  # nameColumns <- c("SAMPLE", "LIBRARY", "READ_GROUP")
  # indexName <- intersect(colnames(metrics), nameColumns)
  # rownames(metrics) <- metrics[ ,indexName]
  # metrics[[indexName]] <- NULL
  file.remove(outputBams, outputFns) ## only remove on success so theat se an debg
  return(metrics)
}
