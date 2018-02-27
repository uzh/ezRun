\name{CollectRnaSeqMetrics}
\alias{CollectRnaSeqMetrics}
\title{
  CollectRnaSeqMetrics
}
\description{
  Produces a summary of RNA-seq metrics from a SAM or BAM file with picard.
}
\usage{
  CollectRnaSeqMetrics(inBams, gtfFn, featAnnoFn,
                       strandMode = c("both", "sense", "antisense"),
                       metricLevel = c("ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"),
                       mc.cores = ezThreads())
}
\arguments{
  \item{inBams}{
    \code{character}(n): the filepaths of input bam files.
  }
  \item{gtfFn}{
    \code{character}(1): the filepath of gtf file.
  }
  \item{featAnnoFn}{
    \code{character}(1): the filepath of feature annotation file.
  }
  \item{strandMode}{
    \code{character}(1):  \dQuote{both}, \dQuote{sense} or \dQuote{antisense}.
  }
  \item{metricLevel}{
    \code{character}(1): the level of summarisation.
  }
  \item{mc.cores}{
    \code{integer}(1): the threads to use.
  }
}
\value{
  A \code{data.frame} of metrics.
}
\references{
  \url{http://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics}
}
\author{
  Ge Tan
}
\seealso{
  \code{\link{CollectAlignmentSummaryMetrics}}
}
\examples{
  \dontrun{
    inBam <- "/srv/gstore/projects/p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5.bam"
    gtfFn <- "/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31/Genes/genes.gtf"
    featAnnoFn <- "/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31/Genes/genes_annotation_byTranscript.txt"
    metricLevel <- "SAMPLE"
    CollectRnaSeqMetrics(inBam, gtfFn, featAnnoFn, strandMode="both", 
                         metricLevel)
  }
}
