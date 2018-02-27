\name{CollectAlignmentSummaryMetrics}
\alias{CollectAlignmentSummaryMetrics}
\title{
  CollectAlignmentSummaryMetrics
}
\description{
  Produces a summary of alignment metrics from a SAM or BAM file with picard.
}
\usage{
  CollectAlignmentSummaryMetrics(inBams, fastaFn,
                                 metricLevel = c("ALL_READS", "SAMPLE",
                                                 "LIBRARY", "READ_GROUP"),
                                 mc.cores=ezThreads())
}
\arguments{
  \item{inBams}{
    \code{character}(n): the filepaths of input bam files.
  }
  \item{fastaFn}{
    \code{character}(1): the filepath of reference genome file for mapping.
  }
  \item{metricLevel}{
    \code{character}(1): the level of summarisation. See details below.
  }
  \item{mc.cores}{
    \code{integer}(1): the threads to use.
  }
}
\details{
  This function relies on the picard tools's CollectAlignmentSummaryMetrics.
  By default, it computes the metric on all reads of the bam file.
  When read group tag is used, it's recommended to summarise over \dQuote{SAMPLE},
  \dQuote{LIBRARY} or \dQuote{READ_GROUP}.
}
\value{
  A \code{data.frame} of metrics.
}
\references{
  \url{http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics}
}
\author{
  Ge Tan
}
\examples{
  \dontrun{
    inBam <- "/srv/gstore/projects/p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5.bam"
    fastaFn <- "/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Sequence/WholeGenomeFasta/genome.fa"
    metricLevel <- "SAMPLE"
    CollectAlignmentSummaryMetrics(inBam, fastaFn, metricLevel)
  }
}
