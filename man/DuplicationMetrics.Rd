\name{DuplicationMetrics}
\alias{DuplicationMetrics}
\title{
  DuplicationMetrics
}
\description{
  Collection DuplicationMetrics with picard.
}
\usage{
  DuplicationMetrics(inBams, mc.cores = ezThreads())
}
\arguments{
  \item{inBams}{
    \code{character}(n): the filepaths of input bam files.
  }
  \item{mc.cores}{
    \code{integer}(1): the threads to use.
  }
}
\value{
   A \code{data.frame} of metrics.
}
\references{
  \url{http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics}
}
\author{
  Ge Tan
}
\seealso{
  \code{\link{CollectRnaSeqMetrics}}
}
\examples{
  \dontrun{
    inBam <- "/srv/gstore/projects/p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5.bam"
    DuplicationMetrics(inBam)
  }
}
