\name{consensusPeaks}
\alias{consensusPeaks}
\title{
  Get consensus peaks
}
\description{
  Get the consensus peaks for ChIP-seq or ATAC-seq peaks.
}
\usage{
  consensusPeaks(x)
}
  
\arguments{
  \item{x}{
    \code{GRangesList}(1): a list of peaks in \code{GRanges}.
    The first \code{GRanges} has to contain the peaks from pooled reads.
    The other GRanges are assumed to be peaks called from each replicate.
  }
}
\details{
  To deal with the replicates and get the consensus peaks,
  the strategy is to pool the reads from replicates and call the peaks.
  
  This function produces a much more conservative set of peaks.
  For each peak from pooled reads, the peak is kept only when the overlap with
  the peak from each replicate is larger than half of either peak.
}
\value{
  \code{GRanges} object of final peaks.
}
\author{
  Ge Tan
}
\examples{
  \dontrun{
    ## consensusPeaks
    bedFns <- file.path("/srv/gstore/projects/p2578/MACS2_22192_2018-01-29--14-07-06",
                        c("A0062P_peaks.bed", "A0101P_peaks.bed", 
                          "A3401P_peaks.bed"))
    bedGRL <- GRangesList(lapply(bedFns, import.Macs2Peaks))
    consensusPeaks(bedGRL)
  }
}
