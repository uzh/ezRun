\name{splitBamByRG}
\alias{splitBamByRG}
\title{
  Split bam file
}
\description{
  Split bam file wit read group.
}
\usage{
  splitBamByRG(inBam, mc.cores = ezThreads())
}
\arguments{
  \item{inBam}{
    \code{character}(1): filename of input sorted bam file.
  }
  \item{mc.cores}{
    \code{integer}(1): the threads to use.
  }
}
\details{
  It relies on samtools split.
}
\value{
  \code{character}(n): the filenames of splitted bam files.
}
\author{
  Ge Tan
}
\examples{
  \dontrun{
    inBam <- "/srv/gstore/projects/p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5.bam"
    inBam <- "/scratch/gtan/p2497-SCCountQC/20171222.A-SiCSeq_SCs_P5_subset.bam"
    splitBamByRG(inBam, mc.cores=8L)
  }
}
