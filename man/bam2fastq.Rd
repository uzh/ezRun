\name{bam2fastq}
\alias{bam2fastq}
\title{
  Bam/Sam file to fastq file conversion
}
\description{
  Convert either paired or unpaired Bam file into fastq files.
}
\usage{
  bam2fastq(bamFns, fastqFns = sub("(\\.bam|\\.sam)$", "fastq", bamFns),
            fastq2Fns=sub("(\\.bam|\\.sam)$", "_R2.fastq", bamFns),
            paired=FALSE)
}
\arguments{
  \item{bamFns}{
    \code{character}(n): paths of input Bam/Sam files
  }
  \item{fastqFns}{
    \code{character}(n): paths of first end of read files.
  }
  \item{fastq2Fns}{
    \code{character}(n): optional paths of second end of paired read files.
  }
  \item{paired}{
    \code{boolean}(1): paired-end or single-end read file.
  }
}
\details{
  The conversion is done with picard's SamToFastq.
  
  Rsamtools::testPairedEndBam cannot be used to test for paired or not.
  It doesn't work on Bam converted from fastq, due to the lack of header.
  \code{paired} has to be specified explicitly.
}
\value{
  invisible \code{fastqFns}.
}
\author{
  Ge Tan
}
\seealso{
  \code{\link{fastqs2bam}}
}
\examples{
  \dontrun{
    fastqFns <- list.files(path="/srv/gstore/projects/p2288/HiSeq2500_20171011_RR99_o3511/dmx",
                           pattern="\\\.fastq\\\.gz$", full.names=TRUE)
    bamFn <- fastqs2bam(fastqFns[1:5], bamFn="20171011.A-C1_HT_24H.bam")
    fastqFn <- bam2fastq(bamFn, paired=FALSE)
  }
}
