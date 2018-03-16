\name{bam2fastq}
\alias{bam2fastq}
\title{
  Bam/Sam file to fastq file conversion
}
\description{
  Convert either paired or unpaired Bam file into fastq files.
}
\usage{
  bam2fastq(bamFn, OUTPUT_PER_RG=TRUE, OUTPUT_DIR=".", paired=FALSE,
            fastqFns=sub("(\\.bam|\\.sam)$", "fastq", bamFn),
            fastq2Fns=sub("(\\.bam|\\.sam)$", "_R2.fastq", bamFn)
            )
}
\arguments{
  \item{bamFn}{
    \code{character}(1): path of input Bam/Sam files
  }
  \item{OUTPUT_PER_RG}{
    \code{boolean}(1): whether or not to output one fastq per RG.
    When TRUE, \code{paired}, \code{fastqFns} and \code{fastq2Fns} are ignoreed.
  }
  \item{OUTPUT_DIR}{
    \code{boolean}(1): the output dir for fastqs.
  }
  \item{fastqFns}{
    \code{character}(1): paths of first end of read file.
  }
  \item{fastq2Fns}{
    \code{character}(2): optional paths of second end of paired read file.
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
  invisible filenames of first read files.
}
\author{
  Ge Tan
}
\seealso{
  \code{\link{fastqs2bam}}
}
\examples{
  \dontrun{
    bamFn <- "/srv/gstore/projects/p2497/HiSeq4000_20171222_RUN420_o3705_fixedRG/20171222.A-SiCSeq_SCs_P5_unmapped.bam"
    bamFn <- "/scratch/gtan/p2497-SCCountQC/20171222.A-SiCSeq_SCs_P5_subset.bam"
    
    ## Split into one Fastq
    fastqFn <- bam2fastq(bamFn, fastqFns=sub("\\\.bam$", ".fastq", basename(bamFn)), 
                         paired=FALSE)
    ## Split by RG
    fastqFn <- bam2fastq(bamFn, OUTPUT_PER_RG=TRUE, OUTPUT_DIR=".", paired=FALSE)
  }
}
