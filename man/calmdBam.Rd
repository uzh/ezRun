\name{calmdBam}
\alias{calmdBam}
\title{
  Calculate MD tag for BAM
}
\description{
  Add the MD tag to BAM file and replace matches with \dQuote{=} in \sQuote{seq}.
}
\usage{
  calmdBam(bamFns, genomeFn, mc.cores=4L)
}
\arguments{
  \item{bamFns}{
    \code{character}(n): filenames of BAM files.
  }
  \item{genomeFn}{
    \code{character}(1): the reference fasta file for the alignment.
  }
  \item{mc.cores}{
    \code{integer}(1): the number of cores to use.
  }
}
\details{
  It uses the \dQuote{calmd} from samtools.
}
\value{
  Filenames of returned BAM files.
}

\examples{
  bamFns <- list.files("/scratch/gtan/R2ErrorRates/GEORNASeqBam",
                       pattern="\\\.bam$", recursive = TRUE, full.names = TRUE)
  calmdBam(bamFns, 
           genomeFn="/scratch/gtan/R2ErrorRates/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26/Genes/BOWTIE2Index/transcripts.fa")
}
