\name{fastqs2bam}
\alias{fastqs2bam}
\title{
  Fastq files to Bam with RG tag
}
\description{
  Convert multiple fastq files into Bam files and merged into one big Bam file
  with filenames as RG tag.
}
\usage{
  fastqs2bam(fastqFns, fastq2Fns=NULL, readGroupNames=NULL, bamFn)
}
\arguments{
  \item{fastqFns}{
    \code{character}(n): paths of input fastq files
  }
  \item{fastq2Fns}{
    \code{character}(n): optional paths of the second read of paired-end data.
  }
  \item{readGroupNames}{
    \code{character}(n): read group IDs. Same length as \code{fastqFns}.
  }
  \item{bamFn}{
    \code{character}(1): filename of merged Bam file.
  }
}
\details{
  Each fastq file first is converted into Bam file with filename
  as \dQuote{SM} and \code{readGroupNames} as \sQuote{ID} in RG tag 
  using picard's FastqToSam.
  Then Bam files are merged into one big Bam file using picard's MergeSamFiles
  with a sort order by queryname.
}
\value{
  invisible \code{bamFn}.
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{bam2fastq}}
}
\examples{
  \dontrun{
    fastqFns <- list.files(path="/srv/gstore/projects/p2288/HiSeq2500_20171011_RR99_o3511/dmx",
                           pattern="\\\.fastq\\\.gz$", full.names=TRUE)
    fastqFns <- fastqFns[1:5]
    fastqs2bam(fastqFns,
               readGroupNames=sub("\\\.(fastq|fq)(\\\.gz){0,1}$", "", 
                                  basename(fastqFns)),
               bamFn="20171011.A-C1_HT_24H_unmapped.bam")
  }
}
