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
  fastqs2bam(fastqFns, bamFn)
}
\arguments{
  \item{fastqFns}{
    \code{character}(n): paths of input fastq files
  }
  \item{bamFn}{
    \code{character}(1): filename of merged Bam file.
  }
}
\details{
  Each fastq file first is converted into Bam file with filename as RG tag
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  \dontrun{
    fastqFns <- list.files(path="/srv/gstore/projects/p2288/HiSeq2500_20171011_RR99_o3511/dmx",
                           pattern="\\.fastq\\.gz$", full.names=TRUE)
    fastqs2bam(fastqFns, bamFn="20171011.A-C1_HT_24H.bam")
  }
}
