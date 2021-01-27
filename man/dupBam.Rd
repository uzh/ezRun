\name{dupBam}
\alias{dupBam}
\title{
  Duplicates operation in bam
}
\description{
  Mark or remove the duplicates in bam file.
}
\usage{
  dupBam(inBam, outBam, operation=c("mark", "remove"))
}
\arguments{
  \item{inBam}{
    \code{character}(1): the input bam file.
  }
  \item{outBam}{
    \code{character}(1): the output bam file.
  }
  \item{operation}{
    \code{character}(1): either to mark or remove the duplicates in bam file.
  }
}
\details{
  It use \dQuote{sambamba} to mark and remove the duplicates,
  as the same criteria of \dQuote{picard}.

  Other alternatives are \dQuote{samtools} with \sQuote{fixmate} and
  \sQuote{markdup}, \dQuote{picard}'s \sQuote{MarkDuplicates}.
}
\value{
  invisible \code{outBam} file.
}
\author{
  Ge Tan
}
\examples{
  bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                         mustWork=TRUE)
  outBam <- tempfile(pattern="test-", fileext=".bam")
  dupBam(inBam=bamFile, outBam=outBam, operation="mark")
  dupBam(inBam=bamFile, outBam=outBam, operation="remove")
}
