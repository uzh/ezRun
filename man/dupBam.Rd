\name{dupBam}
\alias{dupBam}
\title{
  Duplicates operation in bam
}
\description{
  Mark or remove the duplicates in bam file.
}
\usage{
  dupBam(inBam, outBam, operation=c("mark", "remove"),
         program=c("sambamba", "picard"), cores=ezThreads())
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
  \item{program}{
    \code{character}(1): either use \dQuote{sambamba} or \dQuote{picard} to mark the duplicates.
  }
  \item{cores}{
    \code{integer}(1): the number of threads to use.
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
\references{
  \url{http://lomereiter.github.io/sambamba/docs/sambamba-markdup.html}
}
\author{
  Ge Tan
}
\examples{
  bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", 
                         mustWork=TRUE)
  outBam <- tempfile(pattern="test-", fileext=".bam")
  dupBam(inBam=bamFile, outBam=outBam, program="sambamba", operation="mark")
  dupBam(inBam=bamFile, outBam=outBam, program="sambamba", operation="remove")
  dupBam(inBam=bamFile, outBam=outBam, program="picard", operation="mark")
  dupBam(inBam=bamFile, outBam=outBam, program="picard", operation="remove")
}
