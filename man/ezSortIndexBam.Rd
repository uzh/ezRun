\name{ezSortIndexBam}
\alias{ezSortIndexBam}
\title{Sorts and indexes a bam file}
\usage{
ezSortIndexBam(inBam, bam, ram=2, removeBam=TRUE, cores=2,
               method=c("sambamba", "samtools"))
}
\arguments{
\item{inBam}{a character representing the file path of the input bam file to sort.}

\item{bam}{a character representing the file path of the output bam file.}

\item{ram}{max RAM to use in GB}

\item{removeBam}{a logical indicating whether to remove the file specified in \code{inBam}.}

\item{cores}{an integer specifying how many cores to use for the command.}
  \item{method}{
    \code{character}(1): either using \dQuote{sambamba} for better performance
    or \dQuote{samtools}.
  }
}
\description{
  Sorts and indexes a bam file.
}
\examples{
  bamFile <- system.file("extdata", "ex1.bam", package="Rsamtools", 
                         mustWork=TRUE)
  ezSortIndexBam(inBam=bamFile, bam="test.bam", ram=0.5, 
                 removeBam=FALSE, cores = ezThreads())
}
\author{
Rehrauer, Hubert
Schmid, Peter
}
