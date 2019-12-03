\name{bam2bw}
\alias{bam2bw}
\title{
  Bam to bigwig
}
\description{
  Create bigwig file from bam file.
}
\usage{
  bam2bw(file, destination = sub("\\.bam$", ".bw", ignore.case = TRUE), 
         paired = FALSE, method = c("deepTools", "Bioconductor"),
         cores = ezThreads())
}
\arguments{
  \item{file}{
    \code{character}(1): the filepath of input bam file.
  }
  \item{destination}{
    \code{character}(1): the filepath of output bigwig file.
  }
  \item{paired}{
    \code{boolean}(1): the bam file is from paired-end or single-end reads.
  }
  \item{method}{
    \code{character}(1): use the \dQuote{deepTools} or \dQuote{Bioconductor}.
  }
  \item{cores}{
    \code{integer}(1): the number of threads to use.
  }
}
\details{
  The \dQuote{bamCoverage} tool from \dQuote{deepTools} is faster,
  more memory efficient. By default, we export a bigwig file with bin size of 10.
  
  The \dQuote{Bioconductor} approach is more self-contained, slower and consuming
  much more memory. It output single-base resolution.
}
\value{
  An invisible \code{destination}.
}
\references{
  \url{http://deeptools.readthedocs.io/en/latest/index.html}
}
\author{
  Ge Tan
}
\note{
  
}
\seealso{
  \code{\link{ezBam2bigwig}}
}
\examples{
  file <- system.file("extdata", "ex1.bam", package="Rsamtools", 
                      mustWork=TRUE)
  destination <- tempfile(pattern="test-", fileext=".bw")
  bam2bw(file, destination, paired=FALSE, method="Bioconductor")
}
