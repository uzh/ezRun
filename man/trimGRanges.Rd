\name{trimGRanges}
\alias{trimGRanges}
\title{
  trim GRanges object
}
\description{
  trim to get the first or last bases from GRanges.
}
\usage{
  trimGRanges(x, width = 100, start = TRUE)
}
\arguments{
  \item{x}{
    A \code{GRanges} object.
  }
  \item{width}{
    \code{integer}(1): the length of bases to fetch.
  }
  \item{start}{
    \code{logical}(1): fetch from the start or end.
  }
}
\details{
  Image we have exons from the same transcript in the \code{x}.
  We want to fetch 100 bases of CDS from 5' end or 3' end.
  
  
}
\value{
  A \code{GRanges} object.
}
\author{
  Ge Tan
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  gr <- GRanges(seqnames="chr1",
                ranges=IRanges(start=c(2,10,20),
                               end=c(5, 15, 23)),
                strand="+")
  trimGRanges(gr, width=3, start=TRUE)
  trimGRanges(gr, width=4, start=TRUE)
  trimGRanges(gr, width=11, start=TRUE)
  trimGRanges(gr, width=11, start=FALSE)
}
