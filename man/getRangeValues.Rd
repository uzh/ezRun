% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rangesTools.R
\name{getRangeValues}
\alias{getRangeValues}
\title{Gets the values of genomic ranges}
\usage{
getRangeValues(
  cov,
  targetRanges,
  doRev = as.character(strand(targetRanges)) == "-",
  asMatrix = TRUE
)
}
\arguments{
\item{cov}{an object of the class RleList.}

\item{targetRanges}{an object of the class GRanges.}

\item{doRev}{a logical vector indicating whether to reverse ranges. Usually derived from \code{targetRanges}.}

\item{asMatrix}{a logical indicating whether to return the values as a matrix.}
}
\value{
Returns the values or a matrix containing them.
}
\description{
Gets the values of genomic ranges
}
\seealso{
\code{\link[IRanges]{RleViewsList}}

\code{\link[IRanges]{viewApply}}
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
