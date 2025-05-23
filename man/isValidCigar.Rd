% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ngsTools.R
\name{isValidCigar}
\alias{isValidCigar}
\title{Is \code{x} a valid cigar?}
\usage{
isValidCigar(x)
}
\arguments{
\item{x}{a character vector to check.}
}
\value{
Returns a logical vector answering which elements are valid cigars.
}
\description{
Tests whether \code{x} is a valid cigar by checking the format of the character.
}
\examples{
isValidCigar("3M5G")
}
\seealso{
\code{\link[GenomicAlignments]{cigarWidthAlongReferenceSpace}}
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
