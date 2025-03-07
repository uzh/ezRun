% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{ezCollapse}
\alias{ezCollapse}
\title{Collapses a vector in a single character}
\usage{
ezCollapse(x, sep = "; ", na.rm = FALSE, empty.rm = FALSE, uniqueOnly = FALSE)
}
\arguments{
\item{x}{a vector, matrix or list.}

\item{sep}{the separator to use between values.}

\item{na.rm}{a logical specifying whether to remove \code{NA}'s.}

\item{empty.rm}{a logical specifying whether to remove empty values.}

\item{uniqueOnly}{a logical specifying whether to keep only unique values.}
}
\value{
Returns the values collapsed into one character.
}
\description{
This extends the functionality from \code{paste(..., collapse=...)} by optionally removing empty characters, duplicates or NA values
}
\examples{
l1 = list(a=c(1,"",6),c=c("rsrg","yjrt",NA,6))
ezCollapse(l1,sep="_")
ezCollapse(l1,na.rm=T,empty.rm=T,uniqueOnly=T)
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
