% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{ezDuplicated}
\alias{ezDuplicated}
\title{Which values are duplicated?}
\usage{
ezDuplicated(x, mode = "keepFirst")
}
\arguments{
\item{x}{a vector, matrix or list.}

\item{mode}{specifies how to keep duplicated values. Possible modes:
\itemize{
 \item{"keepFirst"}{ uses \code{duplicated()} to keep duplicates in the original order.}
 \item{"keepLast"}{ keeps duplicates in a reversed order.}
 \item{"random"}{ keeps duplicates in a randomized order.}
 \item{"all"}{ keeps all duplicated values independent of the order.}
}}
}
\value{
Returns a logical indicating which values are duplicated according to the specified \code{mode}.
}
\description{
Keeps duplicated values by setting them to TRUE.
}
\examples{
v1 = c(1,2,3,4,5,4,3,2,3,4,5,6,7,8,7)
ezDuplicated(v1)
ezDuplicated(v1,"all")
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
