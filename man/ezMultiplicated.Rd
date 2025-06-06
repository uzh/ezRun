% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{ezMultiplicated}
\alias{ezMultiplicated}
\title{Which values occur at least \code{n} times?}
\usage{
ezMultiplicated(x, n = 2, mode = "keepFirst")
}
\arguments{
\item{x}{a vector, matrix or list.}

\item{n}{a positive integer specifying how many times a value needs to occur to return true.}

\item{mode}{specifies how to keep multiplicated values. Possible modes:
\itemize{
 \item{"keepFirst"}{ keeps values that occur \code{n} times in the original order.}
 \item{"keepLast"}{ keeps values that occur \code{n} times in a reversed order.}
 \item{"random"}{ keeps values that occur \code{n} times in a randomized order.}
 \item{"all"}{ keeps all values that occur \code{n} times independent of the order.}
}}
}
\value{
Returns a logical indicating which values are multiplicated according to the specified \code{mode} and \code{n}.
}
\description{
Keeps multiplicated values by setting them to TRUE.
}
\examples{
v1 = c(1,2,3,4,5,4,3,2,3,4,5,6,7,8,7)
ezMultiplicated(v1)
ezMultiplicated(v1,3)
ezMultiplicated(v1,2,"all")
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
