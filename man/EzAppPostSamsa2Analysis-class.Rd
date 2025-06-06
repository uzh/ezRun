% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/app-PostSamsa2Analysis.R
\docType{class}
\name{EzAppPostSamsa2Analysis-class}
\alias{EzAppPostSamsa2Analysis-class}
\alias{EzAppPostSamsa2Analysis}
\title{The R5 class representing a runnable app using \code{ezMethodPostSamsa2Analysis()()}}
\arguments{
\item{input}{a list, file path or an object of the class EzDataset containing the input.}

\item{output}{a list, file path or an object of the class EzDataset containing the output information.}

\item{param}{a list of parameters to customize the application run.}
}
\description{
Use this reference class to run
}
\section{Fields}{

\describe{
\item{\code{runMethod}}{the function that will be executed in the \code{run} method.}

\item{\code{name}}{the name of the application.}

\item{\code{appDefaults}}{the defaults to run the application with.}
}}

\section{Methods}{

\describe{
\item{\code{initialize()}}{Initializes the application using its specific defaults.}
}}

\section{Functions}{

\itemize{
  \item{\code{ezMethodPostSamsa2Analysis() }}:
  {The function to run this application.}
}
}

\seealso{
\code{\link{EzApp}}

\code{\link{EzDataset}}
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
