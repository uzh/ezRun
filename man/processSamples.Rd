% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/app-exceRpt_report.R
\name{processSamples}
\alias{processSamples}
\title{generate smallRNA counts and QC.}
\usage{
processSamples(samplePaths, outputDir)
}
\arguments{
\item{samplePaths}{<string> vector of all the files that need to be processed together.}

\item{outputDir}{<string> directory name where to output generated objects and report after processing.}

\item{getPlotsObjects}{<bool> TRUE to output plots as an object.}
}
\description{
Obtain counts and QC plots across all samples.
}
\author{
(Rob Kitchen) Miquel Anglada Girotto
}
