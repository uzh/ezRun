% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/app-mapping.R
\name{getBWAReference}
\alias{getBWAReference}
\title{Builds and gets a BWA reference}
\usage{
getBWAReference(param)
}
\arguments{
\item{param}{a list of parameters:
\itemize{
  \item{ezRef@refIndex}{ a character specifying the location of the index that is used in the alignment.}
  \item{ezRef@refBuildDir}{ a character specifying the directory of the reference build.}
  \item{ezRef@refFastaFile}{ a character specifying the file path to the fasta file.}
}}
}
\value{
Returns a character specifying the file path to the BWA indexed reference genome.
}
\description{
Builds and gets a BWA reference by returning the file path to it.
}
\examples{
ls = list('refBuild' = 'Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Version-2013-03-07')
param = ezParam(userParam = ls)
getBWAReference(param)
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
