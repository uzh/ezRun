% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gff.R
\name{getEnsemblTypes}
\alias{getEnsemblTypes}
\title{Gets the ensembl types}
\usage{
getEnsemblTypes(gff)
}
\arguments{
\item{gff}{an annotation data.frame in gtf or gff format.}
}
\value{
Returns a character vector containing the types of the elements or NULL if they are not provided in \code{gff}.
}
\description{
Gets the ensembl types of an annotation data.frame either directly from the source, the gene_biotype or the gene_type.
}
\examples{
param = ezParam()
gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
et = getEnsemblTypes(gtf)
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
