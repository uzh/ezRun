% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf.R
\name{ezFilterVcf}
\alias{ezFilterVcf}
\alias{ezWriteVcf}
\title{Filters a VCF file and writes it back}
\usage{
ezFilterVcf(
  vcfFile,
  vcfFiltFile,
  discardMultiAllelic = TRUE,
  bamDataset = bamDataset,
  param = NULL
)

ezWriteVcf(vcf, vcfFiltFile)
}
\arguments{
\item{vcfFile}{a character representing the file path of the input VCF file}

\item{vcfFiltFile}{a character representing the file path of the output VCF file.}

\item{discardMultiAllelic}{a logical. Removes multi allelic variants if set to TRUE.}

\item{bamDataset}{a data.frame to help renaming and reordering the fields in genotype.}

\item{param}{a list of parameters to extract the value from \code{vcfFilt.minAltCount}.}

\item{vcf}{an object of the class VCF.}
}
\value{
Returns a filtered VCF file.
}
\description{
Filters a VCF file and writes it back to the specified connection.
}
\section{Functions}{
\itemize{
\item \code{ezWriteVcf()}: Writes a new VCF file.

}}
\author{
Rehrauer, Hubert

Schmid, Peter
}
