\name{ezLoadFeatures}
\alias{ezLoadFeatures}
\alias{addPromotersToGff}
\title{Loads annotation features from a file}
\usage{
  ezLoadFeatures(param = NULL, featureFile = param$ezRef["refFeatureFile"], 
                 types = NULL)
  addPromotersToGff(gff, promWidth)
}
\arguments{
\item{param}{contains the feature file and possibly a logical called \code{addPromoters}, which will add promoters if set to true.}

\item{featureFile}{the file to load the features from.}

\item{types}{a character vector containing the types.}
}
\value{
Returns a data.frame of parsed features.
}
\description{
Loads annotation features from a file and returns them as a data.frame.
}
\section{Functions}{
\itemize{
\item \code{ezGffAttributeField}: Gets the attribute from the specified \code{field}.

\item \code{addPromotersToGff}: Adds promoters to the gff list if \code{addPromoters} is specified.
}}

\examples{
  \dontrun{
    param = ezParam()
    gtfFile = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
    gtf = ezLoadFeatures(param, gtfFile)
    attrField = ezGffAttributeField(gtf$attributes, field="transcript_id", 
                                    attrsep="; *", valuesep=" ")
    addPromotersToGff(gtf, promWidth=300L)
  }
}
\seealso{
  \code{\link{ezReadGff}}
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
