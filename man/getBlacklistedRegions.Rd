\name{getBlacklistedRegions}
\alias{getBlacklistedRegions}
\title{
  Get blacklisted genomic regions
}
\description{
  Get the black listed genomic regions for functional genomics analysis.
}
\usage{
getBlacklistedRegions(refBuild = c("hg38", "hg19", "mm10", "mm9", "ce10", "dm3"), 
                      format = c("UCSC", "Ensembl"))
}

\arguments{
  \item{refBuild}{
    \code{character}(1): one of the genome build of \dQuote{hg38}, \dQuote{hg19},
    \dQuote{mm10}, \dQuote{mm9}, \dQuote{ce10}, \dQuote{dm3}.
  }
  \item{format}{
    \code{character}(1): the returned \code{GRanges} with chromosomes in
    \dQuote{UCSC} or \dQuote{Ensembl} format.
  }
}
\details{
  These tracks were generated as part of ENCODE, modENCODE projects.
}
\value{
  \code{GRanges} of blacklisted regions.
}
\references{
  \url{https://sites.google.com/site/anshulkundaje/projects/blacklists}
}
\author{
  Ge Tan
}
\examples{
  getBlacklistedRegions("hg38")
  
  getBlacklistedRegions("hg38", format="Ensembl")
}
