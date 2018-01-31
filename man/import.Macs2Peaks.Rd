\name{import.Macs2Peaks}
\alias{import.Macs2Peaks}
\title{
  Import the MACS2 peak file
}
\description{
  Read the MACS2 peak file into a \code{GRanges}.
}
\usage{
  import.Macs2Peaks(con)
}
\arguments{
  \item{con}{
    \code{character}(1): the filepath of MACS2 peak file.
  }
}
\value{
  A \code{GRanges} object with extra columns.
}
\references{
  \url{https://github.com/taoliu/MACS}
}
\author{
  Ge Tan
}
\examples{
  \dontrun{
    bedFn <- "/srv/gstore/projects/p2578/MACS2_21885_2018-01-29--13-11-27/C8301US_peaks.bed"
    bedGR <- import.Macs2Peaks(bedFn)
  }
}
