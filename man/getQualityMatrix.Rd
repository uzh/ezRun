\name{getQualityMatrix}
\alias{getQualityMatrix}

\title{
  Get quality matrix
}
\description{
  Get quality matrix from FastQ file.
}
\usage{
  getQualityMatrix(fn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{
    \code{character(1)}: path of a fastq file.
  }
}
\details{
  Subsample 1\% of the reads to get the statistics.
}
\value{
  A \code{matrix} of occurency of each quality per base.
  Rows are quality value and columns are base positions.
}
\author{
  Ge Tan
}
\seealso{
  
}
\examples{
  \dontrun{
    fn <- "/srv/gstore/projects/p2438/HiSeq4000_20170601_RUN346_copy_o3365/20170601.B-19TK_d7_R1.fastq.gz"
    getQualityMatrix(fn)
  }
}
