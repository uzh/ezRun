\name{makeExtraControlSeqGR}
\alias{makeExtraControlSeqGR}
\alias{getControlSeqs}
\title{
  Functions for control sequences
}
\description{
  Fetch the control sequences and build a \code{GRanges} object from it.
}
\usage{
  makeExtraControlSeqGR(ids = NULL)
  getControlSeqs(ids = NULL)
}
\arguments{
  \item{ids}{
    A \code{character}(n) obejct: the chromosomes/sequences ids to build.
  }
}
\details{
  We have the control sequences available at \url{https://fgcz-gstore.uzh.ch/reference/controlSeqs.fa}.
}
\value{
  A \code{DNAStringSet} object from \code{getControlSeqs} and
  a \code{GRanges} obejct from \code{makeExtraControlSeqGR}
}
\author{
  Ge Tan
}

\examples{
  makeExtraControlSeqGR()
  getControlSeqs()
}
