\name{writeSCMM}
\alias{writeSCMM}
\alias{readSCMM}
\title{
  Read and write single cell sparse matrix
}
\description{
  Read and write gene cell sparse matrix for single cell experiment.
}
\usage{
  writeSCMM(x, file)
  readSCMM(file)
}
\arguments{
  \item{x}{
    \code{matrix} with rownames as gene ids and colnames as cells
  }
  \item{file}{
    \code{character}(1): the filename of output matrix with \sQuote{.mtx} as extension.
  }
}
\details{
  The corresponding colnames and rownames of the matrix are saved in \sQuote{.colNames}
  and \sQuote{.rowNames} files within the same folder of \sQuote{.mtx} folder.
}
\value{
  Invisible \code{file} for \code{writeSCMM}.
  
  \code{dgCMatrix} of the gene cell count matrix for \code{readSCMM}.
}
\author{
  Ge Tan
}
\examples{
  x <- matrix(c(2,0,0,0,0,1,0,0), ncol=2,
              dimnames=list(c("gene1", "gene2", "gene3", "gene4"),
                            c("cell1", "cell2")))
  file <- tempfile(pattern = "countMatrix", fileext = ".mtx")
  writeSCMM(x, file)
  x2 <- readSCMM(file)
}
