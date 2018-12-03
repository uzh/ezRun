\name{getRpkm}
\alias{getRpkm}
\alias{getTpm}
\alias{getSignal}
\alias{getCpm}
\title{Gets measurements from rawData}
\description{
  Gets the rpkm, tpm, signal or cpm measurement from \code{SummarizedExperiment} or 
  \code{SingleCellExperiment} object.
}

\usage{
  getRpkm(rawData)
  getTpm(rawData)
  getSignal(rawData)
  getCpm(rawData)
}
\arguments{
  \item{rawData}{
    \code{SummarizedExperiment} or \code{SingleCellExperiment} object.
  }
}
\value{
  Returns the \code{matrix} of rpkm, tpm, signal or cpm measurement.
}
\author{
  Ge Tan
}
\references{
  \url{https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/}
}
\examples{
  \dontrun{
    file <- "/srv/gstore/projects/p1688/CountQC_20201_2017-08-15--20-52-44/input_dataset.tsv"
    param <- ezParam(list(dataRoot="/srv/gstore/projects",
                          refBuild="Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31"))
    input <- EzDataset$new(file=file, dataRoot=param$dataRoot)
    rawData <- loadCountDataset(input, param)
    getRpkm(rawData)
    getTpm(rawData)
    getSignal(rawData)
    getCpm(rawData)
  }
}