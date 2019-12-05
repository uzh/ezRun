\name{ezCombineReadDatasets}
\alias{ezCombineReadDatasets}
\title{Combine the reads from multiple datasets into a single dataset}
\usage{
  ezCombineReadDatasets(..., dataRoot = "/srv/gstore/projects", newDsDir = NULL)
}
\arguments{
  \item{...}{\code{tibble} objects from the dataset.tsv files.}
  \item{dataRoot}{a character specifying the file root of the datasets.}
  \item{newDsDir}{a character specifying the directory to save the new dataset in.}
  \item{ds1}{a data.frame from the meta field of an EzDataset.}
  \item{ds2}{a data.frame from the meta field of an EzDataset.}
}
\description{
  Takes the union of the samples in multiple input datasets and generates a new dataset.
}
\details{
  If a sample is present in more than one datasets, 
  the read files are concatenated and a new file is written.
  
  If a sample is present in only one dataset it is simply copied.
  
  The \sQuote{Read Count} column must be present and is updated if two files are combined.
}
\value{
  A \code{tibble} object of new dataset.
}
\author{
  Ge Tan
}
\examples{
  library(readr)
  ds1 <- read_tsv(system.file("extdata/yeast_10k/dataset.tsv",
                              package = "ezRun", mustWork = TRUE))
  ds3 <- ds2 <- ds1
  dataRoot <- path.package("ezRun")
  newDsDir <- "scratch"
  ds <- ezCombineReadDatasets(ds1, ds2, ds3, dataRoot=dataRoot,
                              newDsDir=newDsDir)
}