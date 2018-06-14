\name{countReadsInFastq}
\alias{countReadsInFastq}
\title{Counts reads in FastQ files}
\usage{
countReadsInFastq(fastqFiles)
}
\arguments{
\item{fastqFiles}{a character vector representing file paths to FastQ files.}
}
\value{
Returns a named numeric vector
}
\description{
Counts reads in FastQ files and returns the amount for each file.
}
\examples{
file = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
param = list()
param$dataRoot = system.file(package="ezRun", mustWork = TRUE)
input = EzDataset$new(file=file, dataRoot=param$dataRoot)
fqFiles = input$getFullPaths("Read1")
result = countReadsInFastq(fqFiles)
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
