\name{ezMethodSubsampleReads}
\alias{ezMethodSubsampleReads}
\title{Subsample reads in a fastq dataset file}
\usage{
  ezMethodSubsampleReads(input = NA, output = NA, param = NA)
}
\arguments{
\item{param}{a list of parameters where the following entries are used
\itemize{
  \item{dataRoot} the prefix of the file paths
  \item{nReads} the number of reads to keep; if given will be used to compute subsampleFactor; it is not guaranteed that the number of reads kept is exact
  \item{subsampleFactor} the factor by which subsampling has been done. if \code{nReads} is specified subsampleFactor will not be used
  \item{paired} whether these are paired-end reads
}}
}
\description{
The subsampled reads are equally distributed across the original files
}
\examples{
  inputDatasetFile = system.file(package = "ezRun", "extdata/yeast_10k/dataset.tsv")
  param = ezParam(list(dataRoot=system.file(package = "ezRun"), nReads=100))
  input = EzDataset(file=inputDatasetFile, dataRoot=param$dataRoot)
  xSubsampled = ezMethodSubsampleReads(input=input, param=param)
  file.remove(xSubsampled$getFullPaths("Read1"))
}
