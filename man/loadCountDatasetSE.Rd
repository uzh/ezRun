\name{loadCountDatasetSE}
\alias{loadCountDatasetSE}
\title{
  Loads the count dataset
}
\description{
  Loads the count dataset with the given input and param into a
  \code{SummarizedExperiment} object.
}
\usage{
  loadCountDatasetSE(input, param)
}

\arguments{
  \item{input}{
    A \code{EzDataset} object with \sQuote{dataRoot} and filepath of 
    \sQuote{input_dataset.tsv}.
  }
  \item{param}{
    A \code{list} of parameters:
      \itemize{
        \item{dataRoot}{ the root directory of the files.}
        \item{expressionName}{ if specified, this will be used as the column name...}
        \item{knownExpressionNames}{ ...or otherwise known expression names that occur in the dataset will be used.}
        \item{ezRef@refBuild}{ if specified, the sequence annotation will be extracted from \code{ezFeatureAnnotation()}.}
        \item{transcriptTypes}{ if specified, only the defined transcript type will be used.}
        \item{sigThresh}{ the threshold...}
        \item{useSigThresh}{ ...and whether it should be used.}
        \item{featureLevel}{ if equal to "gene" and the feature level of the dataset to "isoform", the rawdata will be passed to \code{aggregateCountsByGeneSE()} before returning it.}
      }
  }
}
\value{
  A \code{SummarizedExperiment} object with
    \itemize{
      \item{assays: }{counts, presentFlag, RPKM, TPM, (signal)}
      \item{rowData: }{seqAnno}
      \item{colData: }{input$meta}
      \item{metadata: }{isLog, featureLevel, type, countName, param}
    }
}
\references{
  \url{https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html}
}
\author{
  Ge Tan
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  \dontrun{
    file <- "/srv/gstore/projects/p1688/CountQC_20201_2017-08-15--20-52-44/input_dataset.tsv"
    param <- ezParam(list(dataRoot="/srv/gstore/projects",
                          refBuild="Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31"))
    input <- EzDataset$new(file=file, dataRoot=param$dataRoot)
    countsSE1 <- loadCountDatasetSE(input, param)
    
  }
}
