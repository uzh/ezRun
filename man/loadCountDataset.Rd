\name{loadCountDataset}
\alias{loadCountDataset}
\title{
  Loads the count dataset
}
\description{
  Loads the count dataset with the given input and param into a
  \code{SummarizedExperiment} object.
}
\usage{
  loadCountDataset(input, param)
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
        \item{featureLevel}{ if equal to "gene" and the feature level of the dataset to "isoform", the rawdata will be passed to \code{aggregateCountsByGene()} before returning it.}
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
  \code{\link{loadSCCountDataset}}
}
\examples{
  \dontrun{
    ## gene-level counts; gene-level return
    file <- "/srv/gstore/projects/p2828/CountQC_28750_2018-08-13--15-38-33/input_dataset.tsv"
    param <- ezParam(list(dataRoot="/srv/gstore/projects",
                          refBuild="Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26"))
    input <- EzDataset$new(file=file, dataRoot=param$dataRoot)
    rawData <- loadCountDataset(input, param)
    
    ## isoform-level counts; gene-level return
    file <- "/srv/gstore/projects/p2828/CountQC_29754_2018-09-24--15-15-29/input_dataset.tsv"
    param <- ezParam(list(dataRoot="/srv/gstore/projects",
                          refBuild="Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26",
                          featureLevel="gene",
                          expressionName=""))
    input <- EzDataset$new(file=file, dataRoot=param$dataRoot)
    rawData <- loadCountDataset(input, param)
    
    ## isoform-level counts; isoform-level return
    file <- "/srv/gstore/projects/p2828/CountQC_29754_2018-09-24--15-15-29/input_dataset.tsv"
    param <- ezParam(list(dataRoot="/srv/gstore/projects",
                          refBuild="Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26",
                          featureLevel="isoform",
                          expressionName=""))
    input <- EzDataset$new(file=file, dataRoot=param$dataRoot)
    rawData <- loadCountDataset(input, param)
  }
}
