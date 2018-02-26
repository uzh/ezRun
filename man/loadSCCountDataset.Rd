\name{loadSCCountDataset}
\alias{loadSCCountDataset}
\title{
  Loads the count dataset of single cell experiment
}
\description{
  Loads the count dataset with the given input and param into a
  \code{SingleCellExperiment} object.
}
\usage{
  loadSCCountDataset(input, param)
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
\details{
   A \code{SingleCellExperiment} object with
    \itemize{
      \item{assays: }{counts, presentFlag, RPKM, TPM}
      \item{rowData: }{seqAnno}
      \item{colData: }{cell dataset}
      \item{metadata: }{isLog, featureLevel, type, param}
    }
}
\value{
  
}
\references{
  \url{https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html#create-instances-of-singlecellexperiment}
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{loadCountDatasetSE}}  
}
\examples{
  \dontrun{
    library(ezRun)
    param = list()
    param[['cores']] = '1'
    param[['ram']] = '2'
    param[['scratch']] = '10'
    param[['node']] = ''
    param[['process_mode']] = 'DATASET'
    param[['name']] = 'SCCount_QC'
    param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
    param[['refFeatureFile']] = 'genes.gtf'
    param[['featureLevel']] = 'gene'
    param[['mail']] = 'ge.tan@fgcz.ethz.ch'
    param[['dataRoot']] = '/srv/gstore/projects'
    param[['resultDir']] = 'p2497/p2497-SCCountQC'
    output = list()
    output[['Name']] = 'SCCount_QC'
    output[['Species']] = ''
    output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
    output[['Static Report [Link]']] = 'p2497/p2497-SCCountQC/SCCount_QC/00index.html'
    output[['Report [File]']] = 'p2497/p2497-SCCountQC/SCCount_QC'
    input = '/srv/gstore/projects/p2497/SCCountsApp_23745_2018-01-05--11-55-48/dataset.tsv'
    
    input = EzDataset$new(file=input, dataRoot=param$dataRoot)
    param <- ezParam(param)
    output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
    
    sce <- loadSCCountDataset(input, param)
  }
}
