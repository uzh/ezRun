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
        \item{featureLevel}{ if equal to "gene" and the feature level of the dataset to "isoform", the rawdata will be passed to \code{aggregateCountsByGene()} before returning it.}
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
  \code{\link{loadCountDataset}}
}
\examples{
  \dontrun{
    ## smart-Seq2 platform
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
    
    ## 10x
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '100'
param[['scratch']] = '200'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['name']] = 'SCCount_QC'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['refFeatureFile']] = 'genes.gtf'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = ''
param[['min_genes']] = '500'
param[['min_counts']] = '10000'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2886/SCCountQC_29562_2018-09-28--14-43-58'
output = list()
output[['Name']] = 'B9P9_NPCHP11_NDd5'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = ''
output[['Static Report [Link]']] = 'p2886/SCCountQC_29562_2018-09-28--14-43-58/B9P9_NPCHP11_NDd5_SCCountQC/00index.html'
output[['Report [File]']] = 'p2886/SCCountQC_29562_2018-09-28--14-43-58/B9P9_NPCHP11_NDd5_SCCountQC'
input = list()
input[['Name']] = 'B9P9_NPCHP11_NDd5'
input[['ResultDir']] = 'p2886/CellRangerCount_29556_NOV19_o4775_2018-09-16--08-18-12/B9P9_NPCHP11_NDd5'
input[['Report']] = 'p2886/CellRangerCount_29556_NOV19_o4775_2018-09-16--08-18-12/B9P9_NPCHP11_NDd5/outs/web_summary.html'
input[['BAM']] = 'p2886/CellRangerCount_29556_NOV19_o4775_2018-09-16--08-18-12/B9P9_NPCHP11_NDd5/outs/possorted_genome_bam.bam'
input[['BAI']] = 'p2886/CellRangerCount_29556_NOV19_o4775_2018-09-16--08-18-12/B9P9_NPCHP11_NDd5/outs/possorted_genome_bam.bam.bai'
input[['Species']] = 'Homo sapiens (human)'
input[['reference']] = '10X_Ref_Human_GRCh38.p10_20180305_Release_91'
input[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
input[['refFeatureFile']] = 'genes.gtf'
input[['featureLevel']] = 'gene'
input[['CountMatrix']] = 'p2886/CellRangerCount_29556_NOV19_o4775_2018-09-16--08-18-12/B9P9_NPCHP11_NDd5/outs/filtered_gene_bc_matrices/10X_Ref_Human_GRCh38.p10_20180305_Release_91/matrix.mtx'
input[['Sample Id']] = 'bfs_184126'
input[['Order Id']] = '4775'
input = EzDataset$new(meta=input, dataRoot=param$dataRoot)
param <- ezParam(param)
param$scProtocol <- "10x"
output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
    
sce <- loadSCCountDataset(input, param)
  }
}
