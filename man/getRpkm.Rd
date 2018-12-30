\name{getRpkm}
\alias{getRpkm}
\alias{getTpm}
\alias{getSignal}
\alias{getCpm}
\title{Get measurements from rawData}
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
    ## bulk
    file <- "/srv/gstore/projects/p1688/CountQC_20201_2017-08-15--20-52-44/input_dataset.tsv"
    param <- ezParam(list(dataRoot="/srv/gstore/projects",
                          refBuild="Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31"))
    input <- EzDataset$new(file=file, dataRoot=param$dataRoot)
    rawData <- loadCountDataset(input, param)
    getRpkm(rawData)
    getTpm(rawData)
    getSignal(rawData)
    getCpm(rawData)
    
    ## single cell
    param <- ezParam(list(dataRoot="/srv/gstore/projects",
                          refBuild="Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26",
                          scProtocol="smart-Seq2"
                          ))
    input = list()
    input[['Name']] = 'NCC_Trpv5_Tomato_sc_A01'
    input[['Species']] = ''
    input[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
    input[['paired']] = 'false'
    input[['featureLevel']] = 'gene'
    input[['refFeatureFile']] = 'genes.gtf'
    input[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
    input[['ResultDir']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01'
    input[['CellDataset']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01-dataset.tsv'
    input[['CountMatrix']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01-counts.mtx'
    input[['CountMatrixColNames']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01-counts.colNames'
    input[['CountMatrixRowNames']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01-counts.rowNames'
    input[['Stats']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01-stats.txt'
    input[['CellCyclePhase']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01-CellCyclePhase.txt'
    input[['BAM']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01.bam'
    input[['BAI']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01.bam.bai'
    input[['PreprocessingLog']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01_preprocessing.log'
    input[['STARLog']] = 'p2838/SCCountsApp_28469_2018-12-05--15-25-05/NCC_Trpv5_Tomato_sc_A01/NCC_Trpv5_Tomato_sc_A01_STAR.log'
    input <- EzDataset$new(meta=input, dataRoot=param$dataRoot)
    rawData <- loadSCCountDataset(input, param)
    getTpm(rawData)
    getRpkm(rawData)
    getCpm(rawData)
  }
}