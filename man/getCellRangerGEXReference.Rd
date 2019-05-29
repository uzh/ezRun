\name{getCellRangerGEXReference}
\alias{getCellRangerGEXReference}
\alias{getCellRangerVDJReference}
\title{
  Build CellRanger index
}
\description{
  Build CellRanger GEX or VDJ index
}
\usage{
  getCellRangerGEXReference(param)
  getCellRangerVDJReference(param)
}

\arguments{
  \item{param}{
    list from \code{ezParam} .
  }
}
\value{
  \code{character}(1): path of index.
}
\author{
  Ge Tan
}

\examples{
setwd("/scratch/gtan/dev/quickdev")
param = list()
param[['cores']] = '8'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
param[['refFeatureFile']] = 'genes.gtf'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
param[['TenXLibrary']] = 'GEX'
param[['scMode']] = 'SC'
param <- ezParam(param)
getCellRangerGEXReference(param)

setwd("/scratch/gtan/dev/quickdev")
param = list()
param[['cores']] = '8'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
param[['refFeatureFile']] = 'genes.gtf'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
param[['TenXLibrary']] = 'GEX'
param[['scMode']] = 'SN'
param <- ezParam(param)
getCellRangerGEXReference(param)

setwd("/scratch/gtan/dev/quickdev")
param = list()
param[['cores']] = '8'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
param[['refFeatureFile']] = 'genes.gtf'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
param[['TenXLibrary']] = 'VDJ'
param <- ezParam(param)
getCellRangerVDJReference(param)

}
