\name{doGo}
\alias{doGo}
\alias{hasGoAnnotation}
\title{Do GO?}
\description{
Decides whether to do a gene ontologies analysis.
}
\usage{
  doGo(param, seqAnno)
  hasGoAnnotation(seqAnno)
}
\arguments{
  \item{param}{a list of parameters to extract the logical \code{runGO} and the character \code{featureLevel} from.}
  \item{seqAnno}{the sequence annotation.}
}
\value{
  \code{TRUE} or \code{FALSE}.
}

\examples{
  \dontrun{
    ## FALSE
    param = list()
    param[['refBuild']] = 'Apis_Mellifera/HGSC/build4.5/Annotation/Release_4.5-2017-08-15'
    param[['refFeatureFile']] = 'genes.gtf'
    param = ezParam(param)
    seqAnno <- ezFeatureAnnotation(param, dataFeatureType="gene")
    hasGoAnnotation(seqAnno)
    doGo(param, seqAnno)
    
    ## TRUE
    param = list()
    param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
    param[['refFeatureFile']] = 'genes.gtf'
    param = ezParam(param)
    seqAnno <- ezFeatureAnnotation(param, dataFeatureType="gene")
    hasGoAnnotation(seqAnno)
    doGo(param, seqAnno)
  }
}
\author{
Rehrauer, Hubert
Schmid, Peter
}
