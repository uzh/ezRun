\name{ezFeatureAnnotation}
\alias{ezFeatureAnnotation}
\title{
  Feature annotation file IO
}
\description{
  Read feature annotation file, either from transcript/isoform level
  or gene level.
}
\usage{
  ezFeatureAnnotation(param, ids = NULL, 
                      dataFeatureType = c("gene", "transcript", "isoform"))
}
\arguments{
  \item{param}{
    A \code{list} of ezRun parameters or \code{character}(1) of filepath to *_annotation_byTranscript.txt.
  }
  \item{ids}{
    \code{character}(n) of gene_id or transcript_id to fetch.
    When NULL, all annotation are returned.
  }
  \item{dataFeatureType}{
    \code{character}(1): it can be \dQuote{gene}, \dQuote{transcript} or \dQuote{isoform}.
  }
}
\details{
  For each genome build, we have two types of gtf files: genes.gtf and featrues.gtf.
  For each gtf file, we have two correspodning feature annotation file: 
  *_annotation_byGene.txt and *_annotation_byTranscript.txt.
  The former has one gene per row and the latter has one transcript per row.
  
  When \code{dataFeatureType} is \sQuote{gene}, *_annotation_byGene.txt, 
  in the same folder of *_annotation_byTranscript.txt, is used when available.
  Otherwise, *_annotation_byTranscript.txt is used with aggregation.
}
\value{
  \code{data.frame} of feature annotation.
}
\author{
  Ge Tan
}
\examples{
  param <- list()
  param[['refBuild']] <- 'Saccharomyces_cerevisiae/Ensembl/R64/Annotation/Release_98-2019-12-03'
  param[['refFeatureFile']] <- 'genes.gtf'
  param <- ezParam(param)
  txAnno <- ezFeatureAnnotation(param, dataFeatureType="transcript")
  geneAnno <- ezFeatureAnnotation(param, dataFeatureType="gene")
}