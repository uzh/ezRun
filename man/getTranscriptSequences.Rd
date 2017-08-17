\name{getTranscriptSequences}
\alias{getTranscriptSequences}
\title{Gets transcript sequences}
\description{
Gets transcript sequences from annotation (gtf or gff) and sequence (fasta) information.
}
\usage{
  getTranscriptSequences(param=NULL, genomeFn=NULL, featureFn=NULL)
}
\arguments{
  \item{param}{
    the parameters to load the annotation and sequence files from.
  }
  \item{genomeFn}{
    \code{character(1)}: filepath of the genome assembly fasta file.
  }
  \item{featureFn}{
    \code{character(1)}: filepath of the gtf annotation file.
  }
}
\details{
 \sQuote{genomeFn} and \sQuote{featureFn} are determined from \sQuote{param}
 when available.
}
\value{
  Returns an object of the class \code{DNAStringSet}.
}

\examples{
  \dontrun{
    genomeFn <- "/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/WholeGenomeFasta/genome.fa"
    featureFn <- "/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf"
    ts1 <- getTranscriptSequences(genomeFn=genomeFn, featureFn=featureFn)
    
    param = ezParam()
    param$ezRef@refFeatureFile = featureFn
    param$ezRef@refFastaFile = genomeFn
    ts2 <- getTranscriptSequences(param)
  }
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
  