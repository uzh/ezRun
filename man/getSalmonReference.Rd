\name{getSalmonReference}
\alias{getSalmonReference}
\title{Builds and gets a salmon reference}
\usage{
  getSalmonReference(param)
}
\arguments{
  \item{param}{a list of parameters:
    \itemize{
      \item{transcriptFasta}{ an optional character specifying the path to a fasta file. If specified, the reference will be prepared using it.}
      \item{ezRef@refIndex}{ a character specifying the location of the index that is used in the alignment.}
      \item{ezRef@refFeatureFile}{ a character specifying the path to the annotation feature file (.gtf).}
      \item{ezRef@refFastaFile}{ a character specifying the path to the fasta file.}
    }
  }
}
\value{
  Returns a character specifying the file path to the salmon indexed reference genome.
}
\description{
Builds and gets a salmon reference by returning the file path to it.
}
\examples{
  \dontrun{
    GENOMES_ROOT <- "/home/gtan/analysis/gtan/referenceFiles"
    Sys.setenv("PATH"=paste("/usr/local/ngseq/packages/Aligner/Salmon/0.8.2/bin",
                            Sys.getenv("PATH"), sep=":"))
    param <- list()
    param[['transcriptTypes']] = 'protein_coding,rRNA'
    param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-08-11'
    param <- ezParam(param)
    getSalmonReference(param)
  }
}
\author{
Roman Briskine
}
