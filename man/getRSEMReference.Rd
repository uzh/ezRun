\name{getRSEMReference}
\alias{getRSEMReference}
\title{Builds and gets a RSEM reference}
\description{
Builds and gets a RSEM reference by returning the file path to it.
}
\usage{
  getRSEMReference(param)
}
\arguments{
\item{param}{a list of parameters:
\itemize{
  \item{transcriptFasta}{ an optional character specifying the path to a fasta file. If specified, the reference will be prepared using it.}
  \item{ezRef@refIndex}{ a character specifying the location of the index that is used in the alignment.}
  \item{ezRef@refFeatureFile}{ a character specifying the path to the annotation feature file (.gtf).}
  \item{ezRef@refFastaFile}{ a character specifying the path to the fasta file.}
  \item{transcriptTypes}{
  \code{character}(n): an optional list of transcript types to select to build the indices.
  }
}}
}
\value{
Returns a character specifying the file path to the RSEM indexed reference genome.
}

\examples{
  \dontrun{
    GENOMES_ROOT <- "/home/gtan/analysis/gtan/referenceFiles"
    Sys.setenv("PATH"=paste("/usr/local/ngseq/packages/Aligner/RSEM/1.3.0/bin",
                            "/usr/local/ngseq/packages/Aligner/Bowtie/1.1.2/bin",
                            Sys.getenv("PATH"), sep=":"))
    param <- list()
    param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
    param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-08-11'
    param <- ezParam(param)
    getRSEMReference(param)
  }
  ls = list('refBuild' = 'Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Version-2013-03-07')
  param = ezParam(userParam = ls)
  getRSEMReference(param)
}
\author{
Rehrauer, Hubert

Schmid, Peter
}
