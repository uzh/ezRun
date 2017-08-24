\name{getTranscriptCoverage}
\alias{getTranscriptCoverage}
\alias{getRangesCoverage}
\title{Gets transcripts coverage}
\description{
  Gets transcripts coverage.
}
\usage{
  getTranscriptCoverage(chrom, gff, reads, strandMode = "both")
  getRangesCoverage(chrom = NULL, ranges, reads, strandMode = "both")
}
\arguments{
\item{chrom}{a character vector containing chromosome names.}
\item{gff}{an annotation data.frame in gtf or gff format.}
\item{reads}{an object of the class GAlignments.}
\item{strandMode}{a character specifying the mode of the strand.}
\item{ranges}{an object of the class GRanges.}
}
\value{
  Returns a list of transcript coverage.
}
\section{Functions}{
\itemize{
\item \code{getRangesCoverage}: Gets the range coverage.
}}

\author{
Rehrauer, Hubert
Schmid, Peter
}

\examples{
  \dontrun{
    bamFn <- "/srv/gstore/projects/p2438/STAR_18564_2017-06-12--13-46-30/26EV_BL_A.bam"
    reads = ezReadGappedAlignments(bamFn)
    gff <- ezLoadFeatures(featureFile="/srv/GT/reference/Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Genes/genes.gtf")
    gffExon = gff[gff$type == "exon", ]
    gffExon = gffExon[order(gffExon$start), ]
    exonRanges = gffToRanges(gffExon)
    exonRanges <- sort(exonRanges)
    exonsByTx <- GenomicRanges::split(exonRanges, names(exonRanges))
    exonCov <- getRangesCoverage(unlist(exonsByTx, use.names=FALSE),
                                 reads, strandMode="antisense")
    transcriptCov <- getTranscriptCoverage(chrom=NULL, gff, reads,
                                           strandMode="antisense")
  }
}