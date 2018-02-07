\name{getDupRateFromBam}
\alias{getDupRateFromBam}
\title{
  Duplication rate from Bam of RNA-Seq.
}
\description{
  Read in a BAM file and count the tags falling on the features described in the GTF file
}
\usage{
  getDupRateFromBam(bamFile, param = NULL, gtfFn,
                    stranded = c("both", "sense", "antisense"),
                    paired = FALSE, threads = 1)
}
\arguments{
  \item{bamFile}{
    \code{character}(1): filename of bamFile.
  }
  \item{param}{
    \pkg{ezRun} param list.
  }
  \item{gtfFn}{
    \code{character}(1): The gtf file describing the features
  }
  \item{stranded}{
    \code{character}(1): \dQuote{both}, \dQuote{sense} or \dQuote{antisense}
  }
  \item{paired}{
    \code{logical}: Paired end experiment?
  }
  \item{threads}{
    \code{Integer}: The number of threads to be used for counting.
  }
}
\details{
  It uses the implementation of package \pkg{dupRadar}.
}
\value{
  A \code{data.frame} of counts on features, with and without taking into account multimappers/duplicated reads
}
\references{
  Sayols S, Scherzinger D and Klein H (2016). "dupRadar: a Bioconductor package for the assessment of PCR artifacts in RNA-Seq data." BMC Bioinformatics, 17, pp. 428. doi: 10.1186/s12859-016-1276-2, http://dx.doi.org/10.1186/s12859-016-1276-2.
}
\author{
  Ge Tan
}

\examples{
  \dontrun{
    bamFile <- "/srv/gstore/projects/p2438/STAR_18564_2017-06-12--13-46-30/26EV_BL_A.bam"
    gtfFn <- "/srv/GT/reference/Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31/Genes/genes.gtf"
    stranded <- "antisense"
    paired <- FALSE
    threads <- 4
    resBam <- getDupRateFromBam(bamFile=bamFile, gtfFn=gtfFn,
                                    stranded=stranded,
                                    paired=paired, threads=threads)
  }
}
