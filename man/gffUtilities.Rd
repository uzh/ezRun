\name{gff Utilities}
\alias{ezGffAttributeField}
\title{
  gff utilities functions
}
\description{
  gff attributes parsing
}
\usage{
  ezGffAttributeField(x, field, attrsep = ";", valuesep = "=")
}
\arguments{
  \item{x}{
    \code{character}(n): attribute columns of gtf/gff file.
  }
  \item{field}{
    \code{character}(1): the field in attribute to extract
  }
  \item{attrsep}{
    \code{character}(1): the separation character in attribute. 
      \dQuote{;} in gff file and \dQuote{; } in gtf file.
  }
  \item{valuesep}{
    \code{character}(1): the separation character between values in attribute.
      \dQuote{=} in gff file and \dQuote{ } in gtf file.
  }
}
\value{
  A vector of \code{character} for the \code{field}.
}
\references{
  \url{http://www.ensembl.org/info/website/upload/gff.html}
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{ezReadGff}}
}
\examples{
  \dontrun{
    ## gtf
    gtf <- ezReadGff(system.file("extdata/genes.gtf", package="ezRun", 
                                 mustWork=TRUE))
    transcript_id <- ezGffAttributeField(gtf$attributes, field="transcript_id", 
                                         attrsep="; *", valuesep=" ")
    gene_id <- ezGffAttributeField(gtf$attributes, field="gene_id", 
                                   attrsep="; *", valuesep=" ")
    
    ## gff
    gff <- ezReadGff("/srv/GT/reference/Ca.Chlamydia_sanzinia_2742-308/NCBI/CP014639.1/Annotation/Version-2016-11-18/Genes/CP014639.gff3")
    IDs <- ezGffAttributeField(gff$attributes, field="ID", 
                               attrsep=";", valuesep="=")
  }
}
