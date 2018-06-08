\name{fastqs2bam}
\alias{fastqs2bam}
\title{
  Fastq files to Bam with RG tag
}
\description{
  Convert multiple fastq files into Bam files and merged into one big Bam file
  with filenames as RG tag.
}
\usage{
  fastqs2bam(fastqFns, fastq2Fns=NULL, readGroupNames=NULL, bamFn,
             platform="illumina", mc.cores=ezThreads())
}
\arguments{
  \item{fastqFns}{
    \code{character}(n): paths of input fastq files
  }
  \item{fastq2Fns}{
    \code{character}(n): optional paths of the second read of paired-end data.
  }
  \item{readGroupNames}{
    \code{character}(n): read group IDs. Same length as \code{fastqFns}.
  }
  \item{bamFn}{
    \code{character}(1): filename of merged Bam file.
  }
  \item{platform}{
    \code{character}(1): the platform used.
  }
  \item{mc.cores}{
    \code{integer}(1): the threads to run.
  }
}
\details{
  Each fastq file first is converted into Bam file with filename
  as \dQuote{LB} and \code{readGroupNames} as \sQuote{ID} and \sQuote{SM} in RG tag 
  using picard's FastqToSam.
  Then Bam files are merged into one big Bam file using picard's MergeSamFiles
  with a sort order by queryname.
}
\value{
  invisible \code{bamFn}.
}
\author{
  Ge Tan
}

\seealso{
  \code{\link{bam2fastq}}
}
\examples{
  \dontrun{
    fastqFns <- list.files(path="/srv/gstore/projects/p2288/HiSeq2500_20171011_RR99_o3511/dmx",
                           pattern="\\\.fastq\\\.gz$", full.names=TRUE)
    fastqFns <- fastqFns[1:5]
    fastqs2bam(fastqFns,
               readGroupNames=sub("\\\.(fastq|fq)(\\\.gz){0,1}$", "", 
                                  basename(fastqFns)),
               bamFn="20171011.A-C1_HT_24H_unmapped.bam")
    
    ## Make uBAM from cell dataset.tsv
    input <- EzDataset(file="/srv/gstore/projects/p2488/HiSeq4000_20180209_RUN424_o3954/scFastq_dataset.tsv", 
                       dataRoot="/srv/gstore/projects")
    fastqs2bam(fastqFns=input$getFullPaths("Read1"),
               readGroupNames=input$getNames(),
               bamFn="20180209.B-single_cell_hippo_neurons_E1_1_unmapped.bam",
               mc.cores = 4)
  }
}
