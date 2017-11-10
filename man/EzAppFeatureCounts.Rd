\name{EzAppFeatureCounts}
\alias{EzAppFeatureCounts}
\alias{ezMethodFeatureCounts}

\title{
  FeatureCounts App
}
\description{
  gene/Transcripts counting with featurecounts from Rsubread.
}
\usage{
  ezMethodFeatureCounts(input=input, output=output, param=param)
}
\arguments{
  \item{input}{
    \code{EzDataset} object
  }
  \item{output}{
    \code{EzDataset} object
  }
  \item{param}{
    Object returned by \code{ezParam}.
  }
}
\details{
  The features to count can be filtered by \code{param$transcriptTypes}.
  
  A count file with two columns: Identifier and matchCounts.
}
\author{
  Ge Tan
}

\examples{
  \dontrun{
param = list()
param[['cores']] = '8'
param[['ram']] = '20'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
param[['paired']] = 'false'
param[['strandMode']] = 'antisense'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['gtfFeatureType']] = 'exon'
param[['allowMultiOverlap']] = 'true'
param[['countPrimaryAlignmentsOnly']] = 'true'
param[['minFeatureOverlap']] = '10'
param[['minMapQuality']] = '10'
param[['keepMultiHits']] = 'true'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA,long_noncoding,short_noncoding,pseudogene'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2578/FeatureCounts_22125_2017-11-05--17-36-34'
output = list()
output[['Name']] = 'A0062P'
output[['Count [File]']] = 'p2578/FeatureCounts_22125_2017-11-05--17-36-34/A0062P.txt'
output[['Stats [File]']] = 'p2578/FeatureCounts_22125_2017-11-05--17-36-34/A0062P-stats.txt'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
output[['featureLevel']] = 'gene'
output[['refFeatureFile']] = 'genes.gtf'
output[['strandMode']] = 'antisense'
output[['paired']] = 'false'
output[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA,long_noncoding,short_noncoding,pseudogene'
input = list()
input[['Name']] = 'A0062P'
input[['BAM']] = 'p2578/STAR_22114_2017-11-04--20-57-11/A0062P.bam'
input[['BAI']] = 'p2578/STAR_22114_2017-11-04--20-57-11/A0062P.bam.bai'
input[['IGV Starter']] = 'p2578/STAR_22114_2017-11-04--20-57-11/A0062P-igv.jnlp'
input[['Species']] = 'Homo sapiens (human)'
input[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
input[['paired']] = 'false'
input[['refFeatureFile']] = 'genes.gtf'
input[['strandMode']] = 'antisense'
EzAppFeatureCounts$new()$run(input=input, output=output, param=param)

param[['aroundTSSCounting']] = 'true'
param[['upstreamFlanking']] = '500'
param[['downstreamFlanking']] = '500'
EzAppFeatureCounts$new()$run(input=input, output=output, param=param)


  }
}
