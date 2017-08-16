###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodFeatureCounts = function(input=NA, output=NA, param=NA){
  bamFile = input$getFullPaths("BAM")
  localBamFile = .getBamLocally(bamFile)
  outputFile = basename(output$getColumn("Count"))
  statFile = basename(output$getColumn("Stats"))
  
  
  if (ezIsSpecified(param$transcriptTypes)){
    gtfFile = "genes.gtf"
    seqAnno = ezRead.table(param$ezRef@refAnnotationFile)
    transcriptsUse = rownames(seqAnno)[seqAnno$type %in% param$transcriptTypes]
    require(data.table)
    require(stringr)
    ## fread is faster than read_tsv, read.table and import.
    ## some benchmarks for human gtf on Mac i7-3615QM
    ### read_tsv: 27.662s
    ### ezRead.table: 127.461s
    ### import: 33.419s
    gtf <- fread(param$ezRef@refFeatureFile, header=FALSE, sep="\t", 
                 quote="", data.table=FALSE,
                 colClasses=c("character", "character", "character", "integer", 
                              "integer", "character", "character", "character",
                              "character"),
                 col.names=c('seqname','source','feature','start','end',
                             'score','strand','frame','attribute')
                 )
    #transcripts <- sub("\";$", "", 
    #                   sub("^transcript_id \"", "", 
    #                       str_extract(gtf$attribute,
    #                                   "transcript_id \"[[:alnum:][:punct:]]+\";")))
    transcripts <- ezGffAttributeField(gtf$attribute,
                                       field="transcript_id", 
                                       attrsep="; *", valuesep=" ")
    gtf = gtf[transcripts %in% transcriptsUse, ]
    ## write.table is much faster than write_tsv and export.
    ### write_tsv: 38.270s
    ### export: 347.998s
    ### write.table: 8.787s
    write.table(gtf, gtfFile, quote=FALSE, sep="\t", 
                row.names=FALSE, col.names=FALSE)
  } else {
    gtfFile = param$ezRef@refFeatureFile
  }
  
  sink(file="featureCounts-messages.txt")
  if(!is.null(param$aroundTSSCounting) && param$aroundTSSCounting){
    ## TODO:  refactor below; you can not rely that "gene" entries are in the gtf file; GenomicFeatures has functions for flanking
    ## - load the gtf in a txdb, use the promoter function to get the promoters
    ## - read transcript type from the genes_annotation.txt
    ## - filter on transcript type
    ## - provide the promoter coordinates as a data.frame to featureCounts 
    gtf = rtracklayer::import(param$ezRef@refFeatureFile)
    idx = gtf$type == 'gene'
    gtf = gtf[idx]
    
    if (ezIsSpecified(param$transcriptTypes)){
      idx = gtf$gene_biotype %in% param$transcriptTypes
      gtf = gtf[idx]
    }
    
    rtracklayer::export(gtf, 'genes.gtf')
    
    gtf <- fread('genes.gtf', header=FALSE)
    colnames(gtf) = c('seqname','source','feature','start','end','score','strand','frame','attribute')
    for (i in 1:nrow(gtf)){
      if(gtf[['strand']][i] == '+'){
        gtf[['end']][i] = max(0,gtf[['start']][i] + param$downstreamFlanking)
        gtf[['start']][i] = max(0,gtf[['start']][i] - param$upstreamFlanking)
      } else {
        gtf[['start']][i] = max(0,gtf[['end']][i] - param$downstreamFlanking)
        gtf[['end']][i] = max(0,gtf[['end']][i] + param$upstreamFlanking)
      }
    }
    write.table(gtf, gtfFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    countResult = Rsubread::featureCounts(localBamFile, annot.inbuilt=NULL,
                                          annot.ext=gtfFile, isGTFAnnotationFile=TRUE,
                                          GTF.featureType='gene',
                                          GTF.attrType= 'gene_id',                                          
                                          useMetaFeatures=param$useMetaFeatures,
                                          allowMultiOverlap=param$allowMultiOverlap, isPairedEnd=param$paired, 
                                          requireBothEndsMapped=FALSE,
                                          checkFragLength=FALSE,minFragLength=50,maxFragLength=600,
                                          nthreads=param$cores, 
                                          strandSpecific=switch(param$strandMode, "both"=0, "sense"=1, "antisense"=2, stop("unsupported strand mode: ", param$strandMode)),
                                          minMQS=param$minMapQuality,
                                          readExtension5=0,readExtension3=0,read2pos=NULL,
                                          minOverlap=param$minFeatureOverlap,
                                          ignoreDup=param$ignoreDup,
                                          splitOnly=FALSE,
                                          countMultiMappingReads=param$keepMultiHits,
                                          fraction=param$keepMultiHits & !param$countPrimaryAlignmentsOnly,
                                          primaryOnly=param$countPrimaryAlignmentsOnly,
                                          countChimericFragments=TRUE,chrAliases=NULL,reportReads=FALSE)
  }
  else{
    countResult = Rsubread::featureCounts(localBamFile, annot.inbuilt=NULL,
                              annot.ext=gtfFile, isGTFAnnotationFile=TRUE,
                              GTF.featureType=param$gtfFeatureType,
                              GTF.attrType=switch(param$featureLevel,
                                                  "gene"="gene_id",
                                                  "transcript"="transcript_id",
                                                  "isoform"="transcript_id",
                                                  stop("unsupported feature level: ", param$featureLevel)),
                              useMetaFeatures=param$useMetaFeatures,
                              allowMultiOverlap=param$allowMultiOverlap, isPairedEnd=param$paired, 
                              requireBothEndsMapped=FALSE,
                              checkFragLength=FALSE,minFragLength=50,maxFragLength=600,
                              nthreads=param$cores, 
                              strandSpecific=switch(param$strandMode, "both"=0, "sense"=1, "antisense"=2, stop("unsupported strand mode: ", param$strandMode)),
                              minMQS=param$minMapQuality,
                              readExtension5=0,readExtension3=0,read2pos=NULL,
                              minOverlap=param$minFeatureOverlap,
                              ignoreDup=param$ignoreDup,
                              splitOnly=FALSE,
                              countMultiMappingReads=param$keepMultiHits,
                              fraction=param$keepMultiHits & !param$countPrimaryAlignmentsOnly,
                              primaryOnly=param$countPrimaryAlignmentsOnly,
                              countChimericFragments=TRUE,chrAliases=NULL,reportReads=FALSE)
  }
  sink(file=NULL)
  
  colnames(countResult$counts) = "matchCounts"
  ezWrite.table(countResult$counts, file=outputFile)
  colnames(countResult$stat) = c("Status", "Count")
  ezWrite.table(countResult$stat, file=statFile, row.names=FALSE)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFeatureCounts(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppFeatureCounts <-
  setRefClass("EzAppFeatureCounts",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodFeatureCounts
                  name <<- "EzAppFeatureCounts"
                  appDefaults <<- rbind(gtfFeatureType = ezFrame(Type="character",  DefaultValue="exon",  Description="which gtf feature types to use; with Ensembl GTF files; use 'transcript' to count also intronic reads"),
                                        allowMultiOverlap = ezFrame(Type="logical",  DefaultValue="TRUE",  Description="make sure that every read counts only as one"),
                                        countPrimaryAlignmentsOnly = ezFrame(Type="logical",  DefaultValue="TRUE",  Description="count only the primary alignment"),
                                        minFeatureOverlap=ezFrame(Type="integer", DefaultValue="10", Description="the number of bases overlap are need to generate a count"),
                                        useMetaFeatures=ezFrame(Type="logical", DefaultValue="TRUE", Description="should counts be summarized to meta-features"),
                                        ignoreDup=ezFrame(Type="logical", DefaultValue="FALSE", Description="ignore reads marked as duplicates"),
                                        aroundTSSCounting=ezFrame(Type="logical",  DefaultValue="FALSE",  Description="count reads only around TSS of genes"),
                                        downstreamFlanking=ezFrame(Type="integer", DefaultValue="250", Description="the number of bases downstream of TSS"),
                                        upstreamFlanking=ezFrame(Type="integer", DefaultValue="250", Description="the number of bases upstream of TSS"))
                }
              )
  )
