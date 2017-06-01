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
  
  
  if (ezIsSpecified(param$useTranscriptType)){
    gtfFile = "genes.gtf"
    seqAnno = ezRead.table(param$ezRef@refAnnotationFile)
    transcriptsUse = rownames(seqAnno)[seqAnno$type %in% param$useTranscriptType]
    require(readr)
    require(stringr)
    ## read_tsv is faster than read.table and import.
    gtf <- read_tsv(param$ezRef@refFeatureFile, comment="#", quote="", 
                    quoted_na=FALSE, na="", col_names=FALSE, 
                    col_types=cols(col_character(), col_character(),
                                   col_character(), col_character(),
                                   col_character(), col_character(),
                                   col_character(), col_character(),
                                   col_character()
                                   )
                    )
    transcripts <- sub("\";$", "", 
                       sub("^transcript_id \"", "", 
                           str_extract(gtf$X9, 
                                       "transcript_id \"[[:alnum:]]+\";")))
    
    gtf = gtf[transcripts %in% transcriptsUse, ]
    ## as.data.frame() is necessary as read_tsv quotes the attribute column.
    ## write.table is much faster than write_tsv and export.
    write.table(as.data.frame(gtf), gtfFile, quote=FALSE, sep="\t", 
                row.names=FALSE, col.names=FALSE)
  } else {
    gtfFile = param$ezRef@refFeatureFile
  }
  
  sink(file="featureCounts-messages.txt")
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
                                        ignoreDup=ezFrame(Type="logical", DefaultValue="FALSE", Description="ignore reads marked as duplicates"))
                }
              )
  )
