##' @title Gets the transcript annotation
##' @description Gets the transcript annotation from a gtf annotation.
##' @param gtf an annotation data.frame of the gtf format.
##' @param id a column of the gtf data.frame that specifies a factor to group transcripts by.
##' @template types-template
##' @param attributes additional annotation attributes to pass in.
##' @template roxygen-template
##' @return Returns a transposed transcript annotation data.frame.
##' @examples
##' param = ezParam()
##' gtfFile = system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE)
##' gtf = ezLoadFeatures(param, gtfFile)
##' ta = transcriptAnnoFromGtf(gtf)
transcriptAnnoFromGtf = function(gtf, id=gtf$transcript_id, types="exon", attributes=NULL){
  stopifnot(length(id) == nrow(gtf))
  use = gtf$type %in% types
  gtf = gtf[use, ]
  id = id[use]
  seqAnno = data.frame(row.names=unique(id))
  seqAnno$gene_id = tapply(gtf$gene_id, id, ezCollapse, empty.rm=TRUE, uniqueOnly=TRUE)[rownames(seqAnno)] 
  seqAnno$seqid = tapply(gtf$seqid, id, ezCollapse, empty.rm=TRUE, uniqueOnly=TRUE)[rownames(seqAnno)]
  seqAnno$strand = tapply(gtf$strand, id, ezCollapse, empty.rm=TRUE, uniqueOnly=TRUE)[rownames(seqAnno)]
  seqAnno$start = tapply(gtf$start, id, ezCollapse)[rownames(seqAnno)]        
  seqAnno$end = tapply(gtf$end, id, ezCollapse)[rownames(seqAnno)]        
  seqAnno$width = tapply(gtf$end - gtf$start + 1, id, sum)[rownames(seqAnno)]
  if (!is.null(attributes)){
    for (attr in attributes){
      aValues = ezGffAttributeField(x=gtf$attributes, field=attr, attrsep="; *", valuesep=" ")
      seqAnno[[attr]] = tapply(aValues, id, ezCollapse, na.rm=TRUE, empty.rm=TRUE, uniqueOnly=TRUE)[rownames(seqAnno)]
    }
  }
  return(seqAnno)
}

##' @title Adds transcripts from exons to the annotation
##' @description Adds transcripts from exons to the annotation using a gtf or gff data.frame.
##' @param gff the annotation data.frame to get transcripts from.
##' @template roxygen-template
##' @return Returns an annotation data.frame with transcripts added to it.
##' @examples
##' param = ezParam()
##' gtf = ezLoadFeatures(param, system.file("extdata/genes.gtf", package="ezRun", mustWork=TRUE))
##' gtf2 = addTranscriptsToGffExons(gtf)
addTranscriptsToGffExons = function(gff){
  
  gffExon = gff[gff$type == "exon", ]
  if (is.null(gffExon$transcript_id)){
    gffExon$transcript_id = ezGffAttributeField(gffExon$attributes, field="transcript_id", attrsep="; *", valuesep=" ")
  }
  gffExon$id = paste(gffExon$transcript_id, gffExon$seqid, gffExon$strand)
  gffTranscript = gffExon[!duplicated(gffExon$id), ]
  trStart = by(gffExon$start, gffExon$id, min)
  trEnd = by(gffExon$end, gffExon$id, max)
  gffTranscript$start = trStart[gffTranscript$id]
  gffTranscript$end = trEnd[gffTranscript$id]
  gffTranscript$type = "transcript"
  gffTranscript$transcript_id = NULL
  gffTranscript$id = NULL
  gffExon$transcript_id = NULL
  gffExon$id = NULL
  gffAll = rbind(gffTranscript, gffExon, gff[gff$type != "exon", ])
  return(gffAll)
}