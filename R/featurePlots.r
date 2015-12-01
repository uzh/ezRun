###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets transcripts coverages
##' @description Gets transcripts coverages.
##' @param chrom a character vector containing chromosome names.
##' @param gff an annotation data.frame in gtf or gff format.
##' @param reads an object of the class GAlignments.
##' @param strandMode a character specifying the mode of the strand.
##' @param ranges an object of the class GRanges.
##' @template roxygen-template
##' @return Returns a list of transcript coverages.
getTranscriptCoverage = function(chrom, gff, reads, strandMode="both"){
  if (length(chrom) > 1){
    transcriptCov = unlist(lapply(chrom, getTranscriptCoverage, gff, reads, strandMode), recursive=FALSE)
    return(transcriptCov)
  }
  gffExon = gff[gff$type == "exon", ]
  if(!is.null(chrom)){
    gffExon = gffExon[gffExon$seqid == chrom, ]
  }
  gffExon = gffExon[order(gffExon$start), ]
  exonRanges = gffToRanges(gffExon)
  exonCov = getRangesCoverageChrom(chrom, exonRanges, reads, strandMode=strandMode)
  transcriptCov = tapply(exonCov, gffExon$transcript_id, function(exonCovList){
    Rle(values=unlist(lapply(exonCovList, runValue)), lengths=unlist(lapply(exonCovList, runLength)))}, 
                         simplify=FALSE)
  trStrand = gffExon$strand[match(names(transcriptCov), gffExon$transcript_id)]
  indexNegStrand = which(trStrand == "-")
  transcriptCov[indexNegStrand] = lapply(transcriptCov[indexNegStrand], rev)
  return(transcriptCov)
}

##' @describeIn getTranscriptCoverage Gets the range coverages.
getRangesCoverageChrom = function(chrom=NULL, ranges, reads, strandMode="both"){
  if(!is.null(chrom)){
    stopifnot(runValue(seqnames(ranges)) == chrom)
  }
  if (length(ranges) == 0){
    return(list())
  }
  stopifnot(strandMode %in% c("both", "sense", "antisense"))
  rangeCov = vector("list", length(ranges))
  if (strandMode == "both"){
    if (is.null(chrom)){
      covChrom = coverage(reads)
      rangeCov = mapply(function(chr, s, e){covChrom[[chr]][s:e]},
                        as.character(seqnames(ranges)), start(ranges), end(ranges))
    } else {
      # although this can be done in the same way above. But [[chrom]] first can speed up.
      covChrom = coverage(reads)[[chrom]]
      rangeCov = mapply(function(s,e){covChrom[s:e]}, start(ranges), end(ranges))
    }
  } else {
    if (strandMode == "antisense"){
      strand(ranges) = flipStrand(strand(ranges))
    }
    isPos = as.character(strand(ranges)) == "+"
    use = strand(reads)== "+"
    covPos = coverage(reads[use])[[chrom]]
    rangeCov[isPos] = mapply(function(s,e){covPos[s:e]}, start(ranges[isPos]), end(ranges[isPos]))
    use = strand(reads)== "-"
    covNeg = coverage(reads[strand(reads)== "-"])[[chrom]]
    rangeCov[!isPos] = mapply(function(s,e){covNeg[s:e]}, start(ranges[!isPos]), end(ranges[!isPos]))    
  }
  names(rangeCov) = names(ranges)
  return(rangeCov)
}
