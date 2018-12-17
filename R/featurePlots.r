###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch






#' Title
#'
#' @param gRanges 
#' @param bamFiles 
#' @param gtfFile 
#' @param height 
#' @param width 
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' require(GenomicRanges)
#' require(ezRun)
#' gtfFile = "/srv/GT/reference/Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31/Genes/genes.gtf"
#' myGeneRange = GRanges(seqnames="2", ranges=IRanges(start=1, end=1000, name="foo"))
#' plotLocusCoverageProfile(myGeneRange, bamFiles, gtfFile=gtfFile,
#'                          height=20, width=20)
#'}

plotLocusCoverageProfile = function(gRanges, bamFiles, gtfFile=NULL,
                                    height=10, width=20){
  require(Gviz)
  require(GenomicFeatures)
  require(GenomicRanges)
  require(S4Vectors)
  if (!is.null(gtfFile)){
    txdb = makeTxDbFromGFF(gtfFile, dataSource="FGCZ", taxonomyId = 1) #96061234")
    #saveDb(txdb, "tx.db")
    #txdb = loadDb("tx.db")
  } else {
    stop("gtfFile is required")
  }
  if (is.null(names(bamFiles))){
    names(bamFiles) = sub(".bam", "", basename(bamFiles))
  }
  options(ucscChromosomeNames=FALSE)
  
  ### the approach using alTrackObjects
  alTrackList = list()
  for (nm in names(bamFiles)){
    alTrackList[[nm]] <- AlignmentsTrack(bamFiles[nm], name=nm, isPaired = FALSE,
                                         type = c("coverage")) #, "sashimi"))
  }
  pdfFiles = character()
  grList = split(gRanges, seqnames(gRanges))
  chrom = names(grList)[1]
  for (chrom in names(grList)){
    gr = grList[[chrom]]
    if (length(gr) == 0){
      next
    }
    if (is.null(names(gr))){
      stop("genomic ranges must have names")
    }
    geneTrack = GeneRegionTrack(range=txdb, chrom=chrom, name="Gene Model", transcriptAnnotation="symbol")
    nm = names(gr)[1]
    for (i in 1:length(gr)){
      message(names(gr)[i])
      trackList = list(GenomeAxisTrack(), geneTrack)
      for (sm in names(alTrackList)){
        trackList[[sm]] = alTrackList[[sm]]
      }
      pf = paste0(names(gr)[i], "-coverage.pdf")
      pdf(file=pf, width=width, height=height)
      plotTracks(trackList, from=start(gr)[i], to=end(gr)[i], ylim=c(0,40),
                 lwd=1)
      dev.off()
      pdfFiles[names(gr)[i]] = pf
    }
  }
  return(pdfFiles)
}

plotLocusAverageCoverageProfile = function(gRanges, bamFiles, grouping=NULL, gtfFile=NULL,
                                           scalingFactors=NULL,
                                           height=10, width=20){
  require(Gviz)
  require(GenomicFeatures)
  require(GenomicRanges)
  if (!is.null(gtfFile)){
    txdb = makeTxDbFromGFF(gtfFile, dataSource="FGCZ", taxonomyId = 1) #96061234")
    #saveDb(txdb, "tx.db")
    #txdb = loadDb("tx.db")
  } else {
    stop("gtfFile is required")
  }
  if (is.null(names(bamFiles))){
    names(bamFiles) = sub(".bam", "", basename(bamFiles))
  }
  
  options(ucscChromosomeNames=FALSE)
  
  #### compute the coverage manually and normalize between samples
  coverageList = list()
  for (strand in c("+", "-")){
    coverageList[[strand]] = list() 
    for (i in 1:length(gRanges)){
      gene = names(gRanges)[i]
      message(i)
      regionStart = start(gRanges)[i]
      regionEnd = end(gRanges)[i]
      chrom = seqnames(gRanges)[i]
      coverageList[[strand]][[gene]] = ezMatrix(0.0, rows=names(bamFiles), cols=regionStart:regionEnd)
      for (nm in names(bamFiles)){
        ga = ezReadGappedAlignments(bamFiles[nm], seqname=chrom, start=regionStart - 100, end=regionEnd + 100, strand = strand)
        ## TODO average the base positions!!! 
        covVec = as.vector(coverage(ga)[[chrom]][regionStart:regionEnd])
        stopifnot(!is.na(covVec))
        coverageList[[strand]][[gene]][nm, ] = covVec
      }
    }
  }
  
  ### normalize and group the coverages
  groupColors = getSampleColors(unique(grouping), colorNames = unique(grouping))
  names(group) = names(bamFiles)
  normCoverageList = lapply(coverageList, function(cl){
    lapply(cl, function(x){
      xNorm = ezScaleColumns(t(x), scalingFactors[rownames(x)])
      averageColumns(xNorm, group[colnames(xNorm)])
    })
  })
  
  
  pdfFiles = character()
  grList = split(gRanges, seqnames(gRanges))
  chrom = names(grList)[1]
  for (chrom in names(grList)){
    gr = grList[[chrom]]
    if (length(gr) == 0){
      next
    }
    if (is.null(names(gr))){
      stop("genomic ranges must have names")
    }
    geneTrack = GeneRegionTrack(range=txdb, chrom=chrom, name="Gene Model", transcriptAnnotation="symbol")
    ids = symbol(geneTrack)
    for (i in 1:length(gr)){
      message(names(gr)[i])
      trackList = list(GenomeAxisTrack(), geneTrack)
      for (strand in c("+", "-")){
        x = t(normCoverageList[[strand]][[gene]])
        dTrack = DataTrack(
          data = x,  chromosome = chrom, start = as.integer(colnames(x)),
          end = as.integer(colnames(x)), genome="foo",
          name = paste("coverage", strand), col.line=groupColors, groups=rownames(x), type="l")
        trackList[[strand]] = dTrack
      }
      yMax = max(values(trackList[["+"]]), values(trackList[["-"]]))
      pf = paste0(names(gr)[i], "-avgCoverage.pdf")
      pdf(file=pf)
      plotTracks(trackList, from=start(gr)[i], to=end(gr)[i], legend=TRUE,
                 lwd=1, ylim=c(0, yMax))
      dev.off()
      pdfFiles[names(gr)[i]] = pf
    }
  }
  return(pdfFiles)
}

##' @title Gets transcripts coverage
##' @description Gets transcripts coverage.
##' @param chrom a character vector containing chromosome names.
##' @param gff an annotation data.frame in gtf or gff format.
##' @param reads an object of the class GAlignments.
##' @param strandMode a character specifying the mode of the strand.
##' @param ranges an object of the class GRanges.
##' @template roxygen-template
##' @return Returns a list of transcript coverage.
getTranscriptCoverage = function(chrom, gff, reads, strandMode="both"){
  ## chrom can be NULL, single chr or multiple chrs.
  
  #if (length(chrom) > 1){
  #  transcriptCov = unlist(lapply(chrom, getTranscriptCoverage, gff, reads, strandMode), recursive=FALSE)
  #  return(transcriptCov)
  #}
  
  gffExon = gff[gff$type == "exon", ]
  if(!is.null(chrom)){
    gffExon = gffExon[gffExon$seqid %in% chrom, ]
  }
  gffExon = gffExon[order(gffExon$start), ]
  exonRanges = gffToRanges(gffExon)
  ## This sorting is not quite necessary because the order above; 
  ## but for the future, data.frame gff is not ideal.
  exonRanges <- sort(exonRanges) 
  exonsByTx <- GenomicRanges::split(exonRanges, names(exonRanges))
  
  ## For some gene build, especially drosophila, there is trans-splicing genes
  ## The exons within one transcript contains different strands.
  exonsByTx <- exonsByTx[elementNROWS(unique(strand(exonsByTx))) == 1L]
  
  exonCov <- getRangesCoverage(unlist(exonsByTx, use.names=FALSE),
                               reads, strandMode=strandMode)
  #exonCov_old = getRangesCoverage(exonRanges, reads, strandMode=strandMode)
  
  #system.time(transcriptCov2 <- S4Vectors::split(exonCov, names(exonCov)))
  transcriptCov <- relist(exonCov, exonsByTx)
  # system.time(transcriptCov3 <- tapply(exonCov, names(exonCov), unlist))
  ## These two methods are still relative slow. Should be faster than old tapply implementation.
  ## split: 7023.148 seconds.
  ## relist: 6879 seconds.
  ## tapply: 7474.300 seconds.
  
  #system.time(transcriptCov_old <- tapply(exonCov, names(exonCov), function(exonCovList){
  # Rle(values=unlist(lapply(exonCovList, runValue)), lengths=unlist(lapply(exonCovList, runLength)))}, 
  #                       simplify=FALSE))
  # 12939.112 seconds
  
  transcriptCov <- lapply(transcriptCov, unlist, use.names=FALSE)
  # 169.468 seconds
  #if(length(transcriptCov) == 0L){
    ## This can happen when gff on chrom has 0 ranges.
  #  return(list())
  #}else{
  #  transcriptCov <- RleList(transcriptCov)
  #}
  txNegStrand <- unlist(unique(strand(exonsByTx))) == "-"
  stopifnot(length(txNegStrand) == length(exonsByTx)) 
  ## Only one strand from each transcript
  
  txNegStrand <- which(txNegStrand)
  transcriptCov <- revElements(transcriptCov, txNegStrand)
  #trStrand = gffExon$strand[match(names(transcriptCov), gffExon$transcript_id)]
  #indexNegStrand = which(trStrand == "-")
  #transcriptCov[indexNegStrand] = lapply(transcriptCov[indexNegStrand], rev)
  
  return(transcriptCov) ## as a list. RleList is unbelievably slow in loop!
}

##' @describeIn getTranscriptCoverage Gets the range coverage.
getRangesCoverage = function(ranges, reads, strandMode="both"){
  require(GenomicAlignments)
  #if(!is.null(chrom)){
  #  stopifnot(runValue(seqnames(ranges)) == chrom)
  #}
  if (length(ranges) == 0){
    return(RleList())
  }
  stopifnot(strandMode %in% c("both", "sense", "antisense"))
  #rangeCov = vector("list", length(ranges))
  if (strandMode == "both"){
    #if (is.null(chrom)){
    covChrom = coverage(reads)
    rangeCov <- covChrom[ranges]
    #  rangeCov = mapply(function(chr, s, e){covChrom[[chr]][s:e]},
    #                    as.character(seqnames(ranges)), start(ranges), end(ranges))
    #} else {
      # although this can be done in the same way above. But [[chrom]] first can speed up.
    #  covChrom = coverage(reads)[[chrom]]
    #  rangeCov = mapply(function(s,e){covChrom[s:e]}, start(ranges), end(ranges))
    #}
  } else {
    if (strandMode == "antisense"){
      strand(ranges) = flipStrand(strand(ranges))
    }
    #isPos = as.character(strand(ranges)) == "+"
    isPos <- strand(ranges) == "+"
    use = strand(reads) == "+"
    covPos = coverage(reads[use])
    rangeCovPos <- covPos[ranges[isPos]]  ## This is ultra-fast. 1 seconds..
    names(rangeCovPos) <- which(as.logical(isPos))
    #rangeCov[isPos] = mapply(function(chr,s,e){covPos[[chr]][s:e]}, as.character(seqnames(ranges[isPos])), start(ranges[isPos]), end(ranges[isPos]))  ## This is too slow! Takes 1046.812 seconds
    use = strand(reads)== "-"
    covNeg = coverage(reads[use])
    rangeCovNeg <- covNeg[ranges[!isPos]]
    names(rangeCovNeg) <- which(as.logical(!isPos))
    #rangeCov[!isPos] = mapply(function(s,e){covNeg[s:e]}, start(ranges[!isPos]), end(ranges[!isPos]))    
    rangeCov <- c(rangeCovPos, rangeCovNeg)
    ## reorder into the same order in ranges
    rangeCov <- rangeCov[order(as.integer(names(rangeCov)))]
  }
  names(rangeCov) <- names(ranges)
  return(rangeCov)
}
