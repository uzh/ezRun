###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


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


getRangesProfiles = function(ranges, strandMode=NA, bamFiles=NULL){
  sampleList = vector("list", length(bamFiles))
  names(sampleList) = names(bamFiles)
  covList = lapply(1:length(ranges), function(x){sampleList})
  names(covList) = names(ranges)
  .getCov = function(chrom, rgs=NULL, bamFile=NULL, strandMode=NULL){
    message(bamFile, " ", chrom)
    ## TODO: does not consider the parameters minMapQuality, keepMultiHits
    reads = ezReadGappedAlignments(bamFile, seqname=chrom)
    rgCov = getRangesCoverageChrom(chrom, rgs[seqnames(rgs) == chrom], reads, strandMode=strandMode)
    return(rgCov)
  }
  chromNames = intersect(ezBamSeqNames(bamFiles[1]), seqlevels(ranges))
  
  for(nm in names(bamFiles)){
    rgCov = unlist(ezMclapply(chromSet, .getCov, rgs=ranges, bamFile=bamFiles[nm], strandMode=strandMode, mc.preschedule=FALSE))
    covList[names(rgCov)] = mapply(function(x,y){x[[nm]]= y; return(x)}, 
                                   covList[names(rgCov)], rgCov, SIMPLIFY=FALSE, USE.NAMES=TRUE)
  }
  return(covList)
}


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
    if(is.null(chrom)){
      covChrom = coverage(reads)
      rangeCov = mapply(function(chr, s, e){covChrom[[chr]][s:e]},
                        as.character(seqnames(ranges)), start(ranges), end(ranges))
    }else{
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



plotRangesProfiles = function(ranges, rangesProfiles, trdb, pngDir=".", makePng=TRUE, width=50, offset=10, param=param){
  if (makePng){
    if (!file.exists(pngDir)) dir.create(pngDir)
  }
  for (chrom in unique(seqnames(ranges))){
    idx = which(seqnames(ranges) == chrom)
    grTrack = GeneRegionTrack(range=trdb, chrom=chrom)  
    for (i in idx){
      prof = RangesProfiles[[i]]
      dat = sapply(prof, function(x){as.numeric(runmean(x, width))[seq(1, length(x), by=width)]})
      avgDat = rowMeans(dat)
      logRatio = log2((dat + offset) /(avgDat + offset))
      regionStart = start(ranges)[i]
      regionEnd = end(ranges)[i]
      regionPos = seq(regionStart, regionEnd, by=width)    
      trackList = list(grTrack)
      lrShrinked = cbind(seq(-param$logColorRange, param$logColorRange, length.out=ncol(logRatio)), 
                         shrinkToRange(t(logRatio), c(-param$logColorRange, param$logColorRange)))
      dTrack <- DataTrack(data = lrShrinked[rev(colnames(dat)), ], start = c(regionPos[1], regionPos),
                          end = c(regionPos[1], regionPos+width-1), chromosome = chrom, gradient=ezRedBlueScale(256), nrcolors=256,
                          name = names(tssRangesUse[i]), type="heatmap")
      trackList[["heat"]] = dTrack
      for (nm in colnames(dat)){
        dTrack <- DataTrack(data = rbind(avgDat, dat[, nm ]), start = regionPos,
                            end = regionPos+width-1, chromosome = chrom,
                            name = nm, col.line=c("black", "blue"), groups=c("avg", "sample"), type="l")
        trackList[[nm]] = dTrack      
      }
      if (makePng){
        png(file=paste0(pngDir, "/", names(ranges)[i], "-profiles.png"), height=800, width=1600)    
      }
      plotTracks(trackList, from=regionStart, to=regionEnd, lwd=1)
      if (makePng){
        dev.off()
      }
    }
  }
}


getTranscriptProfiles = function(transcriptIds, bamFile, gtf=NULL, strandMode=NA,
                                 nThreads=ezThreads()){
  if (is.null(gtf$transcript_id)){
    gtf$transcript_id = ezGffAttributeField(gtf$attributes, field="transcript_id", attrsep="; *", valuesep=" ")    
  }
  stopifnot(transcriptIds %in% gtf$transcript_id)
  gtfUse = gtf[gtf$transcript_id %in% transcriptIds & gtf$type == "exon", ]
  trByChrom = tapply(gtfUse$transcript_id, gtfUse$seqid, unique)
  .getCov = function(chrom){
    ## TODO: does not consider the parameters minMapQuality, keepMultiHits
    reads = ezReadGappedAlignments(bamFile, seqname=chrom)
    return(getTranscriptCoverage(chrom, gtfUse, reads, strandMode=strandMode))
  }
  trCounts = sapply(trByChrom, length)
  trCounts = sort(trCounts, decreasing=TRUE)
  covByChrom = ezMclapply(names(trCounts), .getCov, mc.cores=nThreads, , mc.preschedule=FALSE)
  return(unlist(covByChrom, recursive=FALSE)[transcriptIds])
}



plotTranscriptProfiles = function(transcriptProfiles, gtf, sampleColors=NULL, makePng=TRUE,
                                  pngDir=".", pngNames=NULL){
  
  if (is.null(pngNames) && makePng){
    pngNames=paste0(pngDir, "/", names(transcriptProfiles), "-profiles.png")
  }
  if (is.null(sampleColors) & length(transcriptProfiles) > 0){
    sampleColors = rainbow(length(transcriptProfiles[[1]]))
    names(sampleColors) = names(transcriptProfiles[[1]])
  }
  gtfUse = gtf[gtf$transcript_id %in% names(transcriptProfiles) & gtf$type == "exon", ]
  ## transcript must not have features on different chromosomes or strands
  stopifnot(tapply(paste(gtfUse$seqid, gtfUse$strand), gtfUse$transcript_id, function(x){length(unique(x))}) == 1)
  if (is.null(gtfUse$width)){
    gtfUse$width = gtfUse$end - gtfUse$start +1
  }
  transcriptWidths = split(gtfUse$width, gtfUse$transcript_id)
  transcriptStrand = tapply(gtfUse$strand, gtfUse$transcript_id, unique)
  transcriptSeqid = tapply(gtfUse$seqid, gtfUse$transcript_id, unique)
  transcriptStart = tapply(gtfUse$start, gtfUse$transcript_id, min)
  transcriptEnd = tapply(gtfUse$end, gtfUse$transcript_id, max)
  for(i in 1:length(transcriptProfiles)){
    tid = names(transcriptProfiles)[i]
    profiles = transcriptProfiles[[i]]
    if (makePng){
      png(file=pngNames[i], width=1200, height=800)
    }
    xlim = c(1, length(profiles[[1]]))
    ylim = c(0, max(sapply(profiles, max)))
    par(mar=c(7,4,4,2)+0.1)
    main = paste(tid, transcriptSeqid[tid], transcriptStrand[tid], transcriptStart[tid], transcriptEnd[tid])
    plot(1, 1, type="n", xlim=xlim, ylim=ylim, main=main, xlab="", ylab="coverage", xaxt="n",cex.lab=1)
    title(xlab="exon start pos",mgp=c(5,1,0), cex.lab=1)
    axisAt = cumsum(transcriptWidths[[tid]])
    abline(v=axisAt)
    axis(side=1, at=axisAt, labels=axisAt, las=2)
    for(nm in names(profiles)){
      lines(profiles[[nm]], col=sampleColors[nm])
    }
    legend("topleft", legend=names(profiles), col=sampleColors[names(profiles)], lty=1)
    if (makePng){
      dev.off()
    }
  }  
}
