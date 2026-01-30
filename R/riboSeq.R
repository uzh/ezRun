###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## see also scripts in the svn for project p1368

getTranslationInitiationSite = function(txdb) {
  cdsRgList = cdsBy(txdb, by = "tx", use.names = TRUE)
  cds = unlist(cdsRgList) ## unlist the GRangesList
  cdsFirst = cds[!duplicated(names(cds))] ## use only the first cds
  tis = resize(cdsFirst, width = 1, fix = "start")

  exonRgList = exonsBy(txdb, by = "tx", use.names = TRUE)
  exonRgList = exonRgList[names(tis)]
  stopifnot(names(exonRgList) %in% names(tis))
  tisMapped = pmapToTranscripts(tis, exonRgList)
  tisPos = start(tisMapped)
  names(tisPos) = names(tisMapped)
  return(tisPos)
}


getCdsProfiles = function(
  bamFile,
  tisPos,
  strand = "+",
  readLength = 32,
  offset = 14
) {
  require(GenomicAlignments)
  aln = readGAlignments(bamFile)
  stopifnot(names(tisPos) %in% levels(seqnames(aln)))
  ezLog(
    "genes with TIS but without alignments: ",
    length(setdiff(names(tisPos), as.vector(seqnames(aln))))
  )
  aln = aln[
    as.vector(seqnames(aln)) %in%
      names(tisPos) &
      as.vector(strand(aln)) == strand &
      qwidth(aln) %in% readLength
  ]
  codonStartBase = start(aln) + offset - tisPos[as.vector(seqnames(aln))] + 1 ## TODO profile may extend beyond end of CDS
  use = codonStartBase > 0
  codonStartRgs = GRanges(
    seqnames = seqnames(aln)[use],
    ranges = IRanges(start = codonStartBase[use], width = 1)
  )
  cdsProfile = coverage(codonStartRgs)
  #   ## check if we picked the right offset
  #   foo = sapply(cdsProfile[sapply(cdsProfile, length) > 20], function(x){as.vector(x[1:20])})
  #   plot(log2(rowSums(foo) + 1), type="both")

  return(cdsProfile)
}


getTranscriptProfiles = function(bamFile, param) {
  # strand="+", readLength=32, readStartOnlyCoverage=TRUE, getCoverageByReadlength=TRUE, minRead){
  require(GenomicAlignments)
  aln = readGAlignments(bamFile)
  aln = aln[qwidth(aln) %in% param$minReadLength:param$maxReadLength]
  if (param$strandMode == "sense") {
    aln = aln[as.vector(strand(aln)) == "+"]
  }
  if (param$strandMode == "antisense") {
    aln = aln[as.vector(strand(aln)) == "-"]
  }
  profileList = list()
  profileList[["all"]] = switch(
    param$coverageType,
    readStart = coverage(GRanges(
      seqnames = seqnames(aln),
      ranges = IRanges(start = start(aln), width = 1)
    )),
    fullRead = coverage(aln)
  )
  if (param$getCoverageByReadLength) {
    stopifnot((param$maxReadLength - param$minReadLength) %in% 1:50)
    for (readLength in param$minReadLength:param$maxReadLength) {
      alnUse = aln[qwidth(aln) == readLength]
      profileList[[paste("length", readLength)]] =
        switch(
          param$coverageType,
          readStart = coverage(GRanges(
            seqnames = seqnames(alnUse),
            ranges = IRanges(start = start(alnUse), width = 1)
          )),
          fullRead = coverage(alnUse)
        )
    }
  }
  return(profileList)
}


getExonProfiles = function(
  bamFile,
  strand = "+",
  readLength = 32,
  offset = 14
) {
  require(GenomicAlignments)
  aln = readGAlignments(bamFile)
  aln = aln[as.vector(strand(aln)) == strand & qwidth(aln) %in% readLength]
  codonStartBase = start(aln) + offset
  codonStartRgs = GRanges(
    seqnames = seqnames(aln),
    ranges = IRanges(start = codonStartBase, width = 1)
  )
  exonProfile = coverage(codonStartRgs)
  #   ## check if we picked the right offset
  #   foo = sapply(cdsProfile[sapply(cdsProfile, length) > 20], function(x){as.vector(x[1:20])})
  #   plot(log2(rowSums(foo) + 1), type="both")

  return(exonProfile)
}


## this function is redundant.
## could be achieved more elegantly be extraing the TIS region from the result of the getExonProfiles function

getAvgTisReadStartProfiles = function(
  bamFile,
  tisPos,
  tisRegion = c(-32, 99),
  maxExpression = NA,
  strand = "+",
  readLengths = 26:35
) {
  require(GenomicAlignments)
  aln = readGAlignments(bamFile)
  stopifnot(names(tisPos) %in% levels(seqnames(aln)))
  ezLog(
    "genes with TIS but without alignments: ",
    length(setdiff(as.vector(seqnames(aln)), names(tisPos)))
  )
  aln = aln[
    as.vector(seqnames(aln)) %in%
      names(tisPos) &
      as.vector(strand(aln)) == strand
  ]
  readStartExon = start(aln) ## in exon coordinates
  readStartCds = readStartExon - tisPos[as.vector(seqnames(aln))] + 1
  readWidth = qwidth(aln)
  ## potential further filters ???
  ## - do deduplicate
  ## - use only genes with low coverage
  #   isDup = duplicated(paste(alGene, alStrand, readWidth, alPosCds))
  #   isHighExpress = alGene %in% names(senseCounts)[senseCounts > 500]
  #   table(isHighExpress)
  #   table(isDup, isNearStart)
  isInTisRegion = readStartCds >= tisRegion[1] & readStartCds <= tisRegion[2]

  tisProfiles = ezMatrix(
    0,
    rows = tisRegion[1]:tisRegion[2],
    cols = readLengths
  )
  for (w in colnames(tisProfiles)) {
    use = readWidth == as.integer(w) & isInTisRegion
    ih = intHist(
      readStartCds[use],
      range = c(tisRegion[1] - 0.5, tisRegion[2] + 0.5),
      plot = FALSE
    )
    tisProfiles[as.character(ih$mids), w] = ih$counts
  }
  return(tisProfiles)
}
