###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCountOverlaps = function(input = NA, output = NA, param = NA) {
  bamFile = input$getFullPaths("BAM")
  outputFile = basename(output$getColumn("Count"))

  require("bitops", warn.conflicts = WARN_CONFLICTS, quietly = !WARN_CONFLICTS)
  gff = ezLoadFeatures(param, types = "exon")
  if (ezIsSpecified(param$countNonredundant) && param$countNonredundant) {
    targets = countNonredundant(bamFile = bamFile, param = param, gff = gff) ## supports both single- and paired-end
    ezWrite.table(
      targets[order(rownames(targets)), , drop = FALSE],
      file = outputFile
    )
  } else {
    if (param$paired) {
      return(.countPairedBamHits(bamFile, outputFile, param))
    }
    if (ezIsSpecified(param$seqNames)) {
      chroms = param$seqNames
    } else {
      chroms = unique(gff$seqid)
    }
    chroms = intersect(ezBamSeqNames(input), chroms)
    names(chroms) = chroms
    message(
      "counting ",
      input,
      "; number of hits [Mio]: ",
      signif(sum(countBam(input)$records) / 1e6, digits = 3)
    )
    targetList = ezMclapply(
      chroms,
      FUN = countBamHitsSingleChrom,
      bamFile = input,
      param = param,
      gff = gff,
      mc.cores = getOption("cores")
    )
    targets = do.call("rbind", targetList)
    targetNames = unlist(lapply(targetList, rownames))
    if (any(duplicated(targetNames))) {
      targets = averageRows(targets, targetNames, func = sum)
    } else {
      rownames(targets) = targetNames ## undo the prepending of the list names
    }
    ezWrite.table(
      targets[order(rownames(targets)), , drop = FALSE],
      file = outputFile
    )
  }
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodCountOverlaps(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @section Functions:
##' \itemize{
##'   \item{\code{countNonredundant(bamFile, param=param, gff=gff): }}
##'   {Counts the non-redundant overlaps.}
##'   \item{\code{countBamHitsSingleChrom(chr, bamFile=NULL, param=NULL, gff=NULL): }}
##'   {Counts the BAM hits for a single chromosome.}
##'   \item{\code{countPairedBamHitsSingleChrom(chr, bamFile=NULL, param=NULL, gff=NULL): }}
##'   {Counts the paired BAM hits for a single chromosome.}
##'   \item{\code{getTargetRanges(gff, param, chrom=NULL): }}
##'   {Gets the target ranges depending on \code{param$featureLevel}.}
##'   \item{\code{getFeatureCounts(chrom, gff, reads, param): }}
##'   {Gets the feature counts from the target ranges.}
##' }
EzAppCountOverlaps <-
  setRefClass(
    "EzAppCountOverlaps",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodCountOverlaps
        name <<- "EzAppCountOverlaps"
        appDefaults <<- rbind(
          countNonredundant = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "make sure that every read counts only as one"
          ),
          minFeatureOverlap = ezFrame(
            Type = "integer",
            DefaultValue = "10",
            Description = "the number of bases overlap are need to generate a count"
          ),
          countTrimmedTranscripts = ezFrame(
            Type = "logical",
            DefaultValue = "FALSE",
            Description = "if counts for trimmed versions of the transcripts should be generated"
          ),
          trimmedMaxLength = ezFrame(
            Type = "integer",
            DefaultValue = '200',
            Description = "for trimmed counting use only that many bases at 3' and 5'-end respectively"
          )
        )
      }
    )
  )

.countPairedBamHits = function(input = NULL, output = NULL, param = NULL) {
  #param = fillWithDefaults(param) ## TODOMF: function doesn't exist
  options(cores = param$cores)
  message("countPairedBamHits")
  require("bitops", warn.conflicts = WARN_CONFLICTS, quietly = !WARN_CONFLICTS)

  gff = ezLoadFeatures(param)
  if (ezIsSpecified(param$seqNames)) {
    chroms = as.character(param$seqNames)
  } else {
    chroms = unique(gff$seqid)
  }
  chroms = intersect(ezBamSeqNames(input), chroms)
  names(chroms) = chroms
  gff = gff[gff$seqid %in% chroms, ]

  message(
    "counting ",
    input,
    "; number of hits [Mio]: ",
    signif(sum(countBam(input)$records) / 1e6, digits = 3)
  )
  targetList = ezMclapply(
    chroms,
    FUN = countPairedBamHitsSingleChrom,
    bamFile = input,
    param = param,
    gff = gff,
    mc.cores = getOption("cores")
  )
  targets = do.call("rbind", targetList) ## allegedly rbindlist would  be faster but it does not seem to exist in R-2.15.2; seems that it is in library(data.table)
  targetNames = unlist(lapply(targetList, rownames))
  if (any(duplicated(targetNames))) {
    targets = averageRows(targets, targetNames, func = sum)
  } else {
    rownames(targets) = targetNames ## undo the prepending of the list names
  }
  ezWrite.table(
    targets[order(rownames(targets)), , drop = FALSE],
    file = output
  )
  return("Success")
}

countNonredundant = function(bamFile, param = param, gff = gff) {
  if (!ezIsSpecified(param$minFeatureOverlap)) {
    param$minFeatureOverlap = 1L
  }
  targetRanges = getTargetRanges(gff, param, chrom = NULL)
  ga = ezReadGappedAlignments(
    bamFile,
    minMapQuality = param$minMapQuality,
    keepMultiHits = param$keepMultiHits
  )
  hits = findOverlaps(
    query = targetRanges,
    subject = ga,
    ignore.strand = param$strandMode == "both",
    minoverlap = param$minFeatureOverlap
  )
  weights = 1 / table(names(ga))
  counts = tapply(
    weights[match(names(ga)[subjectHits(hits)], names(weights))],
    queryHits(hits),
    FUN = sum
  )
  countFrame = data.frame(row.names = names(targetRanges))
  countFrame$matchCounts = counts[match(
    as.character(1:nrow(countFrame)),
    names(counts)
  )]
  if (param$countTrimmedTranscripts) {
    remove(hits)
    gc()
    gffStart = gffTrimTranscripts(
      gff,
      maxLength = param$trimmedMaxLength,
      start = TRUE
    )
    hits = findOverlaps(
      query = getTargetRanges(gffStart, param, chrom = NULL),
      subject = ga,
      ignore.strand = param$strandMode == "both",
      minoverlap = param$minFeatureOverlap
    )
    counts = tapply(
      weights[match(names(ga)[subjectHits(hits)], names(weights))],
      queryHits(hits),
      FUN = sum
    )
    nm = paste("matchCounts [first", param$trimmedMaxLength, "bases]")
    countFrame[[nm]] = counts[match(
      as.character(1:nrow(countFrame)),
      names(counts)
    )]
    remove(hits)
    gc()
    gffEnd = gffTrimTranscripts(
      gff,
      maxLength = param$trimmedMaxLength,
      start = FALSE
    )
    hits = findOverlaps(
      query = getTargetRanges(gffEnd, param, chrom = NULL),
      subject = ga,
      ignore.strand = param$strandMode == "both",
      minoverlap = param$minFeatureOverlap
    )
    counts = tapply(
      weights[match(names(ga)[subjectHits(hits)], names(weights))],
      queryHits(hits),
      FUN = sum
    )
    nm = paste("matchCounts [last", param$trimmedMaxLength, "bases]")
    countFrame[[nm]] = counts[match(
      as.character(1:nrow(countFrame)),
      names(counts)
    )]
  }
  countFrame[is.na(countFrame)] = 0
  return(countFrame)
}

countBamHitsSingleChrom = function(
  chr,
  bamFile = NULL,
  param = NULL,
  gff = NULL
) {
  targets = NULL
  tryCatch(
    {
      ga = ezReadGappedAlignments(
        bamFile,
        seqname = chr,
        minMapQuality = param$minMapQuality,
        keepMultiHits = param$keepMultiHits
      )
      featCounts = getFeatureCounts(chr, gff, ga, param)
      targets = data.frame(row.names = names(featCounts))
      targets$multiMatchCounts = as.numeric(featCounts)
      if (param$countTrimmedTranscripts) {
        gffStart = gffTrimTranscripts(
          gff,
          maxLength = param$trimmedMaxLength,
          start = TRUE
        )
        nm = paste("multiMatchCounts [first", param$trimmedMaxLength, "bases]")
        targets[[nm]] = getFeatureCounts(chr, gffStart, ga, param)[rownames(
          targets
        )]
        gffEnd = gffTrimTranscripts(
          gff,
          maxLength = param$trimmedMaxLength,
          start = FALSE
        )
        nm = paste("multiMatchCounts [last", param$trimmedMaxLength, "bases]")
        targets[[nm]] = getFeatureCounts(chr, gffEnd, ga, param)[rownames(
          targets
        )]
      }
      rm(ga)
      gc()
    },
    error = function(e) {
      ezWrite("Failed on: ", bamFile, " - ", chr)
      ezWrite(e)
    },
    traceback()
  )
  return(targets)
}

countPairedBamHitsSingleChrom = function(
  chr,
  bamFile = NULL,
  param = NULL,
  gff = NULL
) {
  targets = NULL
  tryCatch(
    {
      ga = ezReadPairedAlignments(
        bamFile,
        seqname = chr,
        keepUnpaired = param$keepUnpaired,
        minMapQuality = param$minMapQuality,
        keepMultiHits = param$keepMultiHits
      )
      featCounts = getFeatureCounts(chr, gff, ga, param)
      targets = data.frame(row.names = names(featCounts))
      targets$multiMatchCounts = as.numeric(featCounts)
      if (param$countTrimmedTranscripts) {
        gffStart = gffTrimTranscripts(
          gff,
          maxLength = param$trimmedMaxLength,
          start = TRUE
        )
        nm = paste("multiMatchCounts [first", param$trimmedMaxLength, "bases]")
        targets[[nm]] = getFeatureCounts(chr, gffStart, ga, param)[rownames(
          targets
        )]
        gffEnd = gffTrimTranscripts(
          gff,
          maxLength = param$trimmedMaxLength,
          start = FALSE
        )
        nm = paste("multiMatchCounts [last", param$trimmedMaxLength, "bases]")
        targets[[nm]] = getFeatureCounts(chr, gffEnd, ga, param)[rownames(
          targets
        )]
      }
      rm(ga)
      gc()
    },
    error = function(e) {
      ezWrite("Failed on: ", bamFile, " - ", chr)
      ezWrite(e)
      warnings()
      traceback()
      stop(e)
    }
  )
  return(targets)
}

getTargetRanges = function(gff, param, chrom = NULL) {
  require(S4Vectors)
  stopifnot(gff$type == "exon")
  if (!is.null(chrom)) {
    gff = gff[gff$seqid == chrom, ]
  }
  if (param$strandMode == "antisense") {
    gff$strand = flipStrand(gff$strand)
  }
  targetRanges = NULL
  if (!ezIsSpecified(param$featureLevel)) {
    param$featureLevel == "gene"
  }
  switch(
    param$featureLevel,
    exon = {
      gff$ID = paste(
        gff$transcript_id,
        sprintf("%03d", getExonNumber(gff)),
        sep = ":E"
      )
      targetRanges = gffToRanges(gff)
    },
    intron = {
      trRanges = gffGroupToRanges(gff, gff$transcript_id)
      exonRanges = gffToRanges(gff)
      exonsByTranscript = S4Vectors::split(exonRanges, gff$transcript_id)
      stopifnot(setequal(names(exonsByTranscript), names(trRanges)))
      intronsByTranscript = psetdiff(
        trRanges,
        exonsByTranscript[names(trRanges)]
      )
      targetRanges = unlist(intronsByTranscript)
      ids = names(targetRanges)
      names(targetRanges) = paste(
        ids,
        sprintf("%03d", ezReplicateNumber(ids)),
        sep = ":I"
      )
    },
    isoform = {
      targetRanges = S4Vectors::split(gffToRanges(gff), gff$transcript_id)
    },
    tss = {
      targetRanges = S4Vectors::split(gffToRanges(gff), gff$tss_id)
    },
    gene = {
      targetRanges = S4Vectors::split(gffToRanges(gff), gff$gene_id)
    }
  )
  if (is.null(targetRanges)) {
    stop("unsupported featureLevel: ", param$featureLevel)
  }
  return(targetRanges)
}

getFeatureCounts = function(chrom, gff, reads, param) {
  require(GenomicAlignments)
  if (length(chrom) > 1) {
    featCounts = unlist(
      ezMclapply(chrom, getFeatureCounts, gff, reads, param),
      recursive = FALSE
    )
    return(featCounts)
  }
  targetRanges = getTargetRanges(gff, param, chrom = chrom)

  if (!ezIsSpecified(param$minFeatureOverlap)) {
    param$minFeatureOverlap = 1L
  }
  featCounts = countOverlaps(
    query = targetRanges,
    subject = reads,
    ignore.strand = param$strandMode == "both",
    minoverlap = param$minFeatureOverlap
  )
  return(featCounts)
}

## generates counts when alignment is done to transcript database
## in that case the reference name (rname) is the transcript ID
## does handle paired-end consistently, i.e. a pair is counted as 1 if the ends align to the same  transcript
## and it is counted fractionally if alignments go to multiple transcripts
countTranscriptBam = function(
  bamFile,
  strand = "*",
  isFirstMateRead = NA,
  isSecondMateRead = NA,
  isUnmappedQuery = FALSE,
  isProperPair = NA,
  isSecondaryAlignment = NA
) {
  bam = ezScanBam(
    bamFile,
    what = c("qname", "rname"),
    strand = strand,
    isFirstMateRead = isFirstMateRead,
    isSecondMateRead = isSecondMateRead,
    isUnmappedQuery = isUnmappedQuery,
    isProperPair = isProperPair,
    isSecondaryAlignment = isSecondaryAlignment
  )
  multiMatch = table(bam$qname)
  counts = tapply(1 / multiMatch[bam$qname], bam$rname, sum)
  rm("bam", "multiMatch")
  gc()
  counts = counts[ezBamSeqNames(bamFile)]
  counts[is.na(counts)] = 0
  return(counts)
}
