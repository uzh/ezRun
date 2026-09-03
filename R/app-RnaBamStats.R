###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodRnaBamStats = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  require("GenomicAlignments")
  require("S4Vectors")

  samples = input$getNames()
  files = input$getFullPaths("BAM")
  dataset = input$meta

  setwdNew(basename(output$getColumn("Report")))
  param$featureLevel = "gene"
  param$projectId = sub("\\/.*", "", input$getColumn("BAM")[1]) ## project id is needed for the IGV link

  gff = ezLoadFeatures(param)
  if (!is.null(gff) && nrow(gff) == 0) {
    writeErrorReport(
      htmlFile,
      param = param,
      error = paste(
        "No features found in given feature file:<br>",
        param$ezRef["refFeatureFile"]
      )
    )
    return("Error")
  }

  ## get the RNA_repeats if available
  refRepeatsFeat = file.path(
    GENOMES_ROOT,
    param$ezRef["refBuild"],
    "Repeats/RNA_repeats.gff"
  )
  if (file.exists(refRepeatsFeat)) {
    repeatsGff = ezReadGff(refRepeatsFeat)
  } else {
    repeatsGff = NULL
  }

  if (ezIsSpecified(param$seqNames)) {
    gff = gff[gff$seqid %in% param$seqNames, ]
    if (!is.null(repeatsGff)) {
      repeatsGff = repeatsGff[repeatsGff$seqid %in% param$seqNames, ]
    }
  }
  resultList = list()
  for (sm in samples) {
    ezLog(sm)
    resultList[[sm]] = getStatsFromBam(
      param,
      files[sm],
      sm,
      gff = gff,
      repeatsGff = repeatsGff,
      nReads = dataset[sm, "Read Count"]
    )
    if (isError(resultList[[sm]])) {
      writeErrorReport(htmlFile, param = param, error = resultList[[sm]]$error)
      return()
    }
    print(gc())
  }

  if (is.null(param$posErrorRates) || param$posErrorRates == TRUE) {
    errorRates = ezMclapply(
      files,
      getPosErrorFromBam,
      param,
      mc.preschedule = FALSE,
      mc.cores = min(length(files), ezThreads())
    )
    for (sm in samples) {
      resultList[[sm]][["ErrorRates"]] = errorRates[[sm]]
    }
    rm(errorRates)
    gc()
  }
  ## the junction statistics come from STAR's SJ.out.tab if it is available for every
  ## sample; otherwise they are computed from the bam files with rseqc
  sjFiles = getUpstreamQcFiles(input, "Junctions")
  if (length(sjFiles) == length(samples)) {
    ezLog("using the splice junctions reported by the aligner")
    spliceSites = getAnnotatedSpliceSites(gff) ## must not be a promise; mclapply would force it in every child
    junctionsResults = ezMclapply(
      sjFiles,
      getJunctionStatsFromSJ,
      spliceSites,
      mc.preschedule = FALSE,
      mc.cores = min(length(sjFiles), ezThreads())
    )
  } else {
    junctionsResults = ezMclapply(
      files,
      getJunctionPlotsFromBam,
      param,
      mc.preschedule = FALSE,
      mc.cores = min(length(files), ezThreads())
    )
  }
  for (sm in samples) {
    resultList[[sm]][["Junction"]] = junctionsResults[[sm]]
  }
  rm(junctionsResults)
  gc()

  ## do Assessment of duplication rates from package dupRadar
  if (is.null(param$dupRadar) || param$dupRadar == TRUE) {
    dupRateFiles = getUpstreamQcFiles(input, "DupRate")
    if (
      length(dupRateFiles) == length(samples) &&
        hasSameUpstreamSettings(input, param)
    ) {
      ezLog("using the duplication rates computed by the aligner")
      dupRateResults <- lapply(dupRateFiles, ezRead.table, row.names = NULL)
    } else {
      dupRateResults <- ezMclapply(
        files,
        getDupRateFromBam,
        param,
        mc.preschedule = FALSE,
        mc.cores = min(length(files), ezThreads())
      )
    }
    for (sm in samples) {
      resultList[[sm]][["dupRate"]] = dupRateResults[[sm]]
    }
    rm(dupRateResults)
    gc()
  }

  if (input$hasColumn("StrandFile")) {
    ## checkExists must be off so that the file.exists() test below is reached
    strandFiles <- input$getFullPaths("StrandFile", checkExists = FALSE)
    for (sm in samples) {
      if (ezIsSpecified(strandFiles[sm]) && file.exists(strandFiles[sm])) {
        resultList[[sm]][["Strandness"]] <- ezReadRSeQCStrandness(strandFiles[
          sm
        ])
      }
    }
  }

  makeRmdReport(
    dataset = dataset,
    param = param,
    resultList = resultList,
    rmdFile = "RNABamStats.Rmd",
    reportTitle = "RNA BAM Stats",
    selfContained = TRUE
  )

  rm(resultList)
  gc()
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodRnaBamStats(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
EzAppRnaBamStats <-
  setRefClass(
    "EzAppRnaBamStats",
    contains = "EzApp",
    methods = list(
      ## RSeQC/Rsamtools/GenomicAlignments unconditional; Picard MarkDuplicates +
      ## dupRadar (uses Rsubread::featureCounts internally) gated on param$dupRadar.
      citation = function() {
        c(
          "Wang, L., Wang, S. & Li, W. RSeQC: quality control of RNA-seq experiments. Bioinformatics 28(16), 2184-2185 (2012). https://doi.org/10.1093/bioinformatics/bts356",
          "Lawrence, M. et al. Software for Computing and Annotating Genomic Ranges. PLoS Computational Biology 9(8), e1003118 (2013). https://doi.org/10.1371/journal.pcbi.1003118",
          "Morgan, M. & Pagès, H. Rsamtools: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import. R package version 2.28.0. https://doi.org/10.18129/B9.bioc.Rsamtools",
          "Picard Toolkit. Broad Institute. https://broadinstitute.github.io/picard/",
          "Sayols, S., Scherzinger, D. & Klein, H. dupRadar: a Bioconductor package for the assessment of PCR artifacts in RNA-Seq data. BMC Bioinformatics 17, 428 (2016). https://doi.org/10.1186/s12859-016-1276-2",
          "Liao, Y., Smyth, G.K. & Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30(7), 923-930 (2014). https://doi.org/10.1093/bioinformatics/btt656"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodRnaBamStats
        name <<- "EzAppRnaBamStats"
        appDefaults <<- rbind(
          posErrorRates = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "compute position specific error rates?"
          ),
          dupRadar = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "run dupradar"
          ),
          fragSizeMax = ezFrame(
            Type = "integer",
            DefaultValue = 500,
            Description = "maximum fragment size to plot in fragment size distribution"
          ),
          writeIgvSessionLink = ezFrame(
            Type = "logical",
            DefaultValue = "TRUE",
            Description = "should an IGV link be generated"
          ),
          ignoreDup = ezFrame(
            Type = "logical",
            DefaultValue = "NA",
            Description = "should marked duplicates be ignored?"
          ),
          skipCountQc = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "should we skip the count QC as part of the report"
          )
        )
      }
    )
  )


##' @describeIn computeBamStats Gets the files of an optional input column that were written by
##' an upstream app. Such a column can be missing altogether in datasets from an older app version
##' and the file can be empty if the upstream computation failed.
getUpstreamQcFiles <- function(input, columnName) {
  if (!input$hasColumn(columnName)) {
    return(character(0))
  }
  ## checkExists must be off, the column can name a file that is not there any more
  files <- input$getFullPaths(columnName, checkExists = FALSE)
  use <- file.exists(files)
  use[use] <- file.size(files[use]) > 0
  return(files[use])
}

##' @describeIn computeBamStats Checks whether the annotation and the library settings of the
##' upstream app are the ones requested here. If they differ, its results must not be reused.
hasSameUpstreamSettings <- function(input, param) {
  normalize <- function(x) {
    tolower(trimws(as.character(x)))
  }
  for (columnName in c("refFeatureFile", "strandMode", "paired")) {
    if (!input$hasColumn(columnName)) {
      return(FALSE)
    }
    if (
      !all(
        normalize(input$getColumn(columnName)) == normalize(param[[columnName]])
      )
    ) {
      ezLog("the upstream ", columnName, " differs; recomputing")
      return(FALSE)
    }
  }
  return(TRUE)
}

##' @describeIn computeBamStats Gets the error positions from the BAM file.
getPosErrorFromBam = function(bamFile, param) {
  require("bitops", warn.conflicts = WARN_CONFLICTS, quietly = !WARN_CONFLICTS)
  job = ezJobStart(paste("position error:", bamFile))
  seqLengths = ezBamSeqLengths(bamFile)
  if (ezIsSpecified(param$seqNames)) {
    seqLengths <- seqLengths[param$seqNames]
  }
  chromSel <- names(seqLengths)[which.max(seqLengths)]
  fai = fasta.index(param$ezRef["refFastaFile"])
  targetGenome = readDNAStringSet(fai[
    match(chromSel, sub(" .*", "", fai$desc)),
  ])[[1]] ## the sub clips potentially present description terms from the read id

  result = list()
  ezWriteElapsed(job, "read reference genome done")
  what = c("strand", "cigar", "seq", "rname", "pos")
  if (param$paired) {
    reads = ezScanBam(
      bamFile,
      seqname = chromSel,
      isFirstMateRead = TRUE,
      isSecondMateRead = FALSE,
      isUnmappedQuery = FALSE,
      isDuplicate = ezBamFlagSkip(param$ignoreDup),
      what = what
    )
    result[[paste(
      chromSel,
      "Position Stats of First Read"
    )]] = ezPosSpecErrorRate(reads, targetGenome)
    reads = ezScanBam(
      bamFile,
      seqname = chromSel,
      isFirstMateRead = FALSE,
      isSecondMateRead = TRUE,
      isUnmappedQuery = FALSE,
      isDuplicate = ezBamFlagSkip(param$ignoreDup),
      what = what
    )
    result[[paste(
      chromSel,
      "Position Stats of Second Read"
    )]] = ezPosSpecErrorRate(reads, targetGenome)
  } else {
    reads = ezScanBam(
      bamFile,
      seqname = chromSel,
      isUnmappedQuery = FALSE,
      isDuplicate = ezBamFlagSkip(param$ignoreDup),
      what = what
    )
    result[[paste(chromSel, "Position Stats")]] = ezPosSpecErrorRate(
      reads,
      targetGenome
    )
  }
  ezWriteElapsed(job, "Position Error Rate done")
  rm(reads)
  rm(targetGenome)
  gc()
  return(result)
}

##' @describeIn computeBamStats Calculates the specific error rates for \code{getPosErrorFromBam()}.
ezPosSpecErrorRate = function(bam, ReferenceGenome, nMaxReads = 100000) {
  require("Hmisc", warn.conflicts = WARN_CONFLICTS, quietly = !WARN_CONFLICTS)
  ## remove the reads containing the gaps, insertions, deletions
  hasGap = grepl("N|I|D", bam$cigar)
  readLength <- width(bam$seq)
  isOutOfRange = bam$pos + readLength - 1 > length(ReferenceGenome) |
    bam$pos < readLength ## this is very conservative; needed because there might be clipped bases in the beginning
  if (any(isOutOfRange)) {
    ezWrite("#reads out of range: ", sum(isOutOfRange))
    ezWrite("last pos: ", length(ReferenceGenome))
    idx = which(isOutOfRange)
    idx = idx[1:min(10, length(idx))]
    badAlignments = data.frame(
      pos = bam$pos[idx],
      cigar = bam$cigar[idx],
      start = bam$pos[idx],
      width = readLength[idx]
    )
    print(badAlignments)
  }
  indexKeep = which(!hasGap & !isOutOfRange)
  if (length(indexKeep) > nMaxReads) {
    indexKeep = sample(indexKeep, size = nMaxReads, replace = FALSE)
  }
  ezWrite(
    "#alignments: ",
    length(bam$cigar),
    " #valid alignments: ",
    sum(!hasGap & !isOutOfRange),
    " #used:",
    length(indexKeep)
  )
  if (length(indexKeep) == 0) {
    return(NULL)
  }
  for (nm in setdiff(names(bam), "tag")) {
    ## do the replacement in place in order to save memory
    bam[[nm]] = bam[[nm]][indexKeep]
  }
  ## treat the tag separately
  for (tagName in names(bam$tag)) {
    bam$tag[[tagName]] = bam$tag[[tagName]][indexKeep]
  }

  ## adjust the start POS according to H and/or S
  tempCigar = str_extract(bam$cigar, "^(\\d+H)?(\\d+S)?\\d+M")
  ## get the number of H at the beginning
  clipCigar = str_extract(tempCigar, "^\\d+H")
  noOfH = as.integer(sub("H", "", clipCigar))
  noOfH[is.na(noOfH)] = 0
  ## get the number of S at the beginning
  clipCigar = str_extract(tempCigar, "\\d+S")
  noOfS = as.integer(sub("S", "", clipCigar))
  noOfS[is.na(noOfS)] = 0
  nBeginClipped = noOfH + noOfS
  bam$pos = bam$pos - nBeginClipped

  ## add X to the begin and end of SEQ
  Xbegin = makeNstr("X", noOfH)
  tempCigar = str_extract(bam$cigar, "(\\d+S|\\d+H)$")
  clipCigar = str_extract(tempCigar, "\\d+H$")
  noOfH = as.integer(sub("H", "", clipCigar))
  noOfH[is.na(noOfH)] = 0
  Xend = makeNstr("X", noOfH)
  clipCigar = str_extract(tempCigar, "\\d+S")
  noOfS = as.integer(sub("S", "", clipCigar))
  noOfS[is.na(noOfS)] = 0
  nEndClipped = noOfH + noOfS
  bam$seq = paste0(Xbegin, bam$seq, Xend)

  seqChar = strsplit(bam$seq, "")
  readLength <- lengths(seqChar)
  ## build the reference views object
  maxLength = quantile(readLength, 0.95)
  if (maxLength < max(readLength)) {
    readLength[readLength > maxLength] = maxLength
    seqChar = mapply(
      function(x, l) {
        x[1:l]
      },
      seqChar,
      readLength
    )
  }
  ReferenceViews = Views(ReferenceGenome, start = bam$pos, width = readLength)
  referenceChar = strsplit(as.character(ReferenceViews), "")

  # assuming we have unique read length and set it to the maximal read length here.
  nEndTrimmed = maxLength - readLength
  trimmedMatrix = mapply(
    function(readLength, nEndTrimmed) {
      rep(c(FALSE, TRUE), c(readLength, nEndTrimmed))
    },
    readLength,
    nEndTrimmed,
    SIMPLIFY = FALSE
  )
  ## build a clippedMatrix to record the clipped character

  if (any(nEndClipped >= readLength)) {
    nEndClipped[which(nEndClipped > readLength)] <- readLength[which(
      nEndClipped > readLength
    )] -
      2
  }

  if (any((readLength - nBeginClipped - nEndClipped) <= 0)) {
    pos <- which((readLength - nBeginClipped - nEndClipped) <= 0)

    for (i in 1:length(pos)) {
      res <- readLength[pos[i]] - nBeginClipped[pos[i]] - nEndClipped[pos[i]]
      while (res < 1) {
        nBeginClipped[pos[i]] <- max(0, nBeginClipped[pos[i]] - 1)
        nEndClipped[pos[i]] <- max(0, nEndClipped[pos[i]] - 1)
        res <- readLength[pos[i]] - nBeginClipped[pos[i]] - nEndClipped[pos[i]]
      }
    }
  }

  nNormal = readLength - nBeginClipped - nEndClipped
  clippedMatrix = mapply(
    function(nBeginClipped, nNormal, nEndClipped, nEndTrimmed) {
      rep(
        c(TRUE, FALSE, TRUE, FALSE),
        c(nBeginClipped, nNormal, nEndClipped, nEndTrimmed)
      )
    },
    nBeginClipped,
    nNormal,
    nEndClipped,
    nEndTrimmed,
    SIMPLIFY = FALSE
  )

  matchMatrix = mapply("==", referenceChar, seqChar, SIMPLIFY = FALSE)
  indexNeg = which(bam$strand == "-")
  clippedMatrix[indexNeg] = lapply(clippedMatrix[indexNeg], rev)
  ###  trimmedMatrix[indexNeg] = lapply(trimmedMatrix[indexNeg], rev) ## trimmed matrix must not be reverted!!!
  matchMatrix[indexNeg] = lapply(matchMatrix[indexNeg], rev)
  lengthPadding = maxLength - readLength
  matchMatrix = sapply(matchMatrix, function(x) {
    x[1:maxLength]
  })
  clippedMatrix = matrix(unlist(clippedMatrix), ncol = length(clippedMatrix))
  clippedRate = rowMeans(clippedMatrix, na.rm = TRUE)
  trimmedMatrix = matrix(unlist(trimmedMatrix), ncol = length(trimmedMatrix))
  trimmedRate = rowMeans(trimmedMatrix, na.rm = TRUE)
  ## To distinguish error rate and clipped rate, remove the clipped from the mismatch, make it as match
  matchMatrix[clippedMatrix] = NA
  errorRate = 1 - rowMeans(matchMatrix, na.rm = TRUE)
  names(errorRate) = 1:nrow(matchMatrix)
  names(clippedRate) = 1:nrow(matchMatrix)
  names(trimmedRate) = 1:nrow(matchMatrix)
  return(list(
    trimmedRate = trimmedRate,
    clippedRate = clippedRate,
    errorRate = errorRate
  ))
}

##' @describeIn computeBamStats Gets the result statistics from the BAM file.
getStatsFromBam = function(
  param,
  bamFile,
  sm,
  gff = NULL,
  repeatsGff = NULL,
  nReads = NA
) {
  require("bitops", warn.conflicts = WARN_CONFLICTS, quietly = !WARN_CONFLICTS)
  seqLengths = ezBamSeqLengths(bamFile)

  ## determine the transcripts for which we compute the coverage
  use = rep(TRUE, nrow(gff))
  gff$transcript_type = ezGffAttributeField(
    gff$attributes,
    field = "transcript_type",
    attrsep = "; *",
    valuesep = " "
  )
  if (all(is.na(gff$transcript_type))) {
    gff$transcript_type = ezGffAttributeField(
      gff$attributes,
      field = "transcript_biotype",
      attrsep = "; *",
      valuesep = " "
    )
  }

  if (any(!is.na(gff$transcript_type))) {
    use = use & gff$transcript_type %in% c("protein_coding", "mRNA")
  }

  gff$tsl = ezGffAttributeField(
    gff$attributes,
    field = "transcript_support_level",
    attrsep = "; *",
    valuesep = " "
  )
  if (any(!is.na(gff$tsl))) {
    if (any(gff$tsl %in% "5")) {
      use = use & gff$tsl %in% "5"
    }
  }

  tCount = tapply(gff$transcript_id[use], gff$gene_id[use], function(x) {
    length(unique(x))
  })
  genesWithFewTranscripts <- names(tCount)[tCount <= 2]
  transcriptsForCov = unique(gff$transcript_id[
    gff$gene_id %in% genesWithFewTranscripts
  ])

  if (ezIsSpecified(param$seqNames)) {
    seqLengths = seqLengths[param$seqNames]
  }
  if (is.null(param$splitByChrom) || param$splitByChrom) {
    result = getStatsFromBamParallel(
      seqLengths,
      param,
      bamFile,
      sm,
      gff,
      repeatsGff,
      mc.cores = param$cores,
      nReads = nReads,
      transcriptsForCov = transcriptsForCov
    )
  } else {
    result = getStatsFromBamSingleChrom(
      NULL,
      param,
      bamFile,
      sm,
      gff,
      repeatsGff,
      transcriptsForCov = transcriptsForCov
    )
  }
  gc()
  ## TODO: this getBamMultiMatching should be moved to computeBamStats
  result$multiMatchInFileTable = getBamMultiMatching(param, bamFile, nReads)

  transcriptCov = result$transcriptCov
  transcriptCovRleList <- RleList(transcriptCov)
  # transcriptLengthCov = sapply(transcriptCov, function(x){sum(x>0)}) slow!
  transcriptLengthCov = sum(transcriptCovRleList > 0)
  # transcriptLengthTotal = sapply(transcriptCov, length) slow!
  transcriptLengthTotal <- elementNROWS(transcriptCovRleList)
  percentCovered = transcriptLengthCov / transcriptLengthTotal * 100
  result$TranscriptsCovered = table(ezCut(
    percentCovered,
    breaks = c(0.5, 10, 90, 99.5),
    labels = c("not covered", "<10%", "10 - 90%", ">90%", "fully covered")
  ))

  ## Do the genebody_coverage
  #sampledTranscriptCov = sapply(transcriptCov, ## RleList will be slow. use list
  #                              function(x){as.integer(x[round(seq(1, length(x), length.out=101))])})
  sampledTranscriptCov <- ezMclapply(
    transcriptCov, ## RleList will be slow. use list
    function(x) {
      as.integer(x[round(seq(1, length(x), length.out = 101))])
    },
    mc.preschedule = TRUE,
    mc.cores = param$cores
  )
  sampledTranscriptCov <- do.call(cbind, sampledTranscriptCov)

  trUse = colSums(sampledTranscriptCov) > 0
  sampledTranscriptCov = sampledTranscriptCov[, trUse, drop = FALSE]
  trLength = transcriptLengthTotal[trUse]
  lengthClasses = ezCut(
    trLength,
    breaks = c(599, 1200, 2400),
    labels = c(
      "less than 600nt",
      "600 to 1199nt",
      "1200 to 2400nt",
      "above 2400nt"
    )
  )
  genebody_coverage = list()
  for (lc in levels(lengthClasses)) {
    isInLc = lengthClasses == lc
    if (sum(isInLc) > 40) {
      ltc = sampledTranscriptCov[, isInLc, drop = FALSE]
      avgCov = colMeans(ltc)
      #relativeCov = ezScaleColumns(ltc, 1/colSums(ltc)) ## normalize so that every transcripts adds the same weight
      avgCovQuant = unique(quantile(avgCov, c(0.25, 0.75)))
      if (length(avgCovQuant) == 2) {
        covClasses = ezCut(
          avgCov,
          breaks = avgCovQuant,
          labels = c("low expressed", "medium expressed", "high expressed")
        )
        genebody_coverage[[lc]] = list()
        for (cc in levels(covClasses)) {
          genebody_coverage[[lc]][[cc]] = rowMeans(ltc[,
            covClasses == cc,
            drop = FALSE
          ])
        }
      }
    }
  }
  result$genebody_coverage = genebody_coverage
  result$transcriptCov = NULL
  gc()
  return(result)
}

##' @describeIn computeBamStats Gets parallel by chromosome statistics for \code{getStatsFromBam()} if the logical \code{param$splitByChrom} is true.
getStatsFromBamParallel = function(
  seqLengths,
  param,
  bamFile,
  sm,
  gff = NULL,
  repeatsGff = NULL,
  mc.cores = ezThreads(),
  nReads = NA,
  transcriptsForCov = NULL
) {
  if (!is.na(nReads)) {
    ## heuristic: reduce the number of cores so that we have at least 0.25GB RAM per chromosome per Million Reads in the total bam file
    #reduce the number of threads in case of
    maxCores = ceiling(param$ram / (nReads / 1e6) * 4)
    if (maxCores < mc.cores) {
      ezLog("too many reads --reducing the number of cores: ", maxCores)
      mc.cores = maxCores
    }
  }
  seqNames <- names(sort(seqLengths, decreasing = TRUE)) ## sorting so that longest job starts first
  ## use heuristics to select only chromosomes: NCBI chromosomes start with NC; others have less than 6
  if (sum(grepl("^NC_", seqNames)) > 1) {
    ## NCBI
    seqNames <- grep("^NC_", seqNames, value = TRUE)
  } else {
    seqNames <- seqNames[nchar(seqNames) <= 6] ## remove non-chromosome sequences that usually have long names
  }
  if (length(seqNames) == 0) {
    seqNames <- names(sort(seqLengths, decreasing = TRUE)) %>% head(100)
  }
  names(seqNames) <- seqNames ## set names for lapply
  chromResults = ezMclapply(
    seqNames,
    getStatsFromBamSingleChrom,
    param,
    bamFile,
    sm,
    gff,
    repeatsGff,
    transcriptsForCov,
    mc.preschedule = FALSE,
    mc.cores = mc.cores
  )
  if (param$saveImage) {
    save(chromResults, file = paste0(sm, "-chromResults.RData"))
  }
  gc()
  result = list()

  # merge the fragSizeHist
  idx = which(sapply(chromResults, function(x) {
    !is.null(x$fragSizeHist)
  }))
  if (length(idx) > 0) {
    for (i in idx) {
      if (i == idx[1]) {
        fsh = chromResults[[i]]$fragSizeHist
        counts = fsh$counts
      } else {
        counts = counts + chromResults[[i]]$fragSizeHist$counts
      }
    }
    fsh$counts = counts
    fsh$density = counts / sum(counts)
    result$fragSizeHist = fsh
  }
  # merge the multiMatchTargetTypeCounts
  temp = data.frame(count = integer(0), width = integer(0))
  for (i in 1:length(chromResults)) {
    newResult = chromResults[[i]]$multiMatchTargetTypeCounts
    ## extend the temp data frame if needed
    additionalRows = setdiff(rownames(newResult), rownames(temp))
    if (length(additionalRows) > 0) {
      temp[additionalRows, ] = 0
    }
    temp[rownames(newResult), ] = temp[rownames(newResult), ] + newResult
  }
  if (any(is.na(temp))) {
    ezLog("na counts: ", sum(is.na(temp)))
    temp[is.na(temp)] = 0
  }

  tempNamesOrdered = intersect(
    c(setdiff(rownames(temp), seqNames), seqNames),
    rownames(temp)
  )
  result$multiMatchTargetTypeCounts = temp[tempNamesOrdered, , drop = FALSE]
  rm(temp)
  gc()
  result$seqLengths = seqLengths

  result$transcriptCov = unlist(lapply(chromResults, function(cr) {
    cr$transcriptCov
  }))

  return(result)
}

##' @describeIn computeBamStats Gets the statistics of a single chromosome for \code{getStatsFromBam()}.
getStatsFromBamSingleChrom = function(
  chrom,
  param,
  bamFile,
  sm,
  gff = NULL,
  repeatsGff = NULL,
  transcriptsForCov = NULL
) {
  require("bitops", warn.conflicts = WARN_CONFLICTS, quietly = !WARN_CONFLICTS)
  ezLog("Processing chr ", ifelse(is.null(chrom), "all", chrom))

  result = list()
  seqLengths = ezBamSeqLengths(bamFile)
  if (param$paired) {
    reads = ezReadPairedAlignments(
      bamFile,
      seqname = chrom,
      keepUnpaired = param$keepUnpaired,
      minMapQuality = param$minMapQuality,
      keepMultiHits = param$keepMultiHits
    )
  } else {
    reads = ezReadGappedAlignments(
      bamFile,
      seqname = chrom,
      minMapQuality = param$minMapQuality,
      keepMultiHits = param$keepMultiHits
    )
  }
  if (isError(reads)) {
    return(reads)
  }
  gc()

  if (param$paired && length(reads) > 0) {
    pairedNames = ezScanBam(
      bamFile,
      seqname = chrom,
      isFirstMateRead = TRUE,
      isSecondMateRead = FALSE,
      isProperPair = ezBamFlagKeepOnly(param$keepProperPairsOnly),
      isUnmappedQuery = FALSE,
      isDuplicate = ezBamFlagSkip(param$ignoreDup),
      what = "qname"
    )$qname
    use = names(reads) %in% pairedNames
    result$fragSizeHist = intHist(
      width(reads)[use],
      range = c(-0.5, param$fragSizeMax + 0.5),
      plot = FALSE
    )
    rm(pairedNames)
    gc()
  }

  result$multiMatchTargetTypeCounts = getTargetTypeCounts(
    param,
    gff,
    reads,
    seqid = chrom,
    repeatsGff
  )
  gc()
  ## Do transcripts covered
  result$transcriptCov = getTranscriptCoverage(
    chrom,
    gff[gff$transcript_id %in% transcriptsForCov, ],
    reads,
    strandMode = param$strandMode
  )
  gc()
  rm(reads)
  #rm(isMultiHit)
  gc()
  return(result)
}

##' @describeIn computeBamStats Gets the counts of the target types for \code{getStatsFromBam()}.
getTargetTypeCounts = function(
  param,
  gff,
  rr,
  seqid = NULL,
  repeatsGff = NULL
) {
  require(data.table)
  if (class(rr) == "GRangesList") {
    #sn = unlist(sn, use.names=FALSE)[sn@partitioning@end]
    stop("GRangesList not supported")
  }
  seqNames = names(seqlengths(rr))
  if (!is.null(seqid)) {
    stopifnot(length(seqid) == 1)
    seqNames = intersect(seqNames, seqid)
    readRefIsValid = as.character(seqnames(rr)) %in% seqNames
    if (!all(readRefIsValid)) {
      rr = rr[readRefIsValid]
    }
  }
  #effWidth = sum(as.numeric(seqlengths(rr))) * ifelse(param$isStranded, 2, 1)
  effWidth = (seqlengths(rr) * ifelse(param$strandMode == "both", 1, 2))[
    seqNames
  ]
  result = data.frame(count = 0, width = effWidth, row.names = seqNames)
  #readCounts = table(as.character(seqnames(rr)))
  readCounts <- table(seqnames(rr)) # It also works with different order

  result[seqNames, "count"] = readCounts[seqNames]
  result$count[is.na(result$count)] = 0 ## if a chromosome has no reads the value would be na

  #hasAnyHit = rep(FALSE, length(rr))
  hasAnyHit <- logical(length(rr))
  repeatsRanges = NULL
  gffRanges = NULL

  ## the reads in the repeatsGff
  if (!is.null(repeatsGff)) {
    repeatsGff = repeatsGff[repeatsGff$seqid %in% seqNames, ]
    if (nrow(repeatsGff) > 0) {
      repeatsGff$strand = fixStrand(repeatsGff$strand, param$strandMode)
      repeatsRanges = gffToRanges(repeatsGff)
      classFam = ezGffAttributeField(repeatsGff$attributes, field = "repClass")
      repFamily = ezGffAttributeField(
        repeatsGff$attributes,
        field = "repFamily"
      )
      use = repFamily != classFam
      classFam[use] = paste(classFam[use], repFamily[use], sep = "--")
      for (type in unique(classFam)) {
        use = classFam == type
        targetRanges = gffToRanges(repeatsGff[use, ])
        hitsTarget = overlapsAny(rr, targetRanges, minoverlap = 10)
        result[type, ] = c(
          sum(hitsTarget),
          sum(width(IRanges::reduce(targetRanges)))
        )
        hasAnyHit = hasAnyHit | hitsTarget
      }
    }
  }

  if (!is.null(gff)) {
    gff = gff[gff$seqid %in% seqNames, ]
    if (nrow(gff) > 0) {
      gff$strand = fixStrand(gff$strand, param$strandMode)
      gffRanges = gffToRanges(gff)
      ensemblTypes = getEnsemblTypes(gff)
      if (!is.null(ensemblTypes)) {
        # for (type in unique(ensemblTypes)){
        #   targetRanges = gffRanges[ensemblTypes == type]
        #   hitsTarget = overlapsAny(rr, targetRanges, minoverlap=10)
        #   result[type, ] = c(sum(hitsTarget), sum(width(IRanges::reduce(targetRanges))))
        #   hasAnyHit = hasAnyHit | hitsTarget
        # }
        ## The following code is much faster than the loop above.
        ## The loop: 537.492 seconds;
        ## New implementation: 184 seconds
        hits <- findOverlaps(rr, gffRanges, minoverlap = 10)
        hitsByType <- data.table(
          queryHits = queryHits(hits),
          ensemblTypes = ensemblTypes[subjectHits(hits)]
        )
        hitsByType <- unique(hitsByType)
        countsByType <- hitsByType[, .N, by = ensemblTypes]
        if (nrow(countsByType) == 0L) {
          ## When there is no hit for all reads
          countsByType <- data.table(ensemblTypes = unique(ensemblTypes), N = 0)
        }
        ## some ensemblType has not hits and not included in findOverlaps.
        ## Add them with N=0
        missingTypes <- setdiff(unique(ensemblTypes), countsByType$ensemblTypes)
        if (length(missingTypes) != 0L) {
          countsByType <- rbind(
            countsByType,
            data.table(ensemblTypes = missingTypes, N = 0)
          )
        }

        widthByType <- sum(width(IRanges::reduce(GenomicRanges::split(
          gffRanges,
          ensemblTypes
        ))))

        hasAnyHit[hitsByType$queryHits] <- TRUE
        result[countsByType$ensemblTypes, ] <- cbind(
          countsByType$N,
          widthByType[countsByType$ensemblTypes]
        )

        isMsg = ensemblTypes == "protein_coding" & gff$type == "exon"
        msgRanges = gffGroupToRanges(
          gff[isMsg, ],
          gff$transcript_id[isMsg],
          skipTransSpliced = TRUE
        )
        targetExonRanges = gffRanges[isMsg]
      } else {
        rootTypes = setdiff(gff$type, c("intron", "exon"))
        for (type in rootTypes) {
          use = gff$type == type
          targetRanges = gffRanges[gff$type == type]
          hitsTarget = overlapsAny(rr, targetRanges, minoverlap = 10)
          result[type, ] = c(
            sum(hitsTarget),
            sum(width(IRanges::reduce(targetRanges)))
          )
          hasAnyHit = hasAnyHit | hitsTarget
        }
        isExon = gff$type == "exon"
        msgRanges = gffGroupToRanges(
          gff[isExon, ],
          gff$transcript_id[isExon],
          skipTransSpliced = TRUE
        )
        targetExonRanges = gffRanges[isExon]
      }
      ## Add seqlengths for msgRanges because of potential out-of-bound by flank
      ## TODO: now we use seqlengths from rr, which is not ideal.
      ## gff as data.frame has no seqlengths information.
      seqlengths(msgRanges) <- seqlengths(rr)[names(seqlengths(msgRanges))]

      ## check additionally for intron/exon/prom
      hitsTranscript = overlapsAny(rr, msgRanges, minoverlap = 10)
      hasAnyHit = hasAnyHit | hitsTranscript
      mRnaWidth = sum(width(IRanges::reduce(msgRanges)))
      hitsTargetExons = overlapsAny(
        rr[hitsTranscript],
        targetExonRanges,
        minoverlap = 10
      )
      result["mRNA Exons", ] = c(
        sum(hitsTargetExons),
        sum(width(IRanges::reduce(targetExonRanges)))
      )
      result["mRNA Introns", ] = c(
        sum(!hitsTargetExons),
        mRnaWidth - result["mRNA Exons", "width"]
      )
      ## suppressWarnings for out-of-bound ranges.
      promRanges = trim(suppressWarnings(flank(msgRanges, 2000)))
      hitsTargetProms = overlapsAny(rr, promRanges, minoverlap = 10)
      result["mRNA Promoter 2kb", ] = c(
        sum(hitsTargetProms),
        sum(width(IRanges::reduce(promRanges)))
      )
      hasAnyHit = hasAnyHit | hitsTargetProms
      downRanges = trim(suppressWarnings(flank(msgRanges, 2000, start = FALSE)))
      hitsTargetDown = overlapsAny(rr, downRanges, minoverlap = 10)
      result["mRNA Downstream 2kb", ] = c(
        sum(hitsTargetDown),
        sum(width(IRanges::reduce(promRanges)))
      )
      hasAnyHit = hasAnyHit | hitsTargetDown
      gffRanges = c(gffRanges, promRanges, downRanges)
    }
  }
  allRanges = GRanges()
  if (!is.null(gffRanges)) {
    allRanges = c(allRanges, gffRanges)
  }
  if (!is.null(repeatsGff)) {
    allRanges = c(allRanges, repeatsRanges)
  }
  if (length(allRanges) > 0) {
    annotatedWidth = sum(width(IRanges::reduce(allRanges)))
  } else {
    annotatedWidth = 0
  }
  #result["unannotated", ] = c(sum(!hasAnyHit), result[seqNames, "width"] - annotatedWidth)
  result["unannotated", ] = c(
    sum(!hasAnyHit),
    sum(result[seqNames, "width"]) - annotatedWidth
  )
  result["total", ] = colSums(result[seqNames, ], na.rm = TRUE)
  result = result[c("total", setdiff(rownames(result), "total")), ]
  return(result)
}

##' @describeIn computeBamStats Gets the annotated splice sites from the gff. The positions
##' follow the convention of STAR's SJ.out.tab, i.e. the first and the last base of the intron, 1-based.
getAnnotatedSpliceSites <- function(gff) {
  require(data.table)
  use <- gff$type == "exon" & !is.na(gff$transcript_id)
  exons <- data.table(
    seqid = gff$seqid[use],
    start = as.integer(gff$start[use]),
    end = as.integer(gff$end[use]),
    transcript_id = gff$transcript_id[use]
  )
  setorder(exons, transcript_id, start)
  exons[, intronStart := end + 1L]
  ## shift must be qualified, IRanges masks the data.table version
  exons[,
    intronEnd := data.table::shift(start, type = "lead") - 1L,
    by = transcript_id
  ]
  introns <- exons[!is.na(intronEnd) & intronEnd >= intronStart]
  return(list(
    starts = unique(paste(introns$seqid, introns$intronStart)),
    ends = unique(paste(introns$seqid, introns$intronEnd))
  ))
}

##' @describeIn computeBamStats Gets the junction results from STAR's SJ.out.tab file.
##' This is the cheap alternative to \code{getJunctionPlotsFromBam()} which has to read the
##' whole bam file again.
getJunctionStatsFromSJ = function(sjFile, spliceSites, minIntronSize = 50) {
  require(data.table)
  sj <- fread(
    sjFile,
    header = FALSE,
    col.names = c(
      "seqid",
      "start",
      "end",
      "strand",
      "motif",
      "annotatedInSjdb",
      "nUnique",
      "nMulti",
      "maxOverhang"
    )
  )
  ## the annotatedInSjdb column is deliberately not used: with --twopassMode the junctions
  ## found in the first pass are inserted into the database and reported as annotated
  sj <- sj[end - start + 1L >= minIntronSize & nUnique + nMulti > 0]
  if (nrow(sj) == 0) {
    return(list())
  }
  readCount <- sj$nUnique + sj$nMulti
  ## same definition as in RSeQC junction_annotation.py: the two splice sites are judged
  ## independently, so a novel combination of two known splice sites counts as annotated
  startKnown <- paste(sj$seqid, sj$start) %in% spliceSites$starts
  endKnown <- paste(sj$seqid, sj$end) %in% spliceSites$ends
  annotation <- ifelse(
    startKnown & endKnown,
    "annotated",
    ifelse(!startKnown & !endKnown, "complete_novel", "partial_novel")
  )

  result = list()
  result[["splice_junction"]] = table(annotation) / length(annotation) * 100
  result[["splice_events"]] = as.table(tapply(readCount, annotation, sum)) /
    sum(readCount) *
    100

  ## expected number of distinct junctions when a fraction of the reads is sampled;
  ## this closed form replaces the repeated random subsampling of the individual reads
  quantiles = seq(0.05, 1, by = 0.05)
  saturation <- function(counts) {
    values <- sapply(quantiles, function(q) {
      sum(1 - (1 - q)^counts)
    })
    names(values) <- as.character(quantiles)
    return(values)
  }
  isNovel <- annotation != "annotated"
  result[["junctionSaturation"]] = list(
    "all junctions" = saturation(readCount),
    "known junctions" = saturation(readCount[!isNovel]),
    "novel junctions" = saturation(readCount[isNovel])
  )
  return(result)
}

##' @describeIn computeBamStats Gets the junction results from the BAM file.
getJunctionPlotsFromBam = function(bamFile, param) {
  pngFiles = list()
  ## do the junction annotation
  outputJunction = paste0("junction-", Sys.getpid())
  bed = getReferenceFeaturesBed(param)
  stopifnot(!is.null(bed))
  # junction_annotation.py is in RSeQC package available through Dev/Python2
  cmd = paste(
    "junction_annotation.py",
    "--mapq=1",
    "-i",
    bamFile,
    "-o",
    outputJunction,
    "-r",
    bed
  )
  res = ezSystem(cmd, stopOnFailure = FALSE)
  junctionFile = paste0(outputJunction, ".junction.xls")
  if (res == 0 && length(readLines(junctionFile)) > 1) {
    juncsTable = read.table(junctionFile, header = TRUE)
    junctions = table(juncsTable$annotation) /
      length(juncsTable$annotation) *
      100
    foo = rep(juncsTable$annotation, juncsTable$read_count)
    events = table(foo) / length(foo) * 100
    pngFiles[["splice_events"]] = events
    pngFiles[["splice_junction"]] = junctions

    ## do the junction_saturation
    id = paste(
      juncsTable[["chrom"]],
      juncsTable[["intron_st.0.based."]],
      juncsTable[["intron_end.1.based."]]
    )
    juncReads = rep(id, juncsTable[["read_count"]])
    juncTypes = rep(juncsTable[["annotation"]], juncsTable[["read_count"]])
    juncTypeSet = c("annotated", "complete_novel", "partial_novel")
    nSim = 10
    quantiles = seq(0.05, 1, by = 0.05)
    juncCounts = array(
      0,
      dim = c(length(juncTypeSet), length(quantiles), nSim),
      dimnames = list(
        types = juncTypeSet,
        quantiles = as.character(quantiles),
        sim = 1:nSim
      )
    )
    for (n in 1:nSim) {
      idx = sample(1:length(juncReads), replace = FALSE)
      for (q in 1:length(quantiles)) {
        idxUse = idx[1:round(quantiles[q] * length(idx))]
        cts = tapply(juncReads[idxUse], juncTypes[idxUse], function(x) {
          length(unique(x))
        })
        juncCounts[names(cts), q, n] = as.vector(cts)
      }
    }
    ## TODO why can they be NA????
    juncCounts[is.na(juncCounts)] = 0
    juncCountMeans = apply(juncCounts, c(1, 2), mean)
    juncCountMeans[is.na(juncCountMeans)] = 0
    junctionSaturations = list()
    junctionSaturations[["all junctions"]] = colSums(juncCountMeans[
      c("annotated", "complete_novel", "partial_novel"),
    ])
    junctionSaturations[["known junctions"]] = colSums(juncCountMeans[
      c("annotated"),
      ,
      drop = FALSE
    ])
    junctionSaturations[["novel junctions"]] = colSums(juncCountMeans[
      c("complete_novel", "partial_novel"),
      ,
      drop = FALSE
    ])
    pngFiles[["junctionSaturation"]] = junctionSaturations
  }

  ## do the cleaning
  file.remove(list.files(path = ".", pattern = paste0(outputJunction, ".+")))
  return(pngFiles)
}


### classify reads according to where they match
### - assign reads to transcripts
### - split transcript reads into intronic/exonic
### - use remaining reads and assign to other root features
### - use remaining reads to assign to promoters
### - use remaining reads to downstream
### - use remaining reads to intergenic

## behaviour for strand-preserving library prep:
## - genome target size is 2 * genome size
## - match with ignoreStrand=false

## behaviour for non-strand-preserving library prep:
## - genome target size is 1 * genome size

# getTypePercentTable = function(resultList, name, minPercentage=1){
#   tbl = data.frame(row.names=rownames(resultList[[1]][[name]]))
#   for (sm in names(resultList)){
#     counts = resultList[[sm]][[name]]
#     tbl[sm] = round(counts[ rownames(tbl), "count"] / counts["total", "count"] * 100, digits=3)
#   }
#   useRow = apply(tbl, 1, max) > minPercentage & rownames(tbl) != "total"
#   return(tbl[useRow , ])
# }

getDupRateFromBam <- function(
  bamFile,
  param = NULL,
  gtfFn,
  stranded = c("both", "sense", "antisense"),
  paired = FALSE,
  threads = 1,
  markedBam = bamHasMarkedDuplicates(bamFile),
  ram = NULL
) {
  if (!is.null(param)) {
    gtfFn <- param$ezRef@refFeatureFile
    stranded <- param$strandMode
    paired <- param$paired
    if (missing(threads)) {
      threads <- ceiling(param$cores / 2)
    }
    if (is.null(ram)) {
      ram <- floor(param$ram / param$cores) - 1
    }
  }
  require(dupRadar)

  if (markedBam) {
    ## the duplicates are already flagged upstream, no need to run picard again
    bamDuprmFn <- bamFile
  } else {
    ## Mark the duplicates in bamFile
    inputBam <- paste(Sys.getpid(), basename(bamFile), sep = "-")
    ### The bamFile may not be writable.
    file.symlink(from = bamFile, to = inputBam)

    ## intermediate files
    # picardMetricsFn <- gsub("\\.bam$", "_picard_metrics.txt", inputBam) # picard
    bamDuprmFn <- gsub("\\.bam$", "_duprm.bam", inputBam)
    # bamutilLogFn <- paste0(bamDuprmFn, ".log") # bamutil
    on.exit(file.remove(c(inputBam, bamDuprmFn, paste0(bamDuprmFn, ".bai")))) #, picardMetricsFn, bamutilLogFn)))
    dupBam(inBam = inputBam, outBam = bamDuprmFn, operation = "mark", ram = ram)
  }
  ## Duplication rate analysis
  dm <- analyzeDuprates(
    bam = bamDuprmFn,
    gtf = gtfFn,
    stranded = switch(
      stranded,
      "both" = 0,
      "sense" = 1,
      "antisense" = 2,
      stop("unsupported strand mode: ", stranded)
    ),
    paired = paired,
    threads = threads
  )
  return(dm)
}

ezReadRSeQCStrandness <- function(file) {
  lines <- readLines(file)
  lines <- lines[lines != ""]
  if (length(lines) == 0) {
    return(NULL)
  }

  ## Extract fractions
  failedLine <- grep(
    "Fraction of reads failed to determine:",
    lines,
    value = TRUE
  )
  if (length(failedLine) > 0) {
    failed <- as.numeric(sub(
      "Fraction of reads failed to determine: ",
      "",
      failedLine
    ))
  } else {
    failed <- NA
  }

  explainedLines <- grep("Fraction of reads explained by", lines, value = TRUE)
  if (length(explainedLines) > 0) {
    explainedKeys <- sub(
      "Fraction of reads explained by \"(.*)\": .*",
      "\\1",
      explainedLines
    )
    explainedVals <- as.numeric(sub(
      "Fraction of reads explained by \".*\": ",
      "",
      explainedLines
    ))
    names(explainedVals) <- explainedKeys
  } else {
    explainedVals <- NULL
  }

  return(c(failed = failed, explainedVals))
}
