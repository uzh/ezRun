# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

### Filter ATAC-seq bam file; process by sample;
#### 1. remove duplicates
#### 2. chrM removal
#### 3. low mapping quality filtering (< 10)
#### 4. shift the 5' start
atacBamProcess <- function(input = NA, output = NA, param = NA) {
  require(GenomicAlignments)
  require(Rsamtools)
  require(rtracklayer)

  ## if output is not an EzDataset, set it!
  if (!is(output, "EzDataset")) {
    output <- input$copy()
    output$setColumn(
      "BAM",
      paste0(
        getwd(),
        "/",
        input$getNames(),
        "_processed.bam"
      )
    )
    output$setColumn(
      "BAI",
      paste0(
        getwd(),
        "/",
        input$getNames(),
        "_processed.bam.bai"
      )
    )
    output$dataRoot <- NULL
  }

  bamFile <- input$getFullPaths("BAM")

  ## For now, we only have paired-end ATAC-seq
  stopifnot(param$paired)

  ## remove the duplicates in bam
  noDupBam <- tempfile(pattern = "nodup_", tmpdir = ".", fileext = ".bam")
  ezLog("Remove duplicates...")
  if (param$removeDuplicates) {
    removeDuplicatesFromBam(
      inBam = bamFile,
      outBam = noDupBam,
      ram = param$ram,
      cores = param$cores
    )
  } else {
    outBam = noDupBam
    file.copy(from = bamFile, to = outBam, overwrite = TRUE)
    Rsamtools::indexBam(outBam)
  }
  noDupNoMTNOLowBam <- basename(output$getColumn("BAM"))
  filteroutBam(
    inBam = noDupBam,
    outBam = noDupNoMTNOLowBam,
    cores = param$cores,
    chrs = c("M", "MT", "chrM", "chrMT"),
    mapQ = 10
  )
  file.remove(c(noDupBam, paste0(noDupBam, ".bai")))

  if (param$shiftATAC) {
    cmd <- paste(
      'alignmentSieve --bam',
      noDupNoMTNOLowBam,
      '--outFile',
      paste0(noDupNoMTNOLowBam, '_shifted'),
      '--ATACshift',
      '-p',
      param$cores,
      '--smartLabels',
      '--filterMetrics shift.log'
    )
    ezSystem(cmd)
    cmd <- paste('mv', paste0(noDupNoMTNOLowBam, '_shifted'), noDupNoMTNOLowBam)
    ezSystem(cmd)
    file.remove(c(paste0(noDupNoMTNOLowBam, ".bai")))
    sortBam(
      noDupNoMTNOLowBam,
      paste0(noDupNoMTNOLowBam, '_sorted'),
      maxMemory = 1024 * param$ram
    )
    cmd <- paste(
      'mv',
      paste0(noDupNoMTNOLowBam, '_sorted.bam'),
      noDupNoMTNOLowBam
    )
    ezSystem(cmd)
    indexBam(noDupNoMTNOLowBam)
  }

  if (FALSE) {
    require(ATACseqQC)
    what <- c("qname", "flag", "mapq", "isize", "seq", "qual", "mrnm")
    tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
    flag <- scanBamFlag(
      isSecondaryAlignment = FALSE,
      isUnmappedQuery = FALSE,
      isNotPassingQualityControls = FALSE
    )
    #scanBamParam <- ScanBamParam(flag = flag, tag = tags, what = what, which = GRanges(seqname, IRanges(afrom, ato)))
    scanBamParam <- ScanBamParam(flag = flag, tag = tags, what = what)
    ezLog("Reading nodup noMT bam file...")
    reads <- readGAlignmentPairs(file = noDupNoMTNOLowBam, param = scanBamParam)
    file.remove(c(noDupNoMTNOLowBam, paste0(noDupNoMTNOLowBam, ".bai")))
    ## shiftAlignments
    ezLog("Shifting 5' start...")
    firstReads <- ATACseqQC:::shiftReads(
      first(reads),
      positive = 4L,
      negative = 5L
    )
    lastReads <- ATACseqQC:::shiftReads(
      last(reads),
      positive = 4L,
      negative = 5L
    )
    rm(reads)
    ezLog("Exporting bam file...")
    export(
      GAlignmentPairs(first = firstReads, last = lastReads),
      noDupNoMTNOLowBam
    )
    rm(firstReads)
    rm(lastReads)
    gc()
  }

  return(output)
}

### Make or remove duplicated in bam file
dupBam <- function(
  inBam,
  outBam,
  operation = c("mark", "remove"),
  ram = 20,
  dupDistance = 100
) {
  operation <- match.arg(operation)
  javaCall = paste0(
    "/usr/local/ngseq/packages/Dev/jdk/21/bin/java",
    " -Djava.io.tmpdir=. -Xmx",
    round(0.8 * ram),
    "g"
  )
  metricsFile <- sub('.bam$', '_metrics.txt', outBam)
  cmd <- paste(
    javaCall,
    " -jar ",
    Sys.getenv("Picard_jar"),
    "MarkDuplicates",
    paste0("I=", inBam),
    paste0("O=", outBam),
    paste0("M=", metricsFile),
    paste0("OPTICAL_DUPLICATE_PIXEL_DISTANCE=", dupDistance),
    paste0(
      "REMOVE_DUPLICATES=",
      if_else(operation == "mark", "false", "true")
    ),
    "> /dev/null 2>&1"
  )
  ezSystem(cmd)
  Rsamtools::indexBam(outBam)
  invisible(outBam)
}

##' @title Check whether a BAM already has its duplicates marked
##' @description
##' Inspect the BAM header for a program record left by a duplicate-marking
##' step (Picard MarkDuplicates or samtools markdup). When present, duplicates
##' are already flagged in the file, so a costly re-run of MarkDuplicates can be
##' skipped (e.g. Bowtie2/BWA now mark duplicates during alignment).
##' @param bamFile Path to a BAM file.
##' @return \code{TRUE} when the header shows duplicates were already marked.
##' @template roxygen-template
bamHasMarkedDuplicates <- function(bamFile) {
  if (length(bamFile) != 1 || is.na(bamFile) || !file.exists(bamFile)) {
    return(FALSE)
  }
  header <- tryCatch(
    system2("samtools", c("view", "-H", shQuote(bamFile)),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  )
  pgLines <- grep("^@PG", header, value = TRUE)
  any(grepl("MarkDuplicates|markdup", pgLines, ignore.case = TRUE))
}

##' @title Remove duplicate reads, skipping MarkDuplicates when already done
##' @description
##' Produce a duplicate-free BAM. When the input BAM already has duplicates
##' marked (see \code{\link{bamHasMarkedDuplicates}}), the flagged reads are
##' dropped directly with \code{samtools} (fast, no re-detection); otherwise
##' Picard \code{MarkDuplicates} is run with \code{REMOVE_DUPLICATES=true} via
##' \code{\link{dupBam}}.
##' @param inBam Full path to the input BAM file.
##' @param outBam Output BAM path.
##' @param ram RAM in GB passed to \code{dupBam} when MarkDuplicates is run.
##' @param cores Threads used by \code{samtools} in the fast path.
##' @return Invisibly, the output BAM path.
##' @template roxygen-template
removeDuplicatesFromBam <- function(inBam, outBam, ram = 20, cores = 1) {
  if (bamHasMarkedDuplicates(inBam)) {
    ezLog(paste0(
      "Duplicates are already marked in ", basename(inBam),
      "; dropping flagged reads with samtools (skipping MarkDuplicates)."
    ))
    ## -F 1024 removes reads with the 0x400 (duplicate) flag set.
    ezSystem(paste(
      "samtools view -b", "-@", cores, "-F 1024", "-o", outBam, inBam
    ))
    Rsamtools::indexBam(outBam)
  } else {
    ezLog(paste0(
      "Duplicates are not marked in ", basename(inBam),
      "; running Picard MarkDuplicates to remove them."
    ))
    dupBam(inBam = inBam, outBam = outBam, operation = "remove", ram = ram)
  }
  invisible(outBam)
}

### Filter bam by removing chrs, low mapQ
filteroutBam <- function(
  inBam,
  outBam,
  cores = ezThreads(),
  chrs = NULL,
  mapQ = NULL
) {
  require(Rsamtools)
  cmd <- paste("samtools view -b", "-@", cores, "-o", outBam)

  if (!is.null(mapQ)) {
    cmd <- paste(cmd, "-q", mapQ)
  }

  cmd <- paste(cmd, inBam)

  if (!is.null(chrs)) {
    allChrs <- ezBamSeqNames(inBam)
    keepChrs <- setdiff(allChrs, chrs)
    cmd <- paste(cmd, paste(keepChrs, collapse = " "))
  }
  ezSystem(cmd)

  indexBam(outBam)

  invisible(outBam)
}

### Convert bam to bigwig file
bam2bw <- function(
  file,
  destination = sub("\\.bam$", ".bw", file, ignore.case = TRUE),
  paired = FALSE,
  # mode=c("RNA-seq", "DNA-seq"), ## TODO: add strandness for RNA-seq
  method = c("deepTools", "Bioconductor"),
  cores = ezThreads()
) {
  method <- match.arg(method)

  if (method == "Bioconductor") {
    ## Better compatability; single-base resolutio;
    ## Slower; higher RAM consumption
    require(rtracklayer)
    require(GenomicAlignments)
    if (paired) {
      aligns <- readGAlignmentPairs(file)
    } else {
      aligns <- readGAlignments(file)
    }
    cov <- coverage(aligns)
    export.bw(cov, destination)
  } else if (method == "deepTools") {
    cmd <- paste(
      "bamCoverage",
      "--bam",
      file,
      "-o",
      destination,
      "--binSize 10 -of bigwig",
      "--normalizeUsing CPM",
      "-p",
      cores
    )
    ezSystem(cmd)
  }
  invisible(destination)
}

splitBamByRG <- function(inBam, mc.cores = ezThreads()) {
  # The following Rsamtools requires mapped Bam file.
  # require(Rsamtools)
  # header <- scanBamHeader(inBam)
  # readGroupNames <- header[[inBam]]$text[names(header[[inBam]]$text) == "@RG"]
  # readGroupNames <- sub("ID:", "", sapply(readGroupNames, "[", 1))
  cwd <- getwd()
  inBam <- normalizePath(inBam)

  # This is more general approach for uBam and mapped Bam.
  tempDIR <- file.path(cwd, paste("samtools-split", Sys.getpid(), sep = "-"))
  setwdNew(tempDIR)

  cmd <- paste("samtools split -@", mc.cores, "-f %!", inBam)
  ezSystem(cmd)

  setwd(cwd)

  readGroupNames <- list.files(path = tempDIR)
  bamFns <- paste0(readGroupNames, ".bam")
  file.rename(from = file.path(tempDIR, readGroupNames), to = bamFns)
  names(bamFns) <- readGroupNames

  # Clearning
  unlink(tempDIR, recursive = TRUE)

  return(bamFns)
}

mergeBamAlignments <- function(
  alignedBamFn,
  unmappedBamFn,
  outputBamFn,
  fastaFn,
  keepUnmapped = FALSE
) {
  ## Use . as tmp dir. Big bam generates big tmp files.
  cmd <- paste(
    prepareJavaTools("picard"),
    "MergeBamAlignment",
    paste0("ALIGNED=", alignedBamFn),
    paste0("UNMAPPED=", unmappedBamFn),
    paste0("O=", outputBamFn),
    paste0("R=", fastaFn)
  )
  ezSystem(cmd)
  if (!keepUnmapped) {
    tempBam <- tempfile(
      pattern = "noUnmapped",
      tmpdir = ".",
      fileext = ".bam"
    )
    cmd <- paste(
      "samtools view -F 4 -h -b",
      outputBamFn,
      ">",
      tempBam
    )
    ezSystem(cmd)
    file.remove(outputBamFn)
    file.rename(from = tempBam, to = outputBamFn)
  }
  invisible(outputBamFn)
}

posSpecErrorBam <- function(bamGA, genomeFn) {
  require(GenomicAlignments)
  require(Hmisc)
  require(BSgenome)
  require(Rsamtools)

  nMaxReads <- 100000

  what <- c("qname", "seq")
  stopifnot(all(what %in% colnames(mcols(bamGA))))

  hasGap <- grepl("I|D|N", cigar(bamGA))
  readLength <- qwidth(bamGA)
  genomeFnIndixe <- FaFile(genomeFn)
  genomeLengths <- seqlengths(genomeFnIndixe)

  isOutOfRange <- start(bamGA) + readLength - 1 >
    genomeLengths[as.character(seqnames(bamGA))] |
    start(bamGA) < readLength
  # this is very conservative; needed because there might be clipped bases in the beginning
  if (any(isOutOfRange)) {
    ezWrite("#reads out of range: ", sum(isOutOfRange))
    idx <- which(isOutOfRange)
    idx <- idx[1:min(10, length(idx))]
    badAlignments <- data.frame(
      pos = start(bamGA)[idx],
      cigar = cigar(bamGA)[idx],
      width = readLength[idx]
    )
    print(badAlignments)
  }

  indexKeep <- which(!hasGap & !isOutOfRange)
  # indexKeep <- which(!hasGap)
  if (length(indexKeep) > nMaxReads) {
    indexKeep <- sample(indexKeep, size = nMaxReads, replace = FALSE)
  }
  ezWrite(
    "#alignments: ",
    length(bamGA),
    " #valid alignments: ",
    sum(!hasGap & !isOutOfRange),
    " #used:",
    length(indexKeep)
  )
  if (length(indexKeep) == 0) {
    return(NULL)
  }
  bamGA <- bamGA[indexKeep]
  bamGR <- granges(bamGA) ## used to record the actual mapped coordinates of seq

  ## adjust the start POS according to H and/or S
  tempCigar <- str_extract(cigar(bamGA), "^(\\d+H)?(\\d+S)?\\d+M")
  ## get the number of H at the beginning
  noOfH <- explodeCigarOpLengths(tempCigar, ops = "H")
  stopifnot(!any(lengths(noOfH) > 1))
  noOfH[lengths(noOfH) == 0] <- 0
  noOfH <- as.integer(noOfH)

  ## get the number of S at the beginning
  noOfS <- explodeCigarOpLengths(tempCigar, ops = "S")
  stopifnot(!any(lengths(noOfS) > 1))
  noOfS[lengths(noOfS) == 0] <- 0
  noOfS <- as.integer(noOfS)
  nBeginClipped <- noOfH + noOfS
  start(bamGR) <- start(bamGR) - nBeginClipped

  ## add - to the begin and end of SEQ for revComplement
  Xbegin <- makeNstr("-", noOfH)
  tempCigar <- str_extract(cigar(bamGA), "\\d+[^HS](\\d+S)?(\\d+H)?$")
  noOfH <- explodeCigarOpLengths(tempCigar, ops = "H")
  stopifnot(!any(lengths(noOfH) > 1))
  noOfH[lengths(noOfH) == 0] <- 0
  noOfH <- as.integer(noOfH)

  Xend <- makeNstr("-", noOfH)
  noOfS <- explodeCigarOpLengths(tempCigar, ops = "S")
  stopifnot(!any(lengths(noOfS) > 1))
  noOfS[lengths(noOfS) == 0] <- 0
  noOfS <- as.integer(noOfS)
  nEndClipped <- noOfH + noOfS
  mcols(bamGR)$seq <- DNAStringSet(paste0(Xbegin, mcols(bamGA)$seq, Xend))
  end(bamGR) <- end(bamGR) + nEndClipped

  indexNeg <- as.logical(strand(bamGA) == "-")
  mcols(bamGR)$seq[indexNeg] <- reverseComplement(mcols(bamGR)$seq[indexNeg])

  seqChar <- strsplit(as.character(mcols(bamGR)$seq), "")
  readLength <- lengths(seqChar)

  maxLength <- quantile(readLength, 0.95)
  if (maxLength < max(readLength)) {
    readLength[readLength > maxLength] <- maxLength
    seqChar <- mapply(
      function(x, l) {
        x[1:l]
      },
      seqChar,
      readLength
    )
    bamGR <- resize(bamGR, width = readLength, fix = "start")
  }

  # referenceChar <- genome[bamGR]
  referenceChar <- getSeq(genomeFnIndixe, bamGR)
  stopifnot(all(width(referenceChar) == width(bamGR))) # check the position change is wright
  referenceChar <- strsplit(as.character(referenceChar), "")

  # assuming we have unique read length and set it to the maximal read length here.
  nEndTrimmed <- maxLength - readLength
  trimmedMatrix <- mapply(
    function(readLength, nEndTrimmed) {
      rep(c(FALSE, TRUE), c(readLength, nEndTrimmed))
    },
    readLength,
    nEndTrimmed,
    SIMPLIFY = FALSE
  )
  ## build a clippedMatrix to record the clipped character
  nNormal <- readLength - nBeginClipped - nEndClipped
  clippedMatrixNoTrim <- mapply(
    function(nBeginClipped, nNormal, nEndClipped) {
      rep(
        c(TRUE, FALSE, TRUE),
        c(nBeginClipped, nNormal, nEndClipped)
      )
    },
    nBeginClipped,
    nNormal,
    nEndClipped,
    SIMPLIFY = FALSE
  )

  matchMatrix <- mapply("==", referenceChar, seqChar, SIMPLIFY = FALSE)

  clippedMatrixNoTrim[indexNeg] <- lapply(clippedMatrixNoTrim[indexNeg], rev)
  clippedMatrix <- mapply(
    function(clippedMatrixNoTrim, nEndTrimmed) {
      c(clippedMatrixNoTrim, rep(FALSE, nEndTrimmed))
    },
    clippedMatrixNoTrim,
    nEndTrimmed,
    SIMPLIFY = FALSE
  )

  matchMatrix <- sapply(matchMatrix, function(x) {
    x[1:maxLength]
  })
  clippedMatrix <- matrix(unlist(clippedMatrix), ncol = length(clippedMatrix))
  clippedRate <- rowMeans(clippedMatrix, na.rm = TRUE)
  trimmedMatrix <- matrix(unlist(trimmedMatrix), ncol = length(trimmedMatrix))
  trimmedRate <- rowMeans(trimmedMatrix, na.rm = TRUE)
  matchMatrix[clippedMatrix] <- NA
  errorRate <- 1 - rowMeans(matchMatrix, na.rm = TRUE)
  names(errorRate) <- 1:nrow(matchMatrix)
  names(clippedRate) <- 1:nrow(matchMatrix)
  names(trimmedRate) <- 1:nrow(matchMatrix)
  return(list(
    trimmedRate = trimmedRate,
    clippedRate = clippedRate,
    errorRate = errorRate
  ))
}

### Create the MD tag for BAM file and replace matches with "=" in seq
calmdBam <- function(
  bamFns,
  mdBamFns = sub("\\.bam$", "_md.bam", bamFns),
  genomeFn,
  mc.cores = 4L
) {
  stopifnot(length(bamFns) == length(mdBamFns))

  cmdsCalmd <- paste("samtools calmd -b -e", bamFns, genomeFn, ">", mdBamFns)
  cmdsIndex <- paste("samtools index", mdBamFns)

  indicesExist <- file.exists(mdBamFns)
  if (any(indicesExist)) {
    warning(
      "Some md Bam files already exist!",
      paste(mdBamFns[indicesExist], collapse = " ")
    )
    cmdsCalmd <- cmdsCalmd[!indicesExist]
    cmdsIndex <- cmdsIndex[!indicesExist]
  }
  ezMclapply(cmdsCalmd, ezSystem, mc.preschedule = FALSE, mc.cores = mc.cores)
  ezMclapply(cmdsIndex, ezSystem, mc.preschedule = FALSE, mc.cores = mc.cores)
  return(mdBamFns)
}

### position specific error with md bam
posSpecErrorMDBam <- function(bamGA) {
  hasIndel <- grepl("I|D", cigar(bamGA))
  bamGA <- bamGA[!hasIndel]

  seqs <- mcols(bamGA)$seq
  negStrand <- strand(bamGA) == "-"
  seqs[negStrand] <- reverseComplement(seqs[negStrand])
  consensusCount <- consensusMatrix(seqs)
  errorRate <- 1 - consensusCount["-", ] / colSums(consensusCount)
  return(errorRate)
}
