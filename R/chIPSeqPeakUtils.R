###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## Utilities for the ChIPSeqPeakComparison app.
## Everything here takes plain arguments (no `param` object) so it can be
## unit-tested without a SUSHI job. See app-chIPSeqPeakComparison.R for the
## EzApp wiring and inst/templates/ChIPSeqPeakComparison.Rmd for the report.

## The canonical mcols schema every readPeakFile() result carries. Downstream
## code (harmonizePeaks, QC, consensus) relies on these columns existing, in
## this order, with NA where a format does not supply them.
.peakMcolSchema <- c("name", "score", "signalValue", "pValue", "qValue",
                     "summit", "annotation", "geneSymbol", "distanceToTSS")

## Case-insensitive "first candidate that matches a column name" lookup.
## Returns the matched column name (original casing) or NA_character_.
.matchCol <- function(cols, candidates) {
  lc <- tolower(cols)
  for (cand in tolower(candidates)) {
    hit <- which(lc == cand)
    if (length(hit) > 0L) return(cols[hit[1L]])
  }
  return(NA_character_)
}

## Coerce a vector to numeric, mapping MACS/ENCODE "unavailable" sentinels
## (-1) and non-numeric junk to NA.
.asPeakNumeric <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x[is.finite(x) & x == -1] <- NA_real_
  x
}

## Assemble a GRanges with the canonical schema from already-1-based coordinates.
.buildPeakGRanges <- function(seqnames, start1, end, strand = NULL,
                              name = NULL, score = NULL, signalValue = NULL,
                              pValue = NULL, qValue = NULL, summit = NULL,
                              annotation = NULL, geneSymbol = NULL,
                              distanceToTSS = NULL) {
  n <- length(seqnames)
  fillNum <- function(x) if (is.null(x) || length(x) == 0L) rep(NA_real_, n) else .asPeakNumeric(x)
  fillChr <- function(x) if (is.null(x) || length(x) == 0L) rep(NA_character_, n) else as.character(x)
  if (is.null(strand) || length(strand) == 0L) strand <- rep("*", n)
  strand <- as.character(strand)
  strand[is.na(strand) | strand == "." | !strand %in% c("+", "-", "*")] <- "*"
  gr <- GenomicRanges::GRanges(
    seqnames = as.character(seqnames),
    ranges = IRanges::IRanges(start = as.integer(start1), end = as.integer(end)),
    strand = strand)
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    name          = fillChr(name),
    score         = fillNum(score),
    signalValue   = fillNum(signalValue),
    pValue        = fillNum(pValue),
    qValue        = fillNum(qValue),
    summit        = fillNum(summit),
    annotation    = fillChr(annotation),
    geneSymbol    = fillChr(geneSymbol),
    distanceToTSS = fillNum(distanceToTSS))
  if (n > 0L && (is.null(name) || all(is.na(S4Vectors::mcols(gr)$name)))) {
    S4Vectors::mcols(gr)$name <- paste0("peak_", seq_len(n))
  }
  gr
}

## Empty GRanges carrying the schema (used for empty peak files).
.emptyPeakGRanges <- function() {
  .buildPeakGRanges(character(0), integer(0), integer(0))
}

## Sniff a peak-file format from extension plus a header peek.
## Returns one of: "narrowPeak", "broadPeak", "bed", "macsXls", "xlsx".
.sniffPeakFormat <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "xlsx") return("xlsx")
  if (ext == "narrowpeak") return("narrowPeak")
  if (ext == "broadpeak") return("broadPeak")
  if (ext == "xls") {
    ## MACS .xls is a TSV with '#' comment lines and a 'chr start end ...'
    ## header. A genuine binary .xls would fail this peek -> treat as xlsx.
    con <- file(file, "r"); on.exit(close(con))
    lines <- readLines(con, n = 50L, warn = FALSE)
    hdr <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]
    if (length(hdr) > 0L && grepl("abs_summit|fold_enrichment|start\\s+end", hdr[1L], ignore.case = TRUE)) {
      return("macsXls")
    }
    return("xlsx")
  }
  ## .bed (or unknown text): decide narrow/broad/plain by peeking column count
  ## of the first non-header data line.
  con <- file(file, "r"); on.exit(close(con))
  lines <- readLines(con, n = 200L, warn = FALSE)
  dat <- lines[!grepl("^(#|track|browser)", lines) & nzchar(trimws(lines))]
  if (length(dat) == 0L) return("bed")
  nc <- length(strsplit(dat[1L], "\t")[[1]])
  if (nc == 1L) nc <- length(strsplit(trimws(dat[1L]), "\\s+")[[1]])
  if (nc >= 10L) return("narrowPeak")
  if (nc == 9L)  return("broadPeak")
  return("bed")
}

## Read a BED / narrowPeak / broadPeak file (0-based half-open) into 1-based
## GRanges with the canonical schema.
.readBedLike <- function(file, kind, isZeroBased) {
  con <- file(file, "r")
  lines <- readLines(con, warn = FALSE); close(con)
  keep <- !grepl("^(#|track|browser)", lines) & nzchar(trimws(lines))
  lines <- lines[keep]
  if (length(lines) == 0L) return(.emptyPeakGRanges())
  tab <- data.table::fread(text = paste(lines, collapse = "\n"),
                           header = FALSE, sep = "\t", data.table = FALSE,
                           showProgress = FALSE)
  nc <- ncol(tab)
  zeroBased <- if (is.na(isZeroBased)) TRUE else isZeroBased  # BED-like defaults 0-based
  start1 <- tab[[2]] + (if (zeroBased) 1L else 0L)
  end    <- tab[[3]]
  name   <- if (nc >= 4L) tab[[4]] else NULL
  score  <- if (nc >= 5L) tab[[5]] else NULL
  strand <- if (nc >= 6L) tab[[6]] else NULL
  signalValue <- pValue <- qValue <- summit <- NULL
  if (kind %in% c("narrowPeak", "broadPeak") && nc >= 9L) {
    signalValue <- tab[[7]]
    pValue      <- tab[[8]]
    qValue      <- tab[[9]]
  }
  if (kind == "narrowPeak" && nc >= 10L) {
    offset <- .asPeakNumeric(tab[[10]])          # summit offset from (0-based) start
    summit <- start1 + offset                    # 1-based absolute summit
  }
  .buildPeakGRanges(tab[[1]], start1, end, strand = strand, name = name,
                    score = score, signalValue = signalValue,
                    pValue = pValue, qValue = qValue, summit = summit)
}

## Read a MACS .xls (TSV, 1-based, '#'-commented) into GRanges.
.readMacsXls <- function(file) {
  con <- file(file, "r"); lines <- readLines(con, warn = FALSE); close(con)
  keep <- !grepl("^#", lines) & nzchar(trimws(lines))
  lines <- lines[keep]
  if (length(lines) == 0L) return(.emptyPeakGRanges())
  tab <- data.table::fread(text = paste(lines, collapse = "\n"),
                           header = TRUE, sep = "\t", data.table = FALSE,
                           check.names = FALSE, showProgress = FALSE)
  cn <- colnames(tab)
  getc <- function(cands) { m <- .matchCol(cn, cands); if (is.na(m)) NULL else tab[[m]] }
  .buildPeakGRanges(
    seqnames = getc(c("chr", "chrom", "seqnames")),
    start1   = getc(c("start")),                       # MACS xls is 1-based
    end      = getc(c("end", "stop")),
    name     = getc(c("name", "peak", "peakid")),
    score    = getc(c("pileup")),
    signalValue = getc(c("fold_enrichment", "foldenrichment")),
    pValue   = getc(c("-log10(pvalue)", "-log10pvalue", "log10pvalue")),
    qValue   = getc(c("-log10(qvalue)", "-log10qvalue", "log10qvalue")),
    summit   = getc(c("abs_summit", "summit")))
}

## Read an .xlsx peak table (ChIPpeakAnno / HOMER annotatePeaks / MACS export /
## manual). Detects coordinate columns case-insensitively; 1-based unless told.
.readXlsxPeaks <- function(file, sheet, isZeroBased) {
  tab <- as.data.frame(readxl::read_excel(file, sheet = sheet), stringsAsFactors = FALSE)
  cn <- colnames(tab)
  getc <- function(cands) { m <- .matchCol(cn, cands); if (is.na(m)) NULL else tab[[m]] }
  chrCol   <- getc(c("chr", "chrom", "seqnames", "chromosome"))
  startCol <- getc(c("start", "chromstart", "peakstart"))
  endCol   <- getc(c("end", "stop", "chromend", "peakend"))
  ## HOMER annotatePeaks has PeakID in col1 and Chr/Start/End in 2-4; if the
  ## coordinate columns were not found by name, fall back to that layout.
  if (is.null(chrCol) && ncol(tab) >= 4L &&
      grepl("peak", tolower(cn[1L]))) {
    chrCol <- tab[[2]]; startCol <- tab[[3]]; endCol <- tab[[4]]
  }
  if (is.null(chrCol) || is.null(startCol) || is.null(endCol)) {
    stop("readPeakFile: could not detect chr/start/end columns in xlsx '",
         basename(file), "'. Columns seen: ", paste(cn, collapse = ", "))
  }
  zeroBased <- if (is.na(isZeroBased)) FALSE else isZeroBased  # xlsx defaults 1-based
  start1 <- as.numeric(startCol) + (if (zeroBased) 1L else 0L)
  .buildPeakGRanges(
    seqnames = chrCol, start1 = start1, end = endCol,
    strand = getc(c("strand")),
    name   = getc(c("name", "peakid", "peak", "id")),
    score  = getc(c("score", "pileup")),
    signalValue = getc(c("signalvalue", "fold_enrichment", "foldenrichment", "fold change")),
    pValue = getc(c("-log10(pvalue)", "pvalue", "p-value", "log10pvalue")),
    qValue = getc(c("-log10(qvalue)", "qvalue", "q-value", "fdr", "log10qvalue")),
    summit = getc(c("abs_summit", "summit")),
    annotation = getc(c("annotation", "feature")),
    geneSymbol = getc(c("gene name", "genesymbol", "symbol", "nearest gene", "gene")),
    distanceToTSS = getc(c("distance to tss", "distancetotss", "disttotss", "distance")))
}

#' @title Read a ChIP-seq peak file into a canonical GRanges
#' @description Normalises BED / narrowPeak / broadPeak / MACS `.xls` / `.xlsx`
#'   peak files into a \code{GRanges} carrying a stable \code{mcols} schema
#'   (\code{name, score, signalValue, pValue, qValue, summit, annotation,
#'   geneSymbol, distanceToTSS}), with \code{NA} where a format does not supply
#'   a field. Coordinates are returned 1-based: BED-like inputs are shifted
#'   from their 0-based half-open starts, MACS/xlsx inputs are assumed 1-based.
#'   \code{pValue}/\code{qValue} follow the MACS/ENCODE convention of being
#'   \code{-log10}-scaled.
#' @param file path to the peak file.
#' @param format one of \code{"auto"}, \code{"narrowPeak"}, \code{"broadPeak"},
#'   \code{"bed"}, \code{"macsXls"}, \code{"xlsx"}. \code{"auto"} sniffs by
#'   extension plus a header peek.
#' @param sheet sheet index or name for \code{.xlsx} (default first sheet).
#' @param isZeroBased override coordinate base: \code{NA} uses the per-format
#'   default (BED-like 0-based, MACS/xlsx 1-based); \code{TRUE}/\code{FALSE}
#'   forces it.
#' @return a \code{GRanges} with the canonical peak \code{mcols} schema.
#' @export
readPeakFile <- function(file, format = "auto", sheet = 1, isZeroBased = NA) {
  if (!file.exists(file)) stop("readPeakFile: file does not exist: ", file)
  if (identical(format, "auto")) format <- .sniffPeakFormat(file)
  gr <- switch(format,
    narrowPeak = .readBedLike(file, "narrowPeak", isZeroBased),
    broadPeak  = .readBedLike(file, "broadPeak", isZeroBased),
    bed        = .readBedLike(file, "bed", isZeroBased),
    macsXls    = .readMacsXls(file),
    xlsx       = .readXlsxPeaks(file, sheet, isZeroBased),
    stop("readPeakFile: unknown format '", format, "'"))
  ## de-duplicate names deterministically (some callers emit repeated names)
  nm <- S4Vectors::mcols(gr)$name
  if (anyDuplicated(nm)) S4Vectors::mcols(gr)$name <- make.unique(nm, sep = "_")
  gr
}

#' @title Harmonise peak sets against a reference seqinfo
#' @description Reconciles seqlevels against a reference \code{Seqinfo},
#'   attaches sequence lengths, and drops zero-width, out-of-bounds,
#'   off-reference and (optionally) blacklisted ranges -- recording every drop
#'   rather than discarding silently. Unlike \code{keepStandardChromosomes},
#'   this keeps every contig present in \code{seqinfo} (so bacterial plasmids
#'   and other non-"standard" sequences survive).
#' @param grl a named list or \code{GRangesList} of peak \code{GRanges} (as
#'   returned by \code{\link{readPeakFile}}).
#' @param seqinfo reference \code{Seqinfo} (e.g. built from the refBuild FASTA
#'   index) giving the valid seqlevels and their lengths.
#' @param blacklist optional \code{GRanges} of blacklisted regions; peaks
#'   overlapping it are dropped.
#' @return a list with \code{peaks} (a \code{GRangesList}) and \code{dropLog}
#'   (a \code{data.frame}: sample, reason, n).
#' @export
harmonizePeaks <- function(grl, seqinfo, blacklist = NULL) {
  if (methods::is(grl, "GRanges")) grl <- list(grl)
  nms <- names(grl)
  if (is.null(nms)) nms <- paste0("sample", seq_along(grl))
  refLevels <- GenomeInfoDb::seqnames(seqinfo)
  dropLog <- list()
  addDrop <- function(sample, reason, n) {
    if (n > 0L) dropLog[[length(dropLog) + 1L]] <<-
      data.frame(sample = sample, reason = reason, n = n, stringsAsFactors = FALSE)
  }
  out <- vector("list", length(grl))
  for (i in seq_along(grl)) {
    gr <- grl[[i]]; s <- nms[i]
    n0 <- length(gr)
    if (n0 == 0L) { out[[i]] <- gr; next }
    ## Try a seqlevelsStyle reconciliation, but never let it error the run
    ## (it is unreliable for e.g. bacterial RefSeq accessions).
    if (!all(as.character(GenomeInfoDb::seqnames(gr)) %in% refLevels)) {
      try({
        suppressWarnings(GenomeInfoDb::seqlevelsStyle(gr) <-
          GenomeInfoDb::seqlevelsStyle(seqinfo)[1])
      }, silent = TRUE)
    }
    ## restrict to reference contigs
    offRef <- !as.character(GenomeInfoDb::seqnames(gr)) %in% refLevels
    addDrop(s, "off-reference contig", sum(offRef))
    gr <- gr[!offRef]
    if (length(gr) == 0L) { out[[i]] <- gr; next }
    gr <- GenomeInfoDb::keepSeqlevels(gr, unique(as.character(GenomeInfoDb::seqnames(gr))),
                                      pruning.mode = "coarse")
    GenomeInfoDb::seqlevels(gr) <- refLevels
    GenomeInfoDb::seqinfo(gr) <- seqinfo
    ## zero-width
    zw <- BiocGenerics::width(gr) == 0L
    addDrop(s, "zero-width", sum(zw)); gr <- gr[!zw]
    ## out-of-bounds (start < 1 or end > seqlength)
    if (length(gr) > 0L) {
      slen <- GenomeInfoDb::seqlengths(seqinfo)[as.character(GenomeInfoDb::seqnames(gr))]
      oob <- BiocGenerics::start(gr) < 1L | (!is.na(slen) & BiocGenerics::end(gr) > slen)
      addDrop(s, "out-of-bounds", sum(oob)); gr <- gr[!oob]
    }
    ## blacklist
    if (!is.null(blacklist) && length(gr) > 0L) {
      bl <- IRanges::overlapsAny(gr, blacklist, ignore.strand = TRUE)
      addDrop(s, "blacklisted", sum(bl)); gr <- gr[!bl]
    }
    out[[i]] <- gr
  }
  names(out) <- nms
  dropDf <- if (length(dropLog)) do.call(rbind, dropLog) else
    data.frame(sample = character(0), reason = character(0), n = integer(0))
  list(peaks = GenomicRanges::GRangesList(out), dropLog = dropDf)
}

## ----------------------------------------------------------------------------
## QC layer. Every function tolerates a missing/NULL BAM by returning an
## NA-filled row of the correct shape (never an error), so the report can
## degrade section-by-section.
## ----------------------------------------------------------------------------

.usableFile <- function(f) !is.null(f) && !is.na(f) && nzchar(f) && file.exists(f)

#' @title Fraction of reads in peaks (FRiP)
#' @description Counts reads (or read pairs) overlapping \code{peaks} over the
#'   total, for one BAM.
#' @param bamFile path to an indexed BAM, or NULL.
#' @param peaks a \code{GRanges} of peak regions.
#' @param paired logical; if TRUE reads are counted as pairs.
#' @return a one-row \code{data.frame}: totalReads, readsInPeaks, frip.
#' @export
computeFrip <- function(bamFile, peaks, paired = FALSE) {
  na <- data.frame(totalReads = NA_real_, readsInPeaks = NA_real_, frip = NA_real_)
  if (!.usableFile(bamFile)) return(na)
  ## total mapped reads from the index (no read scan)
  idx <- Rsamtools::idxstatsBam(bamFile)
  total <- sum(idx$mapped)
  if (total == 0L) return(na)
  ## reads in peaks: count records overlapping the reduced (disjoint) peak set
  ## via the BAM index -- memory-safe for multi-GB BAMs.
  red <- GenomicRanges::reduce(peaks, ignore.strand = TRUE)
  common <- intersect(as.character(GenomeInfoDb::seqnames(red)),
                      idx$seqnames[idx$mapped > 0])
  red <- red[as.character(GenomeInfoDb::seqnames(red)) %in% common]
  if (length(red) == 0L) return(data.frame(totalReads = total, readsInPeaks = 0, frip = 0))
  cnt <- Rsamtools::countBam(bamFile,
    param = Rsamtools::ScanBamParam(which = red,
              flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)))
  inPk <- sum(cnt$records)
  data.frame(totalReads = total, readsInPeaks = inPk, frip = inPk / total)
}

## scan a (subsampled) set of BAM fields via a yieldSize-limited BamFile
.scanBamSubsample <- function(bamFile, what, nSub = 1e5) {
  bf <- Rsamtools::BamFile(bamFile, yieldSize = nSub)
  open(bf); on.exit(close(bf))
  Rsamtools::scanBam(bf, param = Rsamtools::ScanBamParam(what = what))[[1]]
}

#' @title Per-sample alignment statistics
#' @description Total/mapped reads (from the BAM index), plus MAPQ>=30 fraction,
#'   duplication rate and median fragment size estimated from a subsample.
#' @param bamFiles named character vector of BAM paths (entries may be NA).
#' @param cores number of cores for \code{ezMclapply}.
#' @param nSub subsample size for the flag/MAPQ/insert scan.
#' @return a \code{data.frame}, one row per input, with NA where a BAM is absent.
#' @export
alignmentStats <- function(bamFiles, cores = 1L, nSub = 1e5) {
  nms <- names(bamFiles); if (is.null(nms)) nms <- as.character(seq_along(bamFiles))
  one <- function(i) {
    f <- bamFiles[[i]]
    row <- data.frame(sample = nms[i], totalReads = NA_real_, mappedReads = NA_real_,
                      mappedFraction = NA_real_, mapq30Fraction = NA_real_,
                      duplicationRate = NA_real_, medianFragmentSize = NA_real_,
                      stringsAsFactors = FALSE)
    if (!.usableFile(f)) return(row)
    idx <- Rsamtools::idxstatsBam(f)
    row$mappedReads <- sum(idx$mapped)
    row$totalReads  <- sum(idx$mapped) + sum(idx$unmapped)
    row$mappedFraction <- if (row$totalReads > 0) row$mappedReads / row$totalReads else NA_real_
    sb <- .scanBamSubsample(f, c("mapq", "flag", "isize"), nSub)
    mapq <- sb$mapq; flag <- sb$flag
    if (length(flag)) {
      mapped <- bitwAnd(flag, 0x4L) == 0L
      row$mapq30Fraction  <- mean(mapq[mapped] >= 30, na.rm = TRUE)
      row$duplicationRate <- mean(bitwAnd(flag[mapped], 0x400L) != 0L)
      isize <- abs(sb$isize[mapped]); isize <- isize[isize > 0 & is.finite(isize)]
      row$medianFragmentSize <- if (length(isize)) stats::median(isize) else NA_real_
    }
    row
  }
  do.call(rbind, ezMclapplyOrLapply(seq_along(bamFiles), one, cores))
}

#' @title Library complexity (NRF, PBC1, PBC2)
#' @description ENCODE library-complexity metrics from (subsampled) uniquely
#'   mapped read 5' positions.
#' @param bamFile path to an indexed BAM, or NULL.
#' @param nSub subsample size.
#' @return a one-row \code{data.frame}: NRF, PBC1, PBC2.
#' @export
libraryComplexity <- function(bamFile, chunkSize = 5e6) {
  na <- data.frame(NRF = NA_real_, PBC1 = NA_real_, PBC2 = NA_real_)
  if (!.usableFile(bamFile)) return(na)
  ## Stream the WHOLE BAM in chunks, accumulating counts per unique 5' position
  ## (rname, pos, strand). Subsampling the first-N reads of a coordinate-sorted
  ## BAM biases these metrics badly (reads cluster on one contig), so we scan
  ## all reads; memory is bounded by the number of distinct positions.
  bf <- Rsamtools::BamFile(bamFile, yieldSize = chunkSize)
  open(bf); on.exit(close(bf))
  agg <- NULL
  repeat {
    sb <- Rsamtools::scanBam(bf, param = Rsamtools::ScanBamParam(
      what = c("rname", "pos", "strand", "flag")))[[1]]
    if (length(sb$pos) == 0L) break
    mapped <- bitwAnd(sb$flag, 0x4L) == 0L & !is.na(sb$pos)
    if (!any(mapped)) next
    dt <- data.table::data.table(r = as.integer(sb$rname[mapped]),
                                 p = sb$pos[mapped],
                                 s = as.integer(sb$strand[mapped]))
    cur <- dt[, .N, by = c("r", "p", "s")]
    agg <- if (is.null(agg)) cur else
      rbind(agg, cur)[, list(N = sum(N)), by = c("r", "p", "s")]
  }
  if (is.null(agg) || nrow(agg) == 0L) return(na)
  total <- sum(agg$N); distinct <- nrow(agg)
  M1 <- sum(agg$N == 1L); M2 <- sum(agg$N == 2L)
  data.frame(NRF = distinct / total,
             PBC1 = if (distinct > 0) M1 / distinct else NA_real_,
             PBC2 = if (M2 > 0) M1 / M2 else NA_real_)
}

#' @title Per-sample peak statistics
#' @param grl a named list / \code{GRangesList} of peak \code{GRanges}.
#' @param genomeSize total genome size in bp (for percent-genome).
#' @return a \code{data.frame}, one row per sample.
#' @export
peakStats <- function(grl, genomeSize = NA_real_) {
  nms <- names(grl); if (is.null(nms)) nms <- as.character(seq_along(grl))
  do.call(rbind, lapply(seq_along(grl), function(i) {
    w <- BiocGenerics::width(grl[[i]]); bp <- sum(as.numeric(w))
    data.frame(sample = nms[i], nPeaks = length(w),
               totalBp = bp,
               pctGenome = if (is.finite(genomeSize) && genomeSize > 0) 100 * bp / genomeSize else NA_real_,
               medianWidth = if (length(w)) stats::median(w) else NA_real_,
               meanWidth = if (length(w)) mean(w) else NA_real_,
               q25Width = if (length(w)) stats::quantile(w, 0.25) else NA_real_,
               q75Width = if (length(w)) stats::quantile(w, 0.75) else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

#' @title Genome-wide BigWig correlation (pure R)
#' @description Bins the genome and computes the mean BigWig score per bin for
#'   each track, then the between-sample correlation matrices. Replaces a
#'   deepTools multiBigwigSummary with a no-temp-file R implementation.
#' @param bwFiles named character vector of BigWig paths (entries may be NA).
#' @param binSize bin width in bp.
#' @param cores cores for \code{ezMclapply}.
#' @return a list: \code{binMatrix} (bins x samples), \code{spearman},
#'   \code{pearson}; or NULL if fewer than two usable tracks.
#' @export
bigwigCorrelation <- function(bwFiles, binSize = 1000L, cores = 1L) {
  ok <- vapply(bwFiles, .usableFile, logical(1))
  bwFiles <- bwFiles[ok]
  if (length(bwFiles) < 2L) return(NULL)
  si <- GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(bwFiles[[1]]))
  tiles <- GenomicRanges::tileGenome(
    GenomeInfoDb::seqlengths(si), tilewidth = binSize,
    cut.last.tile.in.chrom = TRUE)
  scoreOne <- function(f) {
    cov <- rtracklayer::import.bw(f, as = "RleList")
    cov <- cov[GenomeInfoDb::seqlevels(tiles)]
    GenomicRanges::binnedAverage(tiles, cov, "s")$s
  }
  mat <- do.call(cbind, ezMclapplyOrLapply(seq_along(bwFiles),
                                           function(i) scoreOne(bwFiles[[i]]), cores))
  colnames(mat) <- names(bwFiles)
  list(binMatrix = mat,
       spearman = stats::cor(mat, method = "spearman", use = "pairwise.complete.obs"),
       pearson  = stats::cor(mat, method = "pearson",  use = "pairwise.complete.obs"))
}

#' @title Strand cross-correlation QC (native, infra-independent)
#' @description Estimates fragment length and the ENCODE NSC/RSC metrics from
#'   the strand cross-correlation of (subsampled) 5' read positions, binned for
#'   tractability. This is the fallback path when phantompeakqualtools'
#'   \code{run_spp.R} is unavailable; the \code{method} column records which
#'   estimator produced the numbers so they are never silently interchanged.
#' @param bamFile path to an indexed BAM, or NULL.
#' @param binSize resolution of the binned strand coverage.
#' @param maxShift maximum shift (bp) to scan.
#' @param nSub subsample of reads.
#' @return a one-row \code{data.frame}: estFragLen, NSC, RSC, method.
#' @export
crossCorrelationMetrics <- function(bamFile, binSize = 10L, maxShift = 1000L,
                                     nSub = 3e5) {
  na <- data.frame(estFragLen = NA_real_, NSC = NA_real_, RSC = NA_real_,
                   method = NA_character_, stringsAsFactors = FALSE)
  if (!.usableFile(bamFile)) return(na)
  sb <- .scanBamSubsample(bamFile, c("rname", "pos", "strand", "qwidth", "flag"), nSub)
  mapped <- bitwAnd(sb$flag, 0x4L) == 0L & !is.na(sb$pos)
  if (sum(mapped) < 100L) return(na)
  rname <- as.character(sb$rname[mapped]); pos <- sb$pos[mapped]
  strand <- as.character(sb$strand[mapped]); qw <- sb$qwidth[mapped]
  readLen <- stats::median(qw, na.rm = TRUE); if (!is.finite(readLen)) readLen <- 50
  maxBin <- ceiling(maxShift / binSize)
  ## build per-contig binned plus/minus signals, padded by maxBin zeros so the
  ## shifted cross-correlation does not bleed across contigs.
  plus <- list(); minus <- list()
  for (chr in unique(rname)) {
    sel <- rname == chr
    b <- (pos[sel] %/% binSize) + 1L
    nb <- max(b) + maxBin
    p <- tabulate(b[strand[sel] == "+"], nb)
    m <- tabulate(b[strand[sel] == "-"], nb)
    plus[[chr]] <- c(p, integer(maxBin)); minus[[chr]] <- c(m, integer(maxBin))
  }
  P <- unlist(plus, use.names = FALSE); M <- unlist(minus, use.names = FALSE)
  L <- length(P)
  cc <- vapply(0:maxBin, function(d) {
    suppressWarnings(stats::cor(P[1:(L - d)], M[(1 + d):L]))
  }, numeric(1))
  cc[!is.finite(cc)] <- NA_real_
  shifts <- (0:maxBin) * binSize
  ccMin <- min(cc, na.rm = TRUE)
  readBin <- which.min(abs(shifts - readLen))          # phantom peak at read length
  ccRead <- cc[readBin]
  fragRegion <- which(shifts > 2 * readLen)             # look past the phantom peak
  if (!length(fragRegion)) return(na)
  fragBin <- fragRegion[which.max(cc[fragRegion])]
  ccFrag <- cc[fragBin]
  data.frame(estFragLen = shifts[fragBin],
             NSC = ccFrag / ccMin,
             RSC = (ccFrag - ccMin) / (ccRead - ccMin),
             method = "native-strandCC", stringsAsFactors = FALSE)
}

#' @title Load the mark-type-specific QC threshold table
#' @param markType \code{"narrow"} or \code{"broad"}.
#' @return a tidy \code{data.frame}: metric, direction, min, good, note.
#' @export
qcThresholdTable <- function(markType = c("narrow", "broad")) {
  markType <- match.arg(markType)
  f <- getOption("ezRun_chIPThresholdsFile", "")
  if (!nzchar(f) || !file.exists(f))
    f <- system.file("extdata", "chipSeqPeakComparison", "qc_thresholds.tsv",
                     package = "ezRun")
  if (!nzchar(f) || !file.exists(f)) {
    ## not installed (e.g. running under load_all / a bare test session): search
    ## the source tree relative to plausible working directories.
    rel <- file.path("inst", "extdata", "chIPSeqPeakComparison", "qc_thresholds.tsv")
    cand <- c(rel, file.path("..", "..", rel), file.path("..", "..", "..", rel))
    hit <- cand[file.exists(cand)]
    if (!length(hit)) stop("qcThresholdTable: cannot locate qc_thresholds.tsv")
    f <- hit[1L]
  }
  tab <- utils::read.delim(f, stringsAsFactors = FALSE, na.strings = "NA")
  data.frame(metric = tab$metric, direction = tab$direction,
             min = tab[[paste0(markType, "_min")]],
             good = tab[[paste0(markType, "_good")]],
             note = tab$note, stringsAsFactors = FALSE)
}

#' @title Flag QC metrics against mark-type thresholds
#' @description Produces a long \code{data.frame} of (sample, metric, value,
#'   min, good, flag) for colour-coding in the report. \code{flag} is one of
#'   \code{"good"}, \code{"acceptable"}, \code{"low"}, or \code{NA} (no value or
#'   no threshold for this mark type).
#' @param qcWide a wide per-sample \code{data.frame} with a \code{sample}
#'   column and metric columns matching the threshold table.
#' @param markType \code{"narrow"} or \code{"broad"}.
#' @return a long, flagged \code{data.frame}.
#' @export
flagQcMetrics <- function(qcWide, markType = c("narrow", "broad")) {
  thr <- qcThresholdTable(markType)
  rows <- list()
  for (m in thr$metric) {
    if (!m %in% colnames(qcWide)) next
    tmin <- thr$min[thr$metric == m]; tgood <- thr$good[thr$metric == m]
    for (i in seq_len(nrow(qcWide))) {
      v <- suppressWarnings(as.numeric(qcWide[[m]][i]))
      flag <- if (is.na(v) || is.na(tmin)) NA_character_
              else if (!is.na(tgood) && v >= tgood) "good"
              else if (v >= tmin) "acceptable"
              else "low"
      rows[[length(rows) + 1L]] <- data.frame(
        sample = qcWide$sample[i], metric = m, value = v,
        min = tmin, good = tgood, flag = flag, stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(data.frame(sample = character(0), metric = character(0),
                                       value = numeric(0), min = numeric(0),
                                       good = numeric(0), flag = character(0)))
  do.call(rbind, rows)
}

## Use ezMclapply when the ezRun runtime is loaded, plain lapply otherwise
## (keeps the leaf functions runnable in a bare test session).
ezMclapplyOrLapply <- function(x, FUN, cores = 1L) {
  if (cores > 1L && exists("ezMclapply", mode = "function")) {
    ezMclapply(x, FUN, mc.cores = cores)
  } else {
    lapply(x, FUN)
  }
}

## ----------------------------------------------------------------------------
## Overlap / consensus / quantification layer.
## ----------------------------------------------------------------------------

#' @title Pairwise base-pair Jaccard index between peak sets
#' @description \code{J(a,b) = bp(intersect) / bp(union)} on reduced,
#'   strand-agnostic ranges.
#' @param grl named list / \code{GRangesList} of peak \code{GRanges}.
#' @param minOverlapBp reserved for count-based variants; the bp-Jaccard uses
#'   exact overlap widths and ignores it.
#' @return a symmetric numeric matrix with unit diagonal.
#' @export
pairwiseJaccard <- function(grl, minOverlapBp = 1L) {
  nms <- names(grl); if (is.null(nms)) nms <- as.character(seq_along(grl))
  red <- lapply(grl, function(g) GenomicRanges::reduce(g, ignore.strand = TRUE))
  n <- length(red); m <- matrix(1, n, n, dimnames = list(nms, nms))
  bp <- function(g) sum(as.numeric(BiocGenerics::width(g)))
  for (i in seq_len(n)) for (j in seq_len(n)) if (j > i) {
    inter <- bp(GenomicRanges::intersect(red[[i]], red[[j]], ignore.strand = TRUE))
    uni   <- bp(GenomicRanges::union(red[[i]], red[[j]], ignore.strand = TRUE))
    m[i, j] <- m[j, i] <- if (uni > 0) inter / uni else 0
  }
  m
}

#' @title Build a consensus peak set with an occupancy matrix
#' @description Reduces all peaks to candidate consensus regions, records which
#'   samples occupy each region (overlap >= \code{minOverlapBp}), and keeps
#'   regions occupied by at least \code{minSamples} samples.
#' @param grl named list / \code{GRangesList} of peak \code{GRanges}.
#' @param minSamples minimum number of samples occupying a region to keep it.
#' @param minOverlapBp minimum overlap (bp) for a sample to count as occupying.
#' @return a list: \code{consensus} (\code{GRanges} with an \code{occupancy}
#'   mcol) and \code{occupancyMatrix} (regions x samples, 0/1).
#' @export
buildConsensus <- function(grl, minSamples = 1L, minOverlapBp = 1L) {
  nms <- names(grl); if (is.null(nms)) nms <- as.character(seq_along(grl))
  pooled <- suppressWarnings(do.call(c, unname(lapply(grl, function(g) {
    g2 <- g; S4Vectors::mcols(g2) <- NULL; g2 }))))
  cand <- GenomicRanges::reduce(pooled, ignore.strand = TRUE)
  occ <- matrix(0L, length(cand), length(grl), dimnames = list(NULL, nms))
  for (j in seq_along(grl)) {
    hits <- GenomicRanges::findOverlaps(cand, grl[[j]], ignore.strand = TRUE,
                                        minoverlap = minOverlapBp)
    occ[unique(S4Vectors::queryHits(hits)), j] <- 1L
  }
  occCount <- rowSums(occ)
  keep <- occCount >= minSamples
  cons <- cand[keep]
  S4Vectors::mcols(cons)$occupancy <- occCount[keep]
  S4Vectors::mcols(cons)$name <- paste0("consensus_", seq_along(cons))
  list(consensus = cons, occupancyMatrix = occ[keep, , drop = FALSE])
}

## mean BigWig score per region (pure R, via Rle Views)
.scoreRegionsBw <- function(f, regions) {
  cov <- rtracklayer::import.bw(f, as = "RleList")
  res <- rep(NA_real_, length(regions))
  spl <- split(seq_along(regions), as.character(GenomeInfoDb::seqnames(regions)))
  for (chr in names(spl)) {
    if (!chr %in% names(cov)) next
    idx <- spl[[chr]]
    v <- IRanges::Views(cov[[chr]],
                        start = BiocGenerics::start(regions)[idx],
                        end   = pmin(BiocGenerics::end(regions)[idx], length(cov[[chr]])))
    res[idx] <- IRanges::viewMeans(v)
  }
  res
}

#' @title Quantify signal over regions from BigWig or BAM
#' @description Builds a region x sample matrix of signal. BigWig source (the
#'   default) uses mean track score per region; BAM source uses read counts
#'   (\code{summarizeOverlaps}, Union mode). Exactly one source is used: pass
#'   \code{bwFiles} for the BigWig path, or \code{bamFiles} for the BAM path.
#' @param regions a \code{GRanges} (e.g. the consensus set).
#' @param bamFiles named character vector of BAMs, or NULL.
#' @param bwFiles named character vector of BigWigs, or NULL.
#' @param cores cores for \code{ezMclapply}.
#' @param paired whether BAM reads are paired (BAM path only).
#' @return a numeric matrix, regions x samples.
#' @export
quantifyRegions <- function(regions, bamFiles = NULL, bwFiles = NULL,
                            cores = 1L, paired = FALSE) {
  if (!is.null(bwFiles)) {
    ok <- vapply(bwFiles, .usableFile, logical(1)); bwFiles <- bwFiles[ok]
    mat <- do.call(cbind, ezMclapplyOrLapply(seq_along(bwFiles),
      function(i) .scoreRegionsBw(bwFiles[[i]], regions), cores))
    colnames(mat) <- names(bwFiles)
    rownames(mat) <- if (!is.null(S4Vectors::mcols(regions)$name))
      S4Vectors::mcols(regions)$name else as.character(regions)
    return(mat)
  }
  if (!is.null(bamFiles)) {
    ok <- vapply(bamFiles, .usableFile, logical(1)); bamFiles <- bamFiles[ok]
    bfl <- Rsamtools::BamFileList(bamFiles)
    se <- GenomicAlignments::summarizeOverlaps(
      regions, bfl, mode = "Union", ignore.strand = TRUE, singleEnd = !paired)
    mat <- SummarizedExperiment::assay(se)
    colnames(mat) <- names(bamFiles)
    rownames(mat) <- if (!is.null(S4Vectors::mcols(regions)$name))
      S4Vectors::mcols(regions)$name else as.character(regions)
    return(mat)
  }
  stop("quantifyRegions: provide either bwFiles or bamFiles")
}

#' @title Normalise a signal / count matrix
#' @param mat region x sample numeric matrix.
#' @param method \code{"CPM"} (counts per million; needs \code{libSizes}),
#'   \code{"none"}, or \code{"quantile"} (\code{limma::normalizeQuantiles}).
#' @param libSizes per-sample library sizes for CPM; defaults to column sums.
#' @return a normalised matrix of the same shape.
#' @export
normalizeSignal <- function(mat, method = c("CPM", "none", "quantile"),
                            libSizes = NULL) {
  method <- match.arg(method)
  if (method == "none") return(mat)
  if (method == "CPM") {
    if (is.null(libSizes)) libSizes <- colSums(mat, na.rm = TRUE)
    sweep(mat, 2, pmax(libSizes, 1), "/") * 1e6
  } else {
    if (!requireNamespace("limma", quietly = TRUE))
      stop("normalizeSignal: quantile normalisation needs the 'limma' package")
    out <- limma::normalizeQuantiles(mat)
    dimnames(out) <- dimnames(mat); out
  }
}

## ----------------------------------------------------------------------------
## topN coverage layer (plan section 7).
## ----------------------------------------------------------------------------

#' @title Select the top-N peaks of a sample
#' @description Ranks by \code{rankBy} (descending) with ties broken by width,
#'   and returns the top \code{min(n, length(gr))} peaks. Warns when the request
#'   is truncated. \code{qValue}/\code{pValue} are -log10-scaled, so larger is
#'   more significant and the same descending rule applies.
#' @param gr a peak \code{GRanges} (canonical schema).
#' @param n number of peaks to keep.
#' @param rankBy one of \code{"signalValue"}, \code{"qValue"}, \code{"pValue"},
#'   \code{"score"}, \code{"width"}.
#' @return a \code{GRanges} of the selected peaks, ordered best-first.
#' @export
topNPeaks <- function(gr, n = 1000L, rankBy = c("signalValue", "qValue", "pValue",
                                                "score", "width")) {
  rankBy <- match.arg(rankBy)
  if (length(gr) == 0L) return(gr)
  metric <- if (rankBy == "width") BiocGenerics::width(gr) else S4Vectors::mcols(gr)[[rankBy]]
  if (is.null(metric) || all(is.na(metric))) {
    warning("topNPeaks: rankBy '", rankBy, "' is all-NA; falling back to width")
    metric <- BiocGenerics::width(gr)
  }
  ord <- order(metric, BiocGenerics::width(gr), decreasing = TRUE, na.last = TRUE)
  k <- min(n, length(gr))
  if (n > length(gr))
    warning("topNPeaks: requested ", n, " but only ", length(gr), " peaks available")
  gr[ord[seq_len(k)]]
}

## centre of each peak: summit where available, else midpoint (as width-1 GRanges)
.peakCentres <- function(gr) {
  s <- S4Vectors::mcols(gr)$summit
  mid <- as.integer((BiocGenerics::start(gr) + BiocGenerics::end(gr)) / 2)
  pos <- ifelse(is.na(s), mid, as.integer(s))
  GenomicRanges::GRanges(GenomeInfoDb::seqnames(gr),
                         IRanges::IRanges(pos, width = 1L),
                         strand = BiocGenerics::strand(gr),
                         seqinfo = GenomeInfoDb::seqinfo(gr))
}

#' @title Signal profile matrix over peak centres
#' @description Extracts a binned signal matrix (rows = regions, columns =
#'   positions around the centre) from a BigWig using
#'   \code{EnrichedHeatmap::normalizeToMatrix}. Regions are centred on the summit
#'   (or midpoint) and extended by \code{extend} bp each side. Pure R, no temp
#'   files.
#' @param bwFile path to a BigWig, or NULL (returns NULL).
#' @param regions a peak \code{GRanges}.
#' @param extend bp to extend each side of the centre.
#' @param binSize bin width in bp.
#' @return a \code{normalizeToMatrix} object, or NULL if the BigWig is absent.
#' @export
profileMatrix <- function(bwFile, regions, extend = 2000L, binSize = 50L) {
  if (!.usableFile(bwFile) || length(regions) == 0L) return(NULL)
  if (!requireNamespace("EnrichedHeatmap", quietly = TRUE))
    stop("profileMatrix needs the 'EnrichedHeatmap' package")
  centres <- .peakCentres(regions)
  win <- suppressWarnings(GenomicRanges::trim(
    GenomicRanges::resize(centres, width = 2L * extend + 1L, fix = "center")))
  signal <- rtracklayer::import.bw(bwFile, which = win)
  EnrichedHeatmap::normalizeToMatrix(
    signal, centres, value_column = "score", extend = extend, w = binSize,
    mean_mode = "w0", background = 0)
}

#' @title Genomic footprint of a set of peaks
#' @description Per-sample footprint: total bp, percent genome, width summary,
#'   and per-chromosome peak count / bp / bp-per-Mb enrichment. Also returns the
#'   coverage \code{RleList} for downstream ideogram / density plots.
#' @param grl named list / \code{GRangesList} of peak \code{GRanges}.
#' @param seqinfo reference \code{Seqinfo} (for chromosome lengths).
#' @return a list: \code{summary} (per-sample df), \code{perChrom} (long df),
#'   \code{coverage} (named list of \code{RleList}).
#' @export
genomeFootprint <- function(grl, seqinfo) {
  nms <- names(grl); if (is.null(nms)) nms <- as.character(seq_along(grl))
  genomeSize <- sum(as.numeric(GenomeInfoDb::seqlengths(seqinfo)))
  chromLen <- GenomeInfoDb::seqlengths(seqinfo)
  summ <- list(); perChrom <- list(); covs <- list()
  for (i in seq_along(grl)) {
    gr <- grl[[i]]; s <- nms[i]
    GenomeInfoDb::seqlevels(gr, pruning.mode = "coarse") <- GenomeInfoDb::seqnames(seqinfo)
    GenomeInfoDb::seqinfo(gr) <- seqinfo
    w <- BiocGenerics::width(gr); bp <- sum(as.numeric(w))
    summ[[i]] <- data.frame(sample = s, nPeaks = length(gr), totalBp = bp,
                            pctGenome = 100 * bp / genomeSize,
                            medianWidth = if (length(w)) stats::median(w) else NA_real_,
                            stringsAsFactors = FALSE)
    cnt <- table(factor(as.character(GenomeInfoDb::seqnames(gr)),
                        levels = names(chromLen)))
    bpByChr <- tapply(as.numeric(w), factor(as.character(GenomeInfoDb::seqnames(gr)),
                                            levels = names(chromLen)), sum)
    bpByChr[is.na(bpByChr)] <- 0
    perChrom[[i]] <- data.frame(
      sample = s, chrom = names(chromLen), nPeaks = as.integer(cnt),
      bp = as.numeric(bpByChr), chromLen = as.numeric(chromLen),
      bpPerMb = 1e6 * as.numeric(bpByChr) / as.numeric(chromLen),
      stringsAsFactors = FALSE)
    covs[[s]] <- GenomicRanges::coverage(gr)
  }
  list(summary = do.call(rbind, summ),
       perChrom = do.call(rbind, perChrom),
       coverage = covs)
}

#' @title Cumulative-signal concentration curve
#' @description Ranks peaks best-first by \code{rankBy} and returns the
#'   cumulative fraction of total signal vs the fraction of peaks. A flat curve
#'   (signal spread evenly) is itself a QC signal.
#' @param gr a peak \code{GRanges}.
#' @param rankBy ranking metric (see \code{topNPeaks}).
#' @return a \code{data.frame}: rank, fracPeaks, cumFracSignal.
#' @export
cumulativeSignalCurve <- function(gr, rankBy = "signalValue") {
  if (length(gr) == 0L) return(data.frame(rank = integer(0), fracPeaks = numeric(0),
                                          cumFracSignal = numeric(0)))
  v <- if (rankBy == "width") BiocGenerics::width(gr) else S4Vectors::mcols(gr)[[rankBy]]
  if (is.null(v) || all(is.na(v))) v <- BiocGenerics::width(gr)
  v <- sort(v[!is.na(v)], decreasing = TRUE)
  tot <- sum(v)
  data.frame(rank = seq_along(v), fracPeaks = seq_along(v) / length(v),
             cumFracSignal = if (tot > 0) cumsum(v) / tot else NA_real_)
}

#' @title Select the top-n most significant peaks across all samples
#' @description Pools every sample's peaks, ranks by \code{rankBy} (descending;
#'   \code{signalValue} = fold enrichment over background), and greedily keeps
#'   the top \code{n} that do not overlap each other -- i.e. the \code{n}
#'   strongest distinct binding loci in the whole dataset. These are the loci for
#'   which the report draws per-sample coverage tracks.
#' @param grl named list / \code{GRangesList} of peak \code{GRanges}.
#' @param n number of distinct loci to return (small, e.g. 3-10).
#' @param rankBy ranking metric (\code{signalValue}, \code{qValue},
#'   \code{pValue}, \code{score}, \code{width}).
#' @return a \code{GRanges} of \code{<= n} peaks, best-first, with added mcols
#'   \code{sourceSample} and \code{log2Enrichment}.
#' @export
topSignificantPeaks <- function(grl, n = 5L, rankBy = "signalValue") {
  nms <- names(grl); if (is.null(nms)) nms <- as.character(seq_along(grl))
  pooled <- suppressWarnings(do.call(c, unname(lapply(seq_along(grl), function(i) {
    g <- grl[[i]]; if (length(g) == 0L) return(g)
    S4Vectors::mcols(g)$sourceSample <- nms[i]; g
  }))))
  if (length(pooled) == 0L) return(pooled)
  metric <- if (rankBy == "width") BiocGenerics::width(pooled) else S4Vectors::mcols(pooled)[[rankBy]]
  if (is.null(metric) || all(is.na(metric))) metric <- BiocGenerics::width(pooled)
  ord <- order(metric, decreasing = TRUE, na.last = TRUE)
  chosen <- pooled[0]; keep <- integer(0)
  for (idx in ord) {
    if (length(keep) >= n) break
    if (length(chosen) == 0L ||
        !IRanges::overlapsAny(pooled[idx], chosen, ignore.strand = TRUE)) {
      keep <- c(keep, idx); chosen <- c(chosen, GenomicRanges::granges(pooled[idx]))
    }
  }
  res <- pooled[keep]
  m2 <- if (rankBy == "width") BiocGenerics::width(res) else S4Vectors::mcols(res)[[rankBy]]
  res <- res[order(m2, decreasing = TRUE)]
  sv <- S4Vectors::mcols(res)$signalValue
  S4Vectors::mcols(res)$log2Enrichment <- ifelse(is.na(sv) | sv <= 0, NA_real_, log2(sv))
  S4Vectors::mcols(res)$name <- paste0("topPeak_", seq_along(res))
  res
}

#' @title Pre-extract per-sample coverage around loci (for genome-track plots)
#' @description For each peak (extended by \code{flank} bp each side), extracts
#'   the binned mean BigWig coverage of every sample. Returns self-contained data
#'   so the report can draw Gviz coverage tracks without touching the BigWig
#'   files at render time.
#' @param peaks a \code{GRanges} of loci (e.g. \code{topSignificantPeaks} output).
#' @param bwFiles named character vector of BigWig paths (entries may be NA).
#' @param flank bp to extend each side of the peak.
#' @param binSize bin width in bp for the coverage.
#' @return a named list, one entry per peak: \code{peak}, \code{region},
#'   \code{binGR} (bin \code{GRanges}), \code{cov} (bins x samples matrix),
#'   \code{label}.
#' @export
extractLocusCoverage <- function(peaks, bwFiles, flank = 2000L, binSize = 50L) {
  ok <- vapply(bwFiles, .usableFile, logical(1)); bwFiles <- bwFiles[ok]
  if (length(peaks) == 0L || length(bwFiles) == 0L) return(list())
  covRle <- lapply(bwFiles, function(f) rtracklayer::import.bw(f, as = "RleList"))
  names(covRle) <- names(bwFiles)
  out <- lapply(seq_along(peaks), function(i) {
    pk <- peaks[i]
    region <- suppressWarnings(GenomicRanges::trim(GenomicRanges::resize(
      pk, width = BiocGenerics::width(pk) + 2L * flank, fix = "center")))
    chr <- as.character(GenomeInfoDb::seqnames(region))
    binGR <- unlist(GenomicRanges::tile(region, width = binSize))
    ## restrict every returned object to just this contig, so genome-track
    ## renderers (Gviz) do not choke on the other (plasmid/scaffold) seqlevels.
    pk     <- GenomeInfoDb::keepSeqlevels(pk, chr, pruning.mode = "coarse")
    region <- GenomeInfoDb::keepSeqlevels(region, chr, pruning.mode = "coarse")
    binGR  <- GenomeInfoDb::keepSeqlevels(binGR, chr, pruning.mode = "coarse")
    mat <- vapply(names(bwFiles), function(s) {
      rl <- covRle[[s]]
      if (!chr %in% names(rl)) return(rep(0, length(binGR)))
      v <- IRanges::Views(rl[[chr]], start = BiocGenerics::start(binGR),
                          end = pmin(BiocGenerics::end(binGR), length(rl[[chr]])))
      IRanges::viewMeans(v)
    }, numeric(length(binGR)))
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = length(bwFiles),
                                         dimnames = list(NULL, names(bwFiles)))
    sv <- S4Vectors::mcols(pk)$signalValue; l2 <- S4Vectors::mcols(pk)$log2Enrichment
    list(peak = pk, region = region, binGR = binGR, cov = mat,
         label = sprintf("%s  %s:%s-%s  (fold=%.1f, log2=%.1f%s)",
                         S4Vectors::mcols(pk)$name, chr,
                         format(BiocGenerics::start(pk), big.mark = ","),
                         format(BiocGenerics::end(pk), big.mark = ","),
                         ifelse(is.na(sv), NA, sv), ifelse(is.na(l2), NA, l2),
                         if (!is.null(S4Vectors::mcols(pk)$sourceSample))
                           paste0(", from ", S4Vectors::mcols(pk)$sourceSample) else ""))
  })
  names(out) <- S4Vectors::mcols(peaks)$name
  out
}

#' @title Annotate a consensus peak set with ChIPseeker
#' @param gr a \code{GRanges} to annotate (e.g. the consensus set).
#' @param txdb a \code{TxDb} object.
#' @param tssRegion TSS window for \code{annotatePeak}.
#' @return the \code{csAnno} object from \code{ChIPseeker::annotatePeak}, or NULL
#'   if ChIPseeker/txdb are unavailable.
#' @export
annotateConsensus <- function(gr, txdb, tssRegion = c(-3000, 3000)) {
  if (is.null(txdb) || length(gr) == 0L) return(NULL)
  if (!requireNamespace("ChIPseeker", quietly = TRUE)) return(NULL)
  suppressWarnings(ChIPseeker::annotatePeak(gr, TxDb = txdb, tssRegion = tssRegion,
                                            verbose = FALSE))
}
