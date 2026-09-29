###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## resolve the first present column among candidate names (tags already stripped
## from input$colNames, so pass base names like "CalledPeaks", "BigWigFile").
.resolvePeakColumn <- function(input, candidates) {
  for (cc in candidates) if (input$hasColumn(cc)) return(cc)
  NULL
}

## named per-sample full paths for an optional column; NA where the column is
## absent or a cell is empty. Returns NULL if the column does not exist at all.
.optionalFullPaths <- function(input, colName) {
  if (is.null(colName) || !input$hasColumn(colName)) return(NULL)
  p <- input$getFullPaths(colName)
  p[!nzchar(p) | is.na(p)] <- NA_character_
  p
}

## Seqinfo + genome size from the refBuild FASTA index.
.refSeqinfo <- function(param) {
  fasta <- param$ezRef["refFastaFile"]
  fai <- paste0(fasta, ".fai")
  if (!file.exists(fai)) Rsamtools::indexFa(fasta)
  tab <- utils::read.table(fai, stringsAsFactors = FALSE)
  si <- GenomeInfoDb::Seqinfo(seqnames = as.character(tab$V1),
                              seqlengths = as.integer(tab$V2))
  list(seqinfo = si, genomeSize = sum(as.numeric(tab$V2)))
}

ezMethodChIPSeqPeakComparison <- function(input = NA, output = NA, param = NA,
                                          htmlFile = "00index.html") {
  setwdNew(basename(output$getColumn("Report")))
  sampleNames <- input$getNames()

  ## ---- 1. resolve inputs, decide which optional blocks are possible --------
  peakCol <- .resolvePeakColumn(input, c("Peaks", "CalledPeaks", "BED", "MACS"))
  if (is.null(peakCol)) stop("No peak column found (looked for Peaks/CalledPeaks/BED).")
  bwCol  <- .resolvePeakColumn(input, c("BigWig", "BigWigFile"))
  bamCol <- .resolvePeakColumn(input, c("BAM"))
  ctrlCol <- .resolvePeakColumn(input, c("Control BAM", "ControlBAM"))
  condCol <- .resolvePeakColumn(input, c("Condition"))

  peakFiles <- input$getFullPaths(peakCol); names(peakFiles) <- sampleNames
  bwFiles   <- .optionalFullPaths(input, bwCol);  if (!is.null(bwFiles)) names(bwFiles) <- sampleNames
  bamFiles  <- .optionalFullPaths(input, bamCol); if (!is.null(bamFiles)) names(bamFiles) <- sampleNames
  ctrlFiles <- .optionalFullPaths(input, ctrlCol)
  condition <- if (!is.null(condCol)) input$getColumn(condCol) else NULL

  if (is.null(bwFiles) && is.null(bamFiles))
    stop("Coverage plots need at least one of BigWig / BAM; the dataset has neither.")

  ref <- .refSeqinfo(param)
  seqinfo <- ref$seqinfo; genomeSize <- ref$genomeSize
  cores <- as.integer(param$cores)

  hasBam <- !is.null(bamFiles) && any(!is.na(bamFiles))
  hasBw  <- !is.null(bwFiles)  && any(!is.na(bwFiles))
  hasCtrl <- !is.null(ctrlFiles) && any(!is.na(ctrlFiles))
  nSamples <- length(sampleNames)

  ## ---- 2. read + normalise peaks ------------------------------------------
  rawPeaks <- ezMclapply(seq_along(peakFiles), function(i)
    readPeakFile(peakFiles[[i]], format = param$peakFormat %||% "auto",
                 sheet = param$xlsxSheet %||% 1), mc.cores = cores)
  names(rawPeaks) <- sampleNames

  blacklist <- NULL
  if (isTRUE(param$useBlacklist) && nzchar(param$blacklistFile %||% "")) {
    blacklist <- tryCatch(rtracklayer::import(param$blacklistFile), error = function(e) NULL)
  }
  harm <- harmonizePeaks(rawPeaks, seqinfo, blacklist = blacklist)
  peaks <- harm$peaks; dropLog <- harm$dropLog

  ## ---- 3. QC ---------------------------------------------------------------
  pkStats <- peakStats(peaks, genomeSize)
  alnStats <- if (hasBam) alignmentStats(bamFiles, cores = cores) else NULL
  complexity <- if (hasBam) do.call(rbind, ezMclapply(seq_along(bamFiles), function(i)
    cbind(sample = sampleNames[i], libraryComplexity(bamFiles[[i]])),
    mc.cores = cores)) else NULL
  frip <- if (hasBam) do.call(rbind, ezMclapply(seq_along(bamFiles), function(i)
    cbind(sample = sampleNames[i],
          computeFrip(bamFiles[[i]], peaks[[i]], paired = isTRUE(param$paired))),
    mc.cores = cores)) else NULL
  ccMetrics <- if (hasBam) do.call(rbind, ezMclapply(seq_along(bamFiles), function(i)
    cbind(sample = sampleNames[i], crossCorrelationMetrics(bamFiles[[i]])),
    mc.cores = cores)) else NULL
  bwCor <- if (hasBw) bigwigCorrelation(bwFiles, cores = cores) else NULL

  ## assemble a wide QC table for flagging
  qcWide <- pkStats[, c("sample", "nPeaks")]
  merger <- function(base, add, cols) {
    if (is.null(add)) return(base)
    merge(base, add[, c("sample", cols)], by = "sample", all.x = TRUE)
  }
  qcWide <- merger(qcWide, frip, "frip");            if ("frip" %in% colnames(qcWide)) colnames(qcWide)[colnames(qcWide) == "frip"] <- "FRiP"
  qcWide <- merger(qcWide, complexity, c("NRF", "PBC1", "PBC2"))
  qcWide <- merger(qcWide, ccMetrics, c("NSC", "RSC"))
  qcWide <- merger(qcWide, alnStats, c("mappedFraction", "mapq30Fraction"))
  qcFlags <- flagQcMetrics(qcWide, markType = param$markType %||% "narrow")

  ## ---- 4. consensus + overlap ---------------------------------------------
  jaccard <- pairwiseJaccard(peaks, minOverlapBp = param$minOverlapBp %||% 1L)
  cons <- buildConsensus(peaks, minSamples = param$minSamplesForConsensus %||% 1L,
                         minOverlapBp = param$minOverlapBp %||% 1L)
  consensus <- cons$consensus; occupancy <- cons$occupancyMatrix

  ## ---- 5. quantification + clustering -------------------------------------
  quantFrom <- param$quantifyFrom %||% "bigwig"
  signalMat <- NULL
  if (length(consensus) > 0) {
    signalMat <- if (quantFrom == "bam" && hasBam)
      quantifyRegions(consensus, bamFiles = bamFiles, cores = cores,
                      paired = isTRUE(param$paired))
    else if (hasBw)
      quantifyRegions(consensus, bwFiles = bwFiles, cores = cores)
    else if (hasBam)
      quantifyRegions(consensus, bamFiles = bamFiles, cores = cores,
                      paired = isTRUE(param$paired))
    else NULL
  }
  normMat <- if (!is.null(signalMat))
    tryCatch(normalizeSignal(signalMat, method = param$normalization %||% "CPM"),
             error = function(e) signalMat) else NULL

  ## ---- 6. top-n significant peaks -> per-sample coverage tracks -----------
  ## `topN` is small here (e.g. 3-10): the n strongest distinct binding loci in
  ## the whole dataset, for which the report draws genome-browser coverage tracks
  ## across all samples. Peak characterisation (footprint, cumulative signal) is
  ## computed over the FULL peak sets, not these n loci.
  topN <- as.integer(param$topN %||% 5L)
  rankBy <- param$rankBy %||% "signalValue"
  topLoci <- topSignificantPeaks(peaks, n = topN, rankBy = rankBy)
  locusCov <- if (hasBw)
    extractLocusCoverage(topLoci, bwFiles, flank = param$profileExtend %||% 2000L,
                         binSize = param$profileBinSize %||% 50L) else list()
  footprint <- genomeFootprint(peaks, seqinfo)
  cumSignal <- lapply(peaks, cumulativeSignalCurve, rankBy = rankBy)

  ## consensus annotation (best-effort)
  txdb <- tryCatch(GenomicFeatures::makeTxDbFromGFF(param$ezRef["refFeatureFile"]),
                   error = function(e) NULL)
  consAnno <- annotateConsensus(consensus, txdb)

  ## gene models for the coverage-track panels (so the report can show whether a
  ## gene sits near each top peak). Stored as a light GRanges with a symbol mcol.
  genes <- tryCatch({
    imp <- rtracklayer::import(param$ezRef["refFeatureFile"])
    g <- if ("type" %in% colnames(S4Vectors::mcols(imp))) imp[imp$type == "gene"] else imp
    if (length(g) == 0L) g <- imp
    sym <- if (!is.null(g$gene_name)) as.character(g$gene_name) else as.character(g$gene_id)
    sym[is.na(sym) | sym == ""] <- as.character(g$gene_id)[is.na(sym) | sym == ""]
    S4Vectors::mcols(g) <- S4Vectors::DataFrame(gene_id = as.character(g$gene_id), symbol = sym)
    g
  }, error = function(e) NULL)

  ## ---- 7. write result files ----------------------------------------------
  if (length(consensus) > 0) {
    rtracklayer::export(consensus, "consensus_peaks.bed", format = "BED")
    occDf <- data.frame(peak = S4Vectors::mcols(consensus)$name,
                        as.data.frame(occupancy),
                        occupancy = S4Vectors::mcols(consensus)$occupancy,
                        check.names = FALSE)
    ezWrite.table(occDf, "occupancy_matrix.tsv", row.names = FALSE)
  }
  if (!is.null(signalMat)) ezWrite.table(signalMat, "consensus_signal_counts.tsv")
  if (!is.null(normMat))   ezWrite.table(normMat, "consensus_signal_normalized.tsv")
  if (length(topLoci) > 0)
    rtracklayer::export(topLoci, "top_significant_peaks.bed", format = "BED")
  writexl::write_xlsx(list(peakStats = pkStats, qc = qcWide, qcFlags = qcFlags,
                           dropLog = dropLog), "qc_metrics.xlsx")

  ## ---- 8. report -----------------------------------------------------------
  chipData <- list(
    param = param, sampleNames = sampleNames, condition = condition,
    peakCol = peakCol, hasBam = hasBam, hasBw = hasBw, hasCtrl = hasCtrl,
    nSamples = nSamples, seqinfo = seqinfo, genomeSize = genomeSize,
    dropLog = dropLog, peakStats = pkStats, alnStats = alnStats,
    complexity = complexity, frip = frip, ccMetrics = ccMetrics,
    bwCor = bwCor, qcWide = qcWide, qcFlags = qcFlags,
    jaccard = jaccard, consensus = consensus, occupancy = occupancy,
    consAnno = consAnno, signalMat = signalMat, normMat = normMat,
    peaks = peaks, topLoci = topLoci, locusCov = locusCov, genes = genes,
    footprint = footprint, cumSignal = cumSignal)

  ## Quarto report via fgczQuartoTemplate (fgcz_render).
  makeQuartoReport(chipData = chipData, qmdFile = "ChIPSeqPeakComparison.qmd",
                   reportTitle = param$name %||% "ChIPSeqPeakComparison")
  return("Success")
}

## NB: the `%||%` null-coalescing operator is already defined at package level
## (see app-diffShot.R / app-ScMultiOmics.R); appDefaults guarantees the params
## referenced with it are present, so the NULL-only semantics are sufficient.

##' @template app-template
##' @templateVar method ezMethodChIPSeqPeakComparison(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Compares annotated ChIP-seq peaks across samples: ENCODE-style
##'   QC, pairwise overlap and consensus peaks, quantitative signal comparison,
##'   and the signal / genomic footprint of the top-N peaks per sample. Expects
##'   already-called peaks (BED / narrowPeak / MACS xls / xlsx) plus BAM and/or
##'   BigWig. Read-only with respect to the input dataset.
##' @seealso \code{\link{readPeakFile}}, \code{\link{buildConsensus}},
##'   \code{\link{profileMatrix}}
EzAppChIPSeqPeakComparison <-
  setRefClass("EzAppChIPSeqPeakComparison",
    contains = "EzApp",
    methods = list(
      methods_facts = function(param = list()) {
        known <- length(param) > 0
        ## BAM / BigWig availability as the code resolves it (first present column, .usableFile);
        ## NA when the input dataset is unknown
        root <- if (ezIsSpecified(param$dataRoot)) param$dataRoot else "/srv/gstore/projects"
        given <- function(cols) {
          if (is.null(attr(param, "input"))) return(NA)
          for (cc in cols) {
            v <- methodsInput(param, cc)
            if (!is.null(v)) return(any(file.exists(file.path(root, v[!is.na(v) & nzchar(v)]))))
          }
          FALSE
        }
        bam <- given("BAM"); bw <- given(c("BigWig", "BigWigFile"))
        bamWhen <- if (is.na(bam)) "When BAM files were given, " else ""
        c(
          ## ezMethodChIPSeqPeakComparison: .resolvePeakColumn(Peaks/CalledPeaks/BED/MACS); markType only in flagQcMetrics
          paste0("The app called no peaks: it compared the peaks already called upstream, read from the first of the input columns Peaks, CalledPeaks, BED or MACS, and markType",
                 if (known && ezIsSpecified(param$markType)) paste0(" (", param$markType, ")"), " only chose the thresholds of the QC flags."),
          ## harmonizePeaks: off-reference / zero-width / out-of-bounds drops, no keepStandardChromosomes
          "Before any statistic, peaks on contigs absent from the reference genome index, of zero width or extending beyond a contig end were dropped and counted per reason; all reference contigs were kept, with no restriction to standard chromosomes.",
          ## ezMethodChIPSeqPeakComparison: blacklist only when useBlacklist and blacklistFile set; import error -> NULL
          if (known) {
            if (isTRUE(as.logical(param$useBlacklist)) && nzchar(param$blacklistFile %||% ""))
              "Peaks overlapping the regions of blacklistFile were removed (useBlacklist true); had that file not been readable, no peaks would have been removed."
            else paste0("No blacklist filtering was applied (",
                        if (isTRUE(as.logical(param$useBlacklist))) "no blacklistFile" else "useBlacklist false", ").")
          },
          ## computeFrip (idxstatsBam total, countBam over reduced peaks, `paired` unused); libraryComplexity (streams the whole BAM)
          if (!isFALSE(bam)) paste0(bamWhen, "FRiP was the number of BAM records overlapping the merged peaks of the sample over all mapped records in the BAM index (read pairs were not collapsed, whatever paired says), and NRF, PBC1 and PBC2 were computed from the start positions (contig, position, strand) of all mapped records of the BAM."),
          ## crossCorrelationMetrics(binSize 10, maxShift 1000, nSub 3e5), alignmentStats(nSub 1e5): first yieldSize records
          if (!isFALSE(bam)) paste0(bamWhen, "NSC and RSC came from the app's own strand cross-correlation (not phantompeakqualtools) on the first 300,000 records of each BAM in 10 bp bins with shifts up to 1000 bp, and the MAPQ 30 fraction from the first 100,000 records; these are the first records of the file, not a random subsample."),
          ## flagQcMetrics / qcThresholdTable(markType): every thresholded metric comes from the BAM QC
          if (known && ezIsSpecified(param$markType) && !isFALSE(bam)) {
            thr <- tryCatch(qcThresholdTable(param$markType), error = function(e) NULL)
            if (!is.null(thr)) thr <- thr[!is.na(thr$min), , drop = FALSE]
            if (!is.null(thr) && nrow(thr) > 0)
              paste0(bamWhen, "QC metrics were flagged good, acceptable or low against fixed ", param$markType, "-mark thresholds (",
                     paste0(thr$metric, " acceptable at ", thr$min, ", good at ", thr$good, collapse = "; "),
                     "); a metric without a threshold for this mark type was not flagged.")
          },
          ## bigwigCorrelation(binSize = 1000L); the report shows the Spearman matrix
          if (!isFALSE(bw)) paste0(if (is.na(bw)) "When BigWig tracks were given, " else "",
                                   "the genome-wide sample correlation is the Spearman correlation of the mean BigWig score in 1000 bp genome bins."),
          ## pairwiseJaccard (reduce, bp intersect/union; minOverlapBp unused); buildConsensus (reduce pooled, minoverlap, minSamples)
          paste0("Pairwise overlap is the base-pair Jaccard index of the merged peaks of two samples (minOverlapBp is not used for it); the consensus is the merge (strand ignored) of all samples' peaks, keeping regions overlapped by at least ",
                 if (known && length(param$minSamplesForConsensus) == 1 && length(param$minOverlapBp) == 1) paste0(param$minSamplesForConsensus, " sample(s) by at least ", param$minOverlapBp, " bp") else "minSamplesForConsensus samples by at least minOverlapBp bp",
                 ", so with 1 sample it is the union of all peaks."),
          ## ezMethodChIPSeqPeakComparison section 5: quantifyRegions (.scoreRegionsBw viewMeans / summarizeOverlaps Union), normalizeSignal
          if (known && ezIsSpecified(param$quantifyFrom)) {
            counts <- "read counts from the BAM files (summarizeOverlaps, Union mode, strand ignored)"
            means <- "the mean BigWig score"
            src <- if (isTRUE(param$quantifyFrom == "bam")) {
              if (isTRUE(bam)) counts else if (isFALSE(bam)) means else paste0(counts, ", or ", means, " when no BAM was given")
            } else {
              if (isTRUE(bw)) means else if (isFALSE(bw)) counts else paste0(means, ", or ", counts, " when no BigWig was given")
            }
            norm <- switch(as.character(param$normalization %||% ""),
              CPM = "divided by the sample's total over the consensus regions and multiplied by one million (CPM, applied the same way to BigWig mean scores)",
              quantile = "quantile-normalised with limma normalizeQuantiles",
              none = "left unnormalised", NULL)
            if (!is.null(norm)) paste0("The signal per consensus region was ", src, " (quantifyFrom ", param$quantifyFrom,
                                       "), and each sample's signal was ", norm, " (normalization ", param$normalization,
                                       "); had the normalisation failed, the unnormalised matrix would have been used.")
          },
          ## ezMethodChIPSeqPeakComparison section 6: topSignificantPeaks (pooled, greedy non-overlapping); footprint/cumSignal on all peaks
          paste0("topN", if (known && length(param$topN) == 1) paste0(" (", param$topN, ")"), " only picked the loci drawn as per-sample coverage tracks: the highest-",
                 if (known && ezIsSpecified(param$rankBy)) param$rankBy else "rankBy", " peaks pooled over all samples, skipping any that overlaps an already chosen one; every statistic, the consensus and the genomic footprint use all peaks, and rankBy also orders the cumulative-signal curve."),
          ## annotateConsensus(tssRegion = c(-3000, 3000)), tryCatch(makeTxDbFromGFF) -> NULL
          "Consensus regions were annotated with ChIPseeker annotatePeak (TSS region -3000 to +3000 bp) only if a TxDb could be built from the reference GTF; otherwise the report has no consensus annotation and no error was raised.",
          ## runDifferentialBinding is not read; ctrlFiles only feeds hasCtrl
          "No differential binding test was run (runDifferentialBinding is not read by the code), and control (input) BAM files, if listed, were not used in any metric."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodChIPSeqPeakComparison
        name <<- "EzAppChIPSeqPeakComparison"
        appDefaults <<- rbind(
          markType = ezFrame(Type = "character", DefaultValue = "narrow",
                             Description = "narrow (TF) or broad (histone); sets QC thresholds"),
          topN = ezFrame(Type = "numeric", DefaultValue = 5,
                         Description = "number of most-significant peaks to show per-sample coverage tracks for"),
          rankBy = ezFrame(Type = "character", DefaultValue = "signalValue",
                           Description = "peak ranking metric (signalValue = fold enrichment / log2 over background)"),
          minSamplesForConsensus = ezFrame(Type = "numeric", DefaultValue = 1,
                           Description = "min samples occupying a region to keep it in the consensus"),
          minOverlapBp = ezFrame(Type = "numeric", DefaultValue = 1,
                           Description = "min bp overlap to count as occupancy"),
          profileExtend = ezFrame(Type = "numeric", DefaultValue = 2000,
                           Description = "bp extended each side of peak centre for profiles"),
          profileBinSize = ezFrame(Type = "numeric", DefaultValue = 50,
                           Description = "profile bin size in bp"),
          quantifyFrom = ezFrame(Type = "character", DefaultValue = "bigwig",
                           Description = "signal source for the consensus matrix: bigwig or bam"),
          normalization = ezFrame(Type = "character", DefaultValue = "CPM",
                           Description = "CPM, none or quantile"),
          peakFormat = ezFrame(Type = "character", DefaultValue = "auto",
                           Description = "peak file format or auto"),
          useBlacklist = ezFrame(Type = "logical", DefaultValue = TRUE,
                           Description = "drop peaks overlapping the blacklist"),
          blacklistFile = ezFrame(Type = "character", DefaultValue = "",
                           Description = "blacklist BED path; empty = none"),
          runDifferentialBinding = ezFrame(Type = "logical", DefaultValue = FALSE,
                           Description = "stub, off; consensus matrix is the input for a follow-up app")
        )
      }
    )
  )
