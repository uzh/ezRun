###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @template app-template
##' @templateVar method ezMethodTargetedSeqQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Cohort QC for targeted / exome sequencing (WES). Runs Picard
##'   CollectHsMetrics, mosdepth, samtools stats/flagstat/idxstats and optional
##'   CollectInsertSizeMetrics / VerifyBamID2 per sample, then aggregates
##'   everything into a single MultiQC report. This is a modern,
##'   MultiQC-based alternative to \code{EzAppTeqc} (which uses the TEQC package).
##' @section Functions:
EzAppTargetedSeqQC <-
  setRefClass(
    "EzAppTargetedSeqQC",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodTargetedSeqQC
        name <<- "EzAppTargetedSeqQC"
        appDefaults <<- rbind(
          designName = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "capture design: sub-folder name under the target-enrichment design directory, or a full path to a target BED"
          ),
          vendor = ezFrame(
            Type = "character",
            DefaultValue = "auto",
            Description = "auto | illumina | agilent. Illumina uses one BED as bait+target; Agilent uses *_Covered.bed (bait) and *_Regions.bed (target)"
          ),
          mapq = ezFrame(
            Type = "integer",
            DefaultValue = 20,
            Description = "minimum mapping quality"
          ),
          baseq = ezFrame(
            Type = "integer",
            DefaultValue = 20,
            Description = "minimum base quality"
          ),
          coverageThresholds = ezFrame(
            Type = "character",
            DefaultValue = "1,10,20,30,50,100",
            Description = "comma-separated coverage depth thresholds for mosdepth"
          ),
          restrictToMainChromosomes = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "restrict QC to assembled chromosomes (drop scaffolds/decoys)"
          ),
          mainChromMaxNchar = ezFrame(
            Type = "integer",
            DefaultValue = 6,
            Description = "a contig is a 'main chromosome' if its name is at most this many characters (chr1..chr22,chrX,chrY,chrM). Scaffolds like GL000009.2 are longer and excluded"
          ),
          runInsertSize = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "run Picard CollectInsertSizeMetrics"
          ),
          markDuplicates = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "run Picard MarkDuplicates first (only if the BAMs are NOT already duplicate-marked)"
          ),
          verifyBamIDSVD = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "path to a VerifyBamID2 SVD resource prefix to enable contamination (FREEMIX); empty to skip"
          )
        )
      }
    )
  )

ezMethodTargetedSeqQC <- function(input = NA, output = NA, param = NA,
                                  htmlFile = "00index.html") {
  library(stringr)

  ## ---- reference ---------------------------------------------------------
  refFasta <- param$ezRef["refFastaFile"]
  if (length(refFasta) == 0 || is.na(refFasta) || !file.exists(refFasta)) {
    ## fallback for direct/testing invocation
    refFasta <- param$fastaFile
  }
  if (is.null(refFasta) || !file.exists(refFasta)) {
    stop("reference fasta not found; set a valid refBuild (or param$fastaFile)")
  }
  refDict <- sub("\\.(fa|fasta|fna)$", ".dict", refFasta)
  refFai <- paste0(refFasta, ".fai")

  ## ---- capture design (bait vs target) -----------------------------------
  beds <- .resolveDesignBeds(param$designName, param$vendor)
  targetBed <- beds$target
  baitBed <- beds$bait
  ezWrite("Targeted-Seq QC: vendor=", beds$vendor,
          " target=", targetBed, " bait=", baitBed)

  ## ---- output working directory ------------------------------------------
  setwdNew(basename(output$getColumn("Report")))
  reportDir <- getwd()
  perSampleDir <- "per_sample"
  intervalsDir <- "intervals"
  dir.create(perSampleDir, showWarnings = FALSE)
  dir.create(intervalsDir, showWarnings = FALSE)

  cores <- as.integer(param$cores)
  covThresholds <- gsub("\\s", "", param$coverageThresholds)

  ## ---- prepare reference indexes if missing ------------------------------
  if (!file.exists(refFai)) {
    ezSystem(paste("samtools faidx", refFasta))
  }
  if (!file.exists(refDict)) {
    ezSystem(paste(prepareJavaTools("picard"), "CreateSequenceDictionary",
                   paste0("R=", refFasta), paste0("O=", refDict)))
  }

  ## ---- interval lists (once) ---------------------------------------------
  targetIntervals <- file.path(intervalsDir, "target.interval_list")
  ezSystem(paste(prepareJavaTools("picard"), "BedToIntervalList",
                 paste0("I=", targetBed), paste0("O=", targetIntervals),
                 paste0("SD=", refDict)))
  if (identical(baitBed, targetBed)) {
    baitIntervals <- targetIntervals
  } else {
    baitIntervals <- file.path(intervalsDir, "bait.interval_list")
    ezSystem(paste(prepareJavaTools("picard"), "BedToIntervalList",
                   paste0("I=", baitBed), paste0("O=", baitIntervals),
                   paste0("SD=", refDict)))
  }

  ## ---- main-chromosome restriction ---------------------------------------
  samples <- input$getNames()
  bamFiles <- input$getFullPaths("BAM")
  mainChroms <- NULL
  if (isTRUE(as.logical(param$restrictToMainChromosomes))) {
    allContigs <- .bamContigs(bamFiles[1])
    mainChroms <- allContigs[nchar(allContigs) <= as.integer(param$mainChromMaxNchar)]
    if (length(mainChroms) == 0) {
      stop("restrictToMainChromosomes: no contig name is <= ",
           param$mainChromMaxNchar, " characters; adjust mainChromMaxNchar")
    }
    ezWrite("Restricting QC to ", length(mainChroms), " main chromosome(s); ",
            length(allContigs) - length(mainChroms), " scaffold(s) excluded")
  }

  ## ---- per-sample QC (parallel over BAMs) --------------------------------
  ## threads per tool so we do not oversubscribe cores
  threadsPerTool <- max(1L, as.integer(floor(cores / max(1L, length(bamFiles)))))
  threadsPerTool <- max(threadsPerTool, 1L)
  ezMclapply(seq_along(bamFiles), function(i) {
    .runSampleTargetedQC(
      bamFile = bamFiles[i],
      sampleName = samples[i],
      outDir = file.path(perSampleDir, samples[i]),
      refFasta = refFasta,
      baitIntervals = baitIntervals,
      targetIntervals = targetIntervals,
      targetBed = targetBed,
      covThresholds = covThresholds,
      mapq = as.integer(param$mapq),
      baseq = as.integer(param$baseq),
      mainChroms = mainChroms,
      runInsertSize = isTRUE(as.logical(param$runInsertSize)),
      markDuplicates = isTRUE(as.logical(param$markDuplicates)),
      verifyBamIDSVD = param$verifyBamIDSVD,
      threads = threadsPerTool
    )
  }, mc.cores = min(cores, length(bamFiles)), mc.preschedule = FALSE)

  ## ---- MultiQC aggregation -----------------------------------------------
  ## Point TMPDIR at the working dir: on some shared filesystems MultiQC's
  ## end-of-run temp-dir cleanup hits read-only template files and burns
  ## minutes in its quadratic-backoff retry loop before giving up.
  multiqcDir <- "multiqc"
  dir.create(multiqcDir, showWarnings = FALSE)
  multiqcTmp <- file.path(getwd(), ".multiqc_tmp")
  dir.create(multiqcTmp, showWarnings = FALSE)
  multiqcCfg <- system.file("templates", "TargetedSeqQC_multiqc_config.yaml",
                            package = "ezRun")
  cfgArg <- if (nzchar(multiqcCfg) && file.exists(multiqcCfg)) {
    paste("--config", multiqcCfg)
  } else {
    ""
  }
  ezSystem(paste0(
    'TMPDIR="', multiqcTmp, '" multiqc "', perSampleDir, '"',
    ' -o "', multiqcDir, '" -n multiqc_report -f ', cfgArg
  ))
  multiqcHtml <- file.path(multiqcDir, "multiqc_report.html")
  if (!file.exists(multiqcHtml)) {
    stop("MultiQC did not produce a report at ", multiqcHtml)
  }

  ## ---- cohort summary + FGCZ report --------------------------------------
  qcSummary <- .collectHsMetricsSummary(perSampleDir, samples)
  ezWrite.table(qcSummary, file = "hs_metrics_summary.tsv", row.names = FALSE)

  makeRmdReport(
    param = param,
    dataset = input$meta,
    qcSummary = qcSummary,
    multiqcRelPath = multiqcHtml,
    rmdFile = "TargetedSeqQC.Rmd",
    reportTitle = paste("Targeted-Seq QC:", param$name)
  )

  return("Success")
}

##' @describeIn ezMethodTargetedSeqQC Resolve bait/target BEDs from a design
##'   name or path following the vendor convention.
.resolveDesignBeds <- function(designName, vendor = "auto") {
  vendor <- match.arg(tolower(vendor), c("auto", "illumina", "agilent"))
  ## a full path to a single target BED
  if (file.exists(designName) && grepl("\\.bed$", designName, ignore.case = TRUE)) {
    if (vendor == "agilent") {
      stop("vendor='agilent' needs a design folder with *_Covered.bed and *_Regions.bed, not a single BED")
    }
    return(list(vendor = "illumina", target = designName, bait = designName))
  }
  ## otherwise a sub-folder under the design directory
  designDir <- if (dir.exists(designName)) {
    designName
  } else {
    file.path(TEQC_DESIGN_DIR, designName)
  }
  if (!dir.exists(designDir)) {
    stop("capture design not found: ", designName,
         " (looked in ", TEQC_DESIGN_DIR, ")")
  }
  covered <- list.files(designDir, pattern = "Covered\\.bed$", full.names = TRUE)
  regions <- list.files(designDir, pattern = "Regions\\.bed$", full.names = TRUE)
  allBeds <- list.files(designDir, pattern = "\\.bed$", full.names = TRUE)
  hasPair <- length(covered) >= 1 && length(regions) >= 1
  if (vendor == "agilent" || (vendor == "auto" && hasPair)) {
    if (!hasPair) {
      stop("Agilent design requires both *_Covered.bed (bait) and *_Regions.bed (target) in ", designDir)
    }
    return(list(vendor = "agilent", target = regions[1], bait = covered[1]))
  }
  ## Illumina / single-BED: prefer a Covered/Regions file, else the only BED
  single <- if (length(regions) >= 1) {
    regions[1]
  } else if (length(covered) >= 1) {
    covered[1]
  } else if (length(allBeds) == 1) {
    allBeds[1]
  } else if (length(allBeds) >= 1) {
    allBeds[1]
  } else {
    stop("no .bed file found in ", designDir)
  }
  list(vendor = "illumina", target = single, bait = single)
}

##' @describeIn ezMethodTargetedSeqQC Read @SQ contig names from a BAM header.
.bamContigs <- function(bamFile) {
  hdr <- system(paste("samtools view -H", shQuote(bamFile)), intern = TRUE)
  sq <- grep("^@SQ", hdr, value = TRUE)
  sub(".*\tSN:([^\t]+).*", "\\1", sq)
}

##' @describeIn ezMethodTargetedSeqQC Run all per-sample QC steps for one BAM.
.runSampleTargetedQC <- function(bamFile, sampleName, outDir, refFasta,
                                 baitIntervals, targetIntervals, targetBed,
                                 covThresholds, mapq, baseq, mainChroms,
                                 runInsertSize, markDuplicates, verifyBamIDSVD,
                                 threads) {
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(outDir, sampleName)

  ## index if needed
  if (!any(file.exists(c(paste0(bamFile, ".bai"),
                         sub("\\.bam$", ".bai", bamFile))))) {
    ezSystem(paste("samtools index -@", threads, shQuote(bamFile)))
  }

  ## optional MarkDuplicates
  useBam <- bamFile
  if (isTRUE(markDuplicates)) {
    useBam <- paste0(prefix, ".markdup.bam")
    ezSystem(paste(prepareJavaTools("picard"), "MarkDuplicates",
                   paste0("INPUT=", bamFile), paste0("OUTPUT=", useBam),
                   paste0("METRICS_FILE=", prefix, ".dup_metrics.txt"),
                   "CREATE_INDEX=true", "ASSUME_SORT_ORDER=coordinate"))
  }

  ## Picard CollectHsMetrics (core)
  ezSystem(paste(prepareJavaTools("picard"), "CollectHsMetrics",
                 paste0("INPUT=", useBam),
                 paste0("OUTPUT=", prefix, ".hs_metrics.txt"),
                 paste0("REFERENCE_SEQUENCE=", refFasta),
                 paste0("BAIT_INTERVALS=", baitIntervals),
                 paste0("TARGET_INTERVALS=", targetIntervals),
                 paste0("PER_TARGET_COVERAGE=", prefix, ".per_target.txt"),
                 paste0("MINIMUM_MAPPING_QUALITY=", mapq),
                 paste0("MINIMUM_BASE_QUALITY=", baseq)))

  ## mosdepth (uses the target BED directly)
  ezSystem(paste("mosdepth --by", shQuote(targetBed),
                 "--no-per-base --thresholds", covThresholds,
                 "--mapq", mapq, "-t", threads,
                 shQuote(prefix), shQuote(useBam)))

  ## samtools stats/flagstat/idxstats, restricted to main chromosomes if set
  regionArg <- if (!is.null(mainChroms)) paste(mainChroms, collapse = " ") else ""
  ezSystem(paste("samtools stats -@", threads, shQuote(useBam), regionArg,
                 ">", paste0(shQuote(prefix), ".samtools_stats.txt")))
  if (is.null(mainChroms)) {
    ezSystem(paste("samtools flagstat -@", threads, shQuote(useBam),
                   ">", paste0(shQuote(prefix), ".flagstat.txt")))
  } else {
    ## flagstat has no region option -> stream a region-limited view. The pipe
    ## goes into a small shell script because ezSystem() refuses a command that
    ## contains both a '|' and single quotes; running `bash <script>` sidesteps
    ## that (no pipe in the ezSystem command itself).
    flagstatScript <- paste0(prefix, ".flagstat.sh")
    writeLines(c(
      "#!/bin/bash",
      "set -euo pipefail",
      paste("samtools view -u -@", threads, shQuote(useBam), regionArg,
            "| samtools flagstat -@", threads, "-",
            ">", paste0(shQuote(prefix), ".flagstat.txt"))
    ), flagstatScript)
    ezSystem(paste("bash", shQuote(flagstatScript)))
  }
  idxOut <- paste0(prefix, ".idxstats.txt")
  ezSystem(paste("samtools idxstats", shQuote(useBam), ">", shQuote(idxOut)))
  if (!is.null(mainChroms)) {
    lines <- readLines(idxOut, warn = FALSE)
    keep <- sub("\t.*", "", lines) %in% mainChroms
    writeLines(lines[keep], idxOut)
  }

  ## optional insert-size metrics
  if (isTRUE(runInsertSize)) {
    ezSystem(paste(prepareJavaTools("picard"), "CollectInsertSizeMetrics",
                   paste0("I=", useBam),
                   paste0("O=", prefix, ".insert_size_metrics.txt"),
                   paste0("H=", prefix, ".insert_size_histogram.pdf")))
  }

  ## optional contamination (VerifyBamID2)
  if (!is.null(verifyBamIDSVD) && nzchar(verifyBamIDSVD) && Sys.which("VerifyBamID") != "") {
    ezSystem(paste("VerifyBamID", "--SVDPrefix", shQuote(verifyBamIDSVD),
                   "--Reference", shQuote(refFasta),
                   "--BamFile", shQuote(useBam),
                   "--Output", paste0(shQuote(prefix), ".verifybamid")))
  }

  return(outDir)
}

##' @describeIn ezMethodTargetedSeqQC Collect the key WES metrics per sample
##'   from the Picard HsMetrics files into one data frame.
.collectHsMetricsSummary <- function(perSampleDir, samples) {
  keep <- c("PCT_SELECTED_BASES", "MEAN_TARGET_COVERAGE", "MEDIAN_TARGET_COVERAGE",
            "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_50X",
            "FOLD_80_BASE_PENALTY", "ZERO_CVG_TARGETS_PCT",
            "PCT_DUPLICATION", "GC_DROPOUT", "AT_DROPOUT")
  rows <- lapply(samples, function(sm) {
    f <- file.path(perSampleDir, sm, paste0(sm, ".hs_metrics.txt"))
    if (!file.exists(f)) {
      return(NULL)
    }
    m <- ezRead.table(f, comment.char = "#", row.names = NULL, nrows = 1)
    present <- intersect(keep, colnames(m))
    data.frame(Sample = sm, m[, present, drop = FALSE],
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  do.call(rbind, rows)
}
