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
      ## Picard (CollectHsMetrics), mosdepth, samtools stats/flagstat/idxstats and MultiQC run for
      ## every job (.runSampleTargetedQC, ezMethodTargetedSeqQC); the optional MarkDuplicates and
      ## CollectInsertSizeMetrics are Picard too. VerifyBamID2 is not cited: a set verifyBamIDSVD
      ## stops the job early, because none of the app's modules provides a VerifyBamID binary.
      citation = function(param = list()) {
        c(
          "Picard Toolkit. Broad Institute. https://broadinstitute.github.io/picard/",
          "Pedersen, B.S. & Quinlan, A.R. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics 34(5), 867-868 (2018). https://doi.org/10.1093/bioinformatics/btx699",
          "Danecek, P. et al. Twelve years of SAMtools and BCFtools. GigaScience 10(2), giab008 (2021). https://doi.org/10.1093/gigascience/giab008",
          "Ewels, P., Magnusson, M., Lundin, S. & Käller, M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32(19), 3047-3048 (2016). https://doi.org/10.1093/bioinformatics/btw354"
        )
      },
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

  ## VerifyBamID2 is not installed on the cluster: fail now instead of after the
  ## per-sample QC (it used to be skipped silently)
  if (ezIsSpecified(param$verifyBamIDSVD) && !nzchar(Sys.which("VerifyBamID"))) {
    stop("verifyBamIDSVD is set but VerifyBamID (VerifyBamID2) is not on the ",
         "PATH; it is not installed yet. Leave verifyBamIDSVD empty to run ",
         "the QC without contamination estimates.")
  }

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

  ## ---- reference indexes ----------------------------------------------------
  ## never write next to the shared reference: if .fai or .dict is missing,
  ## index a symlinked copy of the FASTA in the job's working directory (next
  ## to, not inside, the report folder that is delivered)
  refFasta <- .targetedQcReference(
    refFasta, localDir = file.path(dirname(reportDir), "targetedSeqQC_ref")
  )
  refDict <- sub("\\.(fa|fasta|fna)$", ".dict", refFasta)

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
    ## restrict the design to the same contigs (HsMetrics, mosdepth)
    sameBeds <- identical(baitBed, targetBed)
    targetBed <- .filterBedToContigs(
      targetBed, mainChroms, file.path(intervalsDir, "target.main.bed")
    )
    baitBed <- if (sameBeds) {
      targetBed
    } else {
      .filterBedToContigs(baitBed, mainChroms, file.path(intervalsDir, "bait.main.bed"))
    }
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

##' @describeIn ezMethodTargetedSeqQC Reference FASTA with .fai and .dict. If
##'   either index is missing next to the shared FASTA, the FASTA is symlinked
##'   into \code{localDir} and the missing indexes are created there (existing
##'   ones are symlinked too); the shared reference is never written to.
.targetedQcReference <- function(refFasta, localDir = "ref") {
  refFai <- paste0(refFasta, ".fai")
  refDict <- sub("\\.(fa|fasta|fna)$", ".dict", refFasta)
  if (file.exists(refFai) && file.exists(refDict)) {
    return(refFasta)
  }
  dir.create(localDir, showWarnings = FALSE)
  localFasta <- file.path(normalizePath(localDir), basename(refFasta))
  localFai <- paste0(localFasta, ".fai")
  localDict <- sub("\\.(fa|fasta|fna)$", ".dict", localFasta)
  if (!file.exists(localFasta)) {
    file.symlink(refFasta, localFasta)
  }
  if (file.exists(refFai)) {
    file.symlink(refFai, localFai)
  } else {
    ezLog("No .fai next to ", refFasta, "; indexing a local copy in ", localDir)
    ezSystem(paste("samtools faidx", shQuote(localFasta)))
  }
  if (file.exists(refDict)) {
    file.symlink(refDict, localDict)
  } else {
    ezLog("No .dict next to ", refFasta, "; creating it in ", localDir)
    ezSystem(paste(prepareJavaTools("picard"), "CreateSequenceDictionary",
                   paste0("R=", localFasta), paste0("O=", localDict)))
  }
  localFasta
}

##' @describeIn ezMethodTargetedSeqQC Keep only the BED records on the given
##'   contigs (track/browser/comment lines are dropped).
.filterBedToContigs <- function(bedFile, contigs, outFile) {
  bedLines <- readLines(bedFile, warn = FALSE)
  bedLines <- bedLines[nzchar(bedLines) & !grepl("^(#|track|browser)", bedLines)]
  keep <- sub("\t.*", "", bedLines) %in% contigs
  if (!any(keep)) {
    stop("No region of ", bedFile, " lies on the main chromosomes (",
         paste(head(contigs, 5), collapse = ", "), ", ...); check the contig ",
         "naming of the design vs the BAM, or set restrictToMainChromosomes = false")
  }
  if (any(!keep)) {
    ezLog(sum(!keep), " of ", length(keep), " regions in ", basename(bedFile),
          " are not on the main chromosomes and were dropped")
  }
  writeLines(bedLines[keep], outFile)
  outFile
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

  ## main-chromosome BAM used by every tool below (HsMetrics, mosdepth,
  ## InsertSize, stats, flagstat, idxstats); removed at the end
  mainBam <- NULL
  if (!is.null(mainChroms)) {
    mainBam <- paste0(prefix, ".mainChroms.bam")
    ezSystem(paste("samtools view -b -@", threads, "-o", shQuote(mainBam),
                   shQuote(useBam), paste(mainChroms, collapse = " ")))
    ezSystem(paste("samtools index -@", threads, shQuote(mainBam)))
    useBam <- mainBam
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

  ## samtools stats/flagstat/idxstats (main-chromosome BAM if restricted)
  ezSystem(paste("samtools stats -@", threads, shQuote(useBam),
                 ">", paste0(shQuote(prefix), ".samtools_stats.txt")))
  ezSystem(paste("samtools flagstat -@", threads, shQuote(useBam),
                 ">", paste0(shQuote(prefix), ".flagstat.txt")))
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

  ## optional contamination (VerifyBamID2); ezMethodTargetedSeqQC stops early
  ## when verifyBamIDSVD is set but the binary is missing
  if (!is.null(verifyBamIDSVD) && nzchar(verifyBamIDSVD)) {
    ezSystem(paste("VerifyBamID", "--SVDPrefix", shQuote(verifyBamIDSVD),
                   "--Reference", shQuote(refFasta),
                   "--BamFile", shQuote(useBam),
                   "--Output", paste0(shQuote(prefix), ".verifybamid")))
  }

  if (!is.null(mainBam)) {
    file.remove(intersect(
      c(mainBam, paste0(mainBam, ".bai")),
      list.files(dirname(mainBam), full.names = TRUE)
    ))
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
