###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodNfCoreAtacSeq <- function(input = NA, output = NA, param = NA) {
  outFolder = output$getColumn("ATAC_Result") |> basename()

  nfSampleFile <- file.path('dataset.csv')
  nfSampleInfo = getAtacSampleSheet(input, param)
  if (isTRUE(param$qcMode)) {
    nfSampleInfo <- subsampleAtacFastqs(nfSampleInfo, input, param)
  }
  write_csv(nfSampleInfo, nfSampleFile)
  prepNFCoreEnv()
  configFile <- writeNextflowLimits(param)
  ezSystem(buildNfCoreAtacCmd(param, nfSampleFile, outFolder, configFile))
  ezSystem(paste('mv', configFile, outFolder))
  if (isTRUE(param$qcMode)) {
    writeAtacQcModeInfo(nfSampleInfo, param, outFolder)
  }

  ## multiple fastq files per library have been merged by the processing (if any)
  ## now we work with the library names and reduce the dataset
  nfSampleInfo$libName <- paste0(
    nfSampleInfo$sample,
    "_REP",
    nfSampleInfo$replicate
  )
  nfSampleInfo <- nfSampleInfo[!duplicated(nfSampleInfo$sid), ]

  writePerSampleCountPeaksFiles(
    nfSampleInfo,
    countDir = paste0(
      outFolder,
      "/bwa/merged_library/macs2/",
      param$peakStyle,
      "_peak/consensus/"
    )
  )
  renameAtacBigwigs(
    nfSampleInfo,
    bigwigDir = file.path(outFolder, "bwa/merged_library/bigwig")
  )

  jsonFile <- writeAtacIgvSession(
    param,
    outFolder,
    jsonFileName = paste0(outFolder, "/igv_session.json"),
    bigwigRelPath = "/bwa/merged_library/bigwig/",
    baseUrl = file.path(PROJECT_BASE_URL, output$getColumn("ATAC_Result"))
  )
  writeNfCoreIgvHtml(
    param,
    jsonFile,
    title = "NfCoreAtacSeq MultiSample Coverage Tracks",
    htmlTemplate = "templates/igvNfCoreTemplate.html",
    htmlFileName = paste0(outFolder, "/igv_session.html")
  )

  dirsToRemove <- c("genome", "trimgalore", "fastqc", "igv")
  cleanupAtacOutFolder(outFolder, dirsToRemove, isTRUE(param$keepBams))

  return("Success")
}

EzAppNfCoreAtacSeq <- setRefClass(
  "EzAppNfCoreAtacSeq",
  contains = "EzApp",
  methods = list(
    initialize = function() {
      "Initializes the application using its specific defaults."
      runMethod <<- ezMethodNfCoreAtacSeq
      name <<- "EzAppNfCoreAtacSeq"
      ## minimum nf-core parameters
      appDefaults <<- rbind(
        peakStyle = ezFrame(
          Type = "character",
          DefaultValue = "broad",
          Description = "Run MACS2 in broadPeak mode, otherwise in narrowPeak mode"
        ),
        varStabilizationMethod = ezFrame(
          Type = "character",
          DefaultValue = "vst",
          Description = "Use rlog transformation or vst (DESeq2)"
        ),
        keepBams = ezFrame(
          Type = "logical",
          DefaultValue = FALSE,
          Description = "Should bam files be stored"
        ),
        pipelineVersion = ezFrame(
          Type = "character",
          DefaultValue = '2.1.2',
          Description = "specify pipeline version"
        ),
        qcMode = ezFrame(
          Type = "logical",
          DefaultValue = FALSE,
          Description = "QC run: process only the first qcReadsPerSample reads (pairs) of each sample"
        ),
        qcReadsPerSample = ezFrame(
          Type = "numeric",
          DefaultValue = 1e7,
          Description = "number of reads (pairs) per sample used in QC mode"
        )
      )
    }
  )
)

##' @description build the nextflow command line for nf-core/atacseq
buildNfCoreAtacCmd <- function(param, nfSampleFile, outFolder, configFile) {
  fullGenomeSize <- param$ezRef@refFastaFile %>%
    Rsamtools::FaFile() %>%
    GenomeInfoDb::seqlengths() %>%
    sum()
  effectiveGenomeSize <- (fullGenomeSize * 0.8) %>% round()
  ## reuse the prebuilt reference index instead of building it in every run
  bwaIndexDir <- file.path(param$ezRef@refBuildDir, "Sequence/BWAIndex")
  hasBwaIndex <- file.exists(file.path(bwaIndexDir, "genome.fa.bwt")) &&
    !file.exists(file.path(bwaIndexDir, "lock"))

  paste(
    "/usr/local/ngseq/src/nextflow/nextflow run nf-core/atacseq",
    ## i/o
    "--input",
    nfSampleFile,
    "--outdir",
    outFolder,
    ## genome files
    "--fasta",
    param$ezRef@refFastaFile,
    "--gtf",
    param$ezRef@refFeatureFile,
    "--gene_bed",
    str_replace(
      param$ezRef@refAnnotationFile,
      basename(param$ezRef@refAnnotationFile),
      'genes.bed'
    ),
    if (hasBwaIndex) paste("--bwa_index", bwaIndexDir) else "",
    ## parameters
    "--macs_gsize",
    sprintf("%.0f", effectiveGenomeSize),
    if (param[['peakStyle']] == 'broad') "" else "--narrow_peak",
    if (param[['varStabilizationMethod']] != 'vst') {
      "--deseq2_vst false"
    } else {
      ""
    },
    ## configuration
    "-work-dir nfatacseq_work",
    "-profile apptainer",
    "-r",
    param$pipelineVersion,
    "-c",
    configFile,
    param$cmdOptions
  )
}

##' @description get an nf-core/atacseq-formatted sample sheet.
##' nf-core uses the group (sample) and replicate columns to merge libraries and
##' replicates; samples without a group are processed as their own group.
getAtacSampleSheet <- function(input, param) {
  sampleNames <- input$getNames()
  if (ezIsSpecified(param$grouping) && input$hasColumn(param$grouping)) {
    groups <- as.character(input$getColumn(param$grouping))
  } else {
    ezLog(
      "grouping column '",
      param$grouping,
      "' not available; every sample is processed as its own group",
      level = "warn"
    )
    groups <- sampleNames
  }
  isMissing <- is.na(groups) | trimws(groups) %in% c("", "NA")
  if (any(isMissing)) {
    ezLog(
      "no ",
      param$grouping,
      " for sample(s) ",
      paste(sampleNames[isMissing], collapse = ", "),
      "; using the sample name as group",
      level = "warn"
    )
    groups[isMissing] <- sampleNames[isMissing]
  }
  ## nf-core requires sample (group) names without whitespace
  groups <- gsub("[^[:alnum:]_.-]", "_", groups)

  listFastq1 <- input$getFullPathsList("Read1")
  if (isTRUE(param$paired)) {
    fastq2 <- unlist(input$getFullPathsList("Read2"))
  } else {
    fastq2 <- ""
  }

  nfSampleInfo <- ezFrame(
    sample = rep(groups, lengths(listFastq1)),
    fastq_1 = unlist(listFastq1),
    fastq_2 = fastq2,
    replicate = rep(ezReplicateNumber(groups), lengths(listFastq1)),
    sid = rep(sampleNames, lengths(listFastq1))
  )
  return(nfSampleInfo)
}

##' @description QC mode: keep only the first nReads reads (pairs) of each
##' sample; multiple fastq files of a sample are concatenated.
subsampleAtacFastqs <- function(
  nfSampleInfo,
  input,
  param,
  outDir = "qc_fastq"
) {
  nReads <- as.numeric(param$qcReadsPerSample)
  stopifnot(length(nReads) == 1, !is.na(nReads), nReads > 0)
  dir.create(outDir, showWarnings = FALSE)
  if (input$hasColumn("Read Count")) {
    readCounts <- input$getColumn("Read Count")
    readCounts <- setNames(as.numeric(readCounts), names(readCounts))
  } else {
    readCounts <- setNames(rep(NA, input$getLength()), input$getNames())
  }
  paired <- isTRUE(param$paired)
  sids <- unique(nfSampleInfo$sid)
  toSubsample <- sids[is.na(readCounts[sids]) | readCounts[sids] > nReads]
  if (length(toSubsample) < length(sids)) {
    ezLog(
      "QC mode: sample(s) with <= ",
      format(nReads, big.mark = ",", scientific = FALSE),
      " reads are used completely: ",
      paste(setdiff(sids, toSubsample), collapse = ", ")
    )
  }

  headFastq <- function(inFiles, outFile) {
    ## head closes the pipe early; ezSystem runs with pipefail, so the SIGPIPE
    ## of the decompressor is ignored here and the result is validated below
    ezSystem(paste(
      "(pigz -dc",
      paste(inFiles, collapse = " "),
      "|| true) | head -n",
      sprintf("%.0f", 4 * nReads),
      "| pigz -p 2 >",
      outFile
    ))
    as.numeric(ezSystem(paste("pigz -dc", outFile, "| wc -l"), intern = TRUE))
  }

  newRows <- parallel::mclapply(
    toSubsample,
    function(sid) {
      rows <- nfSampleInfo[nfSampleInfo$sid == sid, , drop = FALSE]
      r1 <- file.path(outDir, paste0(sid, "_R1.fastq.gz"))
      nLines <- headFastq(rows$fastq_1, r1)
      r2 <- ""
      if (paired) {
        r2 <- file.path(outDir, paste0(sid, "_R2.fastq.gz"))
        nLines2 <- headFastq(rows$fastq_2, r2)
        if (nLines2 != nLines) {
          stop(
            "QC mode: R1 and R2 of ",
            sid,
            " differ after subsampling (",
            nLines,
            " vs ",
            nLines2,
            " lines)"
          )
        }
      }
      if (nLines == 0 || nLines %% 4 != 0) {
        stop("QC mode: invalid subsampled fastq for ", sid, ": ", nLines, " lines")
      }
      rows <- rows[1, , drop = FALSE]
      rows$fastq_1 <- normalizePath(r1)
      rows$fastq_2 <- if (paired) normalizePath(r2) else ""
      rows$qcReads <- nLines / 4
      rows
    },
    mc.cores = max(1, floor(as.numeric(param$cores) / 4)),
    mc.preschedule = FALSE
  )
  isError <- sapply(newRows, inherits, "try-error")
  if (any(isError)) {
    stop(as.character(newRows[[which(isError)[1]]]))
  }
  newRows <- do.call(rbind, newRows)

  keepRows <- nfSampleInfo[!nfSampleInfo$sid %in% toSubsample, , drop = FALSE]
  if (nrow(keepRows) > 0) {
    keepRows$qcReads <- readCounts[keepRows$sid]
  }
  result <- rbind(newRows, keepRows)
  result <- result[order(match(result$sid, nfSampleInfo$sid)), , drop = FALSE]
  for (i in which(!duplicated(result$sid))) {
    ezLog(
      "QC mode: ",
      result$sid[i],
      " uses ",
      format(result$qcReads[i], big.mark = ",", scientific = FALSE),
      " reads"
    )
  }
  result$qcReads <- NULL
  rownames(result) <- NULL
  return(result)
}

##' @description document in the result folder that this is a QC run
writeAtacQcModeInfo <- function(nfSampleInfo, param, outFolder) {
  lines <- c(
    paste(
      "QC run: only the first",
      sprintf("%.0f", as.numeric(param$qcReadsPerSample)),
      "reads (pairs) of each sample were processed."
    ),
    "Samples with fewer reads were processed completely.",
    paste("Samples:", paste(unique(nfSampleInfo$sid), collapse = ", "))
  )
  writeLines(lines, file.path(outFolder, "qc_mode.txt"))
}

writePerSampleCountPeaksFiles <- function(nfSampleInfo, countDir = ".") {
  libColumnNames <- paste0(nfSampleInfo$libName, ".mLb.clN.sorted.bam")
  sampleNames <- nfSampleInfo$sid
  sampleCountFiles <- paste0(countDir, "/", sampleNames, ".txt")
  annoColumnNames <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  x <- data.table::fread(
    file.path(countDir, "consensus_peaks.mLb.clN.featureCounts.txt"),
    data.table = FALSE
  )
  for (i in 1:nrow(nfSampleInfo)) {
    xSel <- x[, c(annoColumnNames, libColumnNames[i])] |>
      dplyr::rename("matchCounts" := !!libColumnNames[i])
    ezWrite.table(xSel, file = sampleCountFiles[i], row.names = FALSE)
  }
  return(sampleCountFiles)
}


##' @description nf-core names the bigwigs by library (<group>_REP<n>); rename
##' them to the sample names that SUSHI links in the grandchild datasets
renameAtacBigwigs <- function(nfSampleInfo, bigwigDir) {
  from <- file.path(bigwigDir, paste0(nfSampleInfo$libName, ".mLb.clN.bigWig"))
  to <- file.path(bigwigDir, paste0(nfSampleInfo$sid, ".bigWig"))
  isMissing <- !file.exists(from)
  if (any(isMissing)) {
    ezLog(
      "bigwig files not found: ",
      paste(basename(from[isMissing]), collapse = ", "),
      level = "warn"
    )
  }
  file.rename(from[!isMissing], to[!isMissing])
  return(invisible(to[!isMissing]))
}

##' @description clean up NfCoreAtacSeq_result directory
cleanupAtacOutFolder <- function(outFolder, dirsToRemove, keepBams = TRUE) {
  if (!keepBams) {
    bamPath <- paste0(outFolder, "/bwa/")
    bamsToDelete <- dir(
      path = bamPath,
      pattern = "*.bam(.bai)?$",
      recursive = TRUE
    )
    file.remove(file.path(bamPath, bamsToDelete))
    ezLog("Deleted .bam and .bam.bai files from the bwa directory.")
  }
  absolutePaths <- paste(outFolder, dirsToRemove, sep = "/")
  unlink(absolutePaths, recursive = TRUE)
  ezLog(paste0("Deleted subdirectory: ", paste(dirsToRemove, collapse = ',')))
}

##' @description write IGV session in json format
writeAtacIgvSession <- function(
  param,
  outFolder,
  jsonFileName,
  bigwigRelPath,
  baseUrl
) {
  refBuildName = param$ezRef@refBuildName
  refUrlBase = file.path(REF_HOST, param$ezRef@refBuild)
  fastaUrl = sub(
    "Annotation.*",
    "Sequence/WholeGenomeFasta/genome.fa",
    refUrlBase
  )
  faiUrl = paste0(fastaUrl, ".fai")

  bigwigPath = file.path(outFolder, bigwigRelPath)
  bigwigFiles <- dir(path = bigwigPath, pattern = "\\.bigWig$")
  trackNames <- sub("(\\.mLb\\.clN)?\\.bigWig$", "", bigwigFiles)
  bigwigTracks <- lapply(seq_along(bigwigFiles), function(i) {
    list(
      id = trackNames[[i]],
      url = paste0(baseUrl, file.path(bigwigRelPath, bigwigFiles[[i]])),
      format = "bigWig",
      name = trackNames[[i]]
    )
  })
  annotationTracks <- list(
    list(
      id = "genes",
      url = file.path(
        REF_HOST,
        param$ezRef@refBuild,
        'Genes/transcripts.only.gtf'
      ),
      format = "gtf",
      type = "annotation",
      name = "genes"
    ),
    list(
      id = "exons",
      url = file.path(REF_HOST, param$ezRef@refBuild, 'Genes/genes.bed'),
      format = "bed",
      type = "annotation",
      name = "exons"
    )
  )
  tracks <- c(list(list(type = "sequence")), bigwigTracks, annotationTracks)
  jsonLines <- list(
    version = "3.5.3",
    showSampleNames = FALSE,
    reference = list(id = refBuildName, fastaUrl = fastaUrl, indexURL = faiUrl),
    tracks = tracks
  )
  jsonFile <- rjson::toJSON(jsonLines, indent = 5, method = "C")
  write(jsonFile, jsonFileName)
  return(jsonFile)
}

##' @description write html wrapper for IGV session
writeNfCoreIgvHtml = function(
  param,
  jsonFile,
  title,
  htmlTemplate,
  htmlFileName
) {
  htmlLines = readLines(system.file(
    htmlTemplate,
    package = "ezRun",
    mustWork = TRUE
  ))
  htmlLines = gsub("TITLE", title, htmlLines)
  htmlLines = gsub("IGV_JSON_CONTENT", jsonFile, htmlLines)
  writeLines(htmlLines, htmlFileName)
}
