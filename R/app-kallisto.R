###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodKallisto = function(input = NA, output = NA, param = NA) {
  sampleName = input$getNames()
  ## gpukallisto's GPU quant path never writes bootstraps (abundance.h5), so
  ## force bootstrap-samples to 0 in GPU mode; the flag is then dropped below.
  if (!is.null(param$gpu) && param$gpu > 0) {
    param$"bootstrap-samples" = 0
    ezLogMethods(
      "kallisto ran in GPU mode, without bootstrap samples."
    )
  }
  ref = getKallistoReference(param)
  refIdx = paste0(ref, ".idx")
  if (!param$paired) {
    if (param$"fragment-length" == 0) {
      param$"fragment-length" = 180
    }
    if (param$sd == 0) {
      param$sd = 50
    }
    ezLogMethods(
      "Single-end reads were quantified assuming a mean fragment length of ",
      param$"fragment-length",
      " and a standard deviation of ",
      param$sd,
      "."
    )
  }

  iftrue <- function(p, yes, no = "") {
    if (!is.null(p) && p) {
      yes
    } else {
      no
    }
  }
  condCharAdd <- function(name, val) {
    if (ezIsSpecified(val)) {
      paste(name, val)
    } else {
      ""
    }
  }
  condNumAdd <- function(name, val) {
    if (!is.null(val) && val != 0) {
      paste(name, val)
    } else {
      ""
    }
  }

  # Specifying all options
  opt = paste(
    "-i",
    refIdx,
    "-o",
    param$outputDir,
    "-t",
    ezThreads(),
    condNumAdd("--bootstrap-samples", param$"bootstrap-samples"),
    condCharAdd("--seed", param$seed),
    iftrue(param$paired, "", "--single"),
    strandOpt = switch(
      param$strandMode,
      "sense" = "--fr-stranded",
      "antisense" = "--rf-stranded",
      "both" = ""
    ),
    condNumAdd("--fragment-length", param$"fragment-length"),
    condNumAdd("--sd", param$sd)
  )
  ezLogMethods(
    "Transcript abundances were quantified with kallisto quant ",
    switch(
      as.character(param$strandMode)[1],
      "sense" = "strand-specifically, with the first read on the transcript strand",
      "antisense" = "strand-specifically, with the first read on the strand opposite to the transcript",
      "both" = "without strand information",
      "with the library strandedness"
    ),
    "."
  )
  if (
    isTRUE(as.logical(param$paired)) &&
      isTRUE(as.numeric(param$"fragment-length") == 0) &&
      isTRUE(as.numeric(param$sd) == 0)
  ) {
    ezLogMethods(
      "For paired-end reads, kallisto estimated the fragment length distribution from the data."
    )
  }

  trimmedInput = ezMethodFastpTrim(input = input, param = param)

  pathFastqFiles = paste(
    trimmedInput$getColumn("Read1"),
    iftrue(param$paired, trimmedInput$getColumn("Read2"))
  )

  pathAbundance.h5 = file.path(param$outputDir, "abundance.h5")
  pathAbundance.tsv = file.path(param$outputDir, "abundance.tsv")
  pathRunInfo = file.path(param$outputDir, "run_info.json")

  cmd = paste(
    "export HDF5_DISABLE_VERSION_CHECK=1 kallisto; kallisto",
    "quant",
    opt,
    pathFastqFiles
    #  capture kallisto.stderr that contains fragment length estimate; that differs from the fastp estimate
    # "2> kallisto.stderr",
    # " > kallisto.stdout"
  )
  ezSystem(cmd)

  ezSystem(paste("mv", pathAbundance.tsv, basename(output$getColumn("Count"))))
  ## In GPU mode kallisto does not produce bootstraps, so the output dataset
  ## has no bootstrappedCount column; only move the h5 when it is expected.
  if (
    !is.null(param$"bootstrap-samples") &&
      param$"bootstrap-samples" > 0 &&
      output$hasColumn("bootstrappedCount")
  ) {
    ezSystem(paste(
      "mv",
      pathAbundance.h5,
      basename(output$getColumn("bootstrappedCount"))
    ))
  }
  ezSystem(paste("mv", pathRunInfo, basename(output$getColumn("runInfo"))))

  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodKallisto(input=NA, output=NA, param=NA)
##' @description Use this reference class to run kallisto
##' @seealso \code{\link{getKallistoReference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
##' @author Roman Briskine
EzAppKallisto <-
  setRefClass(
    "EzAppKallisto",
    contains = "EzApp",
    methods = list(
      citation = function(param = list()) {
        c(
          "Bray, N.L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 34, 525-527 (2016). https://doi.org/10.1038/nbt.3519",
          "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodKallisto
        name <<- "EzAppKallisto"
        appDefaults <<- rbind(
          gpu = ezFrame(
            Type = "integer",
            DefaultValue = 0,
            Description = "Set to 1 to run gpukallisto (CUDA); no bootstrap output (abundance.h5) is produced and bootstrap-samples is forced to 0"
          ),
          "bootstrap-samples" = ezFrame(
            Type = "integer",
            DefaultValue = 10,
            Description = "Number of bootstrap samples"
          ),
          seed = ezFrame(
            Type = "numeric",
            DefaultValue = 42,
            Description = "seed for the bootstrap sampling"
          ),
          "fragment-length" = ezFrame(
            Type = "integer",
            DefaultValue = 0,
            Description = "estimated average fragment length (required for single-end reads but should be set to 0 for paired-end reads)"
          ),
          sd = ezFrame(
            Type = "numeric",
            DefaultValue = 0,
            Description = 'estimated fragment length standard deviation (required for single-end reads but should be set to 0 for paired-end reads)'
          ),
          outputDir = ezFrame(
            Type = "character",
            DefaultValue = ".",
            Description = "Output directory"
          ),
          secondRef = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "Path to fasta file with additional sequences to add to reference"
          )
        )
      }
    )
  )


getKallistoReference <- function(param) {
  ## Get kallisto version for index versioning
  versionOutput <- system2(
    "kallisto",
    args = "version",
    stdout = TRUE,
    stderr = TRUE
  )[1]
  ## Extract version number (e.g., "kallisto, version 0.51.0" -> "v51-0")
  versionMatch <- regmatches(
    versionOutput,
    regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", versionOutput)
  )
  if (length(versionMatch) > 0) {
    versionParts <- strsplit(versionMatch, "\\.")[[1]]
    kallistoVersion <- paste0(
      "_v",
      versionParts[1],
      versionParts[2],
      "-",
      versionParts[3]
    )
  } else {
    kallistoVersion <- ""
  }

  ## default values
  pathTranscripts <- NULL
  if (ezIsSpecified(param$transcriptFasta)) {
    refBase <- file.path(
      getwd(),
      paste0("kallistoIndex", kallistoVersion, "/transcripts")
    )
    # paste0(file_path_sans_ext(param$trinityFasta), "_kallistoIndex/transcripts")
  } else {
    if (ezIsSpecified(param$transcriptTypes)) {
      kallistoBase <- paste(sort(param$transcriptTypes), collapse = "-")
      ## This is a combination of transcript types to use.
    } else {
      kallistoBase <- ""
    }
    refBase <- ifelse(
      param$ezRef["refIndex"] == "",
      sub(
        ".gtf$",
        paste0(
          "_",
          kallistoBase,
          "_kallistoIndex",
          kallistoVersion,
          "/transcripts"
        ),
        param$ezRef["refFeatureFile"]
      ),
      param$ezRef["refIndex"]
    )
  }

  ## update if secondRef exists - build a temporary index
  if (ezIsSpecified(param$secondRef)) {
    stopifnot(file.exists(param$secondRef))
    refBase <- file.path(
      getwd(),
      paste0("Custom_kallistoIndex", kallistoVersion, "/transcripts")
    )
  }

  ## A finished index is recognised by its .idx file only. The .idx file is
  ## moved into place as the very last step of a build, so it is never seen
  ## half-written.
  ## The lock is a directory: dir.create() creates it or fails if it already
  ## exists in one atomic step, so exactly one job can hold it, also across
  ## nodes. The former lock file (file.exists() followed by a write) let jobs
  ## that started at the same moment all pass the check and build in parallel
  ## in the same directory.
  idxFile <- paste0(refBase, ".idx")
  lockDir <- paste0(dirname(refBase), ".lock")
  dir.create(dirname(lockDir), recursive = TRUE, showWarnings = FALSE)
  hasLock <- FALSE
  i <- 0
  while (!file.exists(idxFile)) {
    if (dir.create(lockDir, showWarnings = FALSE)) {
      hasLock <- TRUE
      ## only the job that created the lock ever removes it
      on.exit(unlink(lockDir, recursive = TRUE), add = TRUE)
      writeLines(
        paste0(
          "host=", Sys.info()[["nodename"]],
          " pid=", Sys.getpid(),
          " slurmJob=", Sys.getenv("SLURM_JOB_ID", "none"),
          " started=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        ),
        file.path(lockDir, "owner.txt")
      )
      break
    }
    ## somebody else builds and we wait
    ownerFile <- file.path(lockDir, "owner.txt")
    lockOwner <- if (file.exists(ownerFile)) {
      paste(readLines(ownerFile, warn = FALSE), collapse = " ")
    } else {
      "unknown"
    }
    if (i >= INDEX_BUILD_TIMEOUT) {
      stop(
        "kallisto index building still in progress after ",
        INDEX_BUILD_TIMEOUT, " min; lock: ", lockDir, " (", lockOwner, ")",
        ". If that job no longer runs, remove the lock and resubmit."
      )
    }
    if (i == 0) {
      ezLog("waiting for kallisto index built by another job: ", lockOwner)
    }
    Sys.sleep(60)
    i <- i + 1
  }
  if (hasLock && file.exists(idxFile)) {
    ## another job finished the index between our check and taking the lock
    hasLock <- FALSE
  }
  if (!hasLock) {
    ## the index is built and complete
    if (
      !ezIsSpecified(param$transcriptFasta) &&
        !ezIsSpecified(param$secondRef) &&
        !ezIsSpecified(param$ezRef["refIndex"])
    ) {
      ezLogMethods(
        "A kallisto index of the transcript sequences of the reference annotation",
        if (ezIsSpecified(param$transcriptTypes)) {
          paste0(
            ", restricted to transcripts of type ",
            paste(sort(param$transcriptTypes), collapse = ", ")
          )
        },
        ", was used."
      )
    }
    return(refBase)
  }
  ## we hold the lock and have to build the index
  ## Build in a private directory next to the final one (same filesystem, so
  ## the final file.rename is atomic). This keeps the kallisto temp directory
  ## and any partial files out of the shared reference directory.
  buildDir <- paste0(
    dirname(refBase), ".build_", Sys.info()[["nodename"]], "_", Sys.getpid()
  )
  unlink(buildDir, recursive = TRUE)
  dir.create(buildDir)
  on.exit(unlink(buildDir, recursive = TRUE), add = TRUE, after = FALSE)
  buildBase <- file.path(buildDir, basename(refBase))

  job <- ezJobStart("kallisto index")
  if (ezIsSpecified(param$transcriptFasta)) {
    pathTranscripts <- param$transcriptFasta
  } else {
    pathTranscripts <- paste0(buildBase, ".fa")
    transcripts <- getTranscriptSequences(param)
    if (ezIsSpecified(param$transcriptTypes)) {
      ## Prepare the subset of transcipts sequence
      seqAnno <- ezFeatureAnnotation(
        param$ezRef@refAnnotationFile,
        dataFeatureType = "transcript"
      )
      transcriptsUse <- rownames(seqAnno)[
        seqAnno$type %in% param$transcriptTypes
      ]
      transcriptsUse <- intersect(transcriptsUse, names(transcripts))
      transcripts <- transcripts[unique(transcriptsUse)]
    }

    ## If secondRef is specified, add those sequences
    if (ezIsSpecified(param$secondRef)) {
      require(Biostrings)
      secondRefSeqs <- readDNAStringSet(param$secondRef)
      transcripts <- c(transcripts, secondRefSeqs)
    }

    writeXStringSet(transcripts, pathTranscripts)
  }
  ## kallisto >= 0.50 writes temporary files to ./tmp unless -T is given
  tmpOption <- ""
  if (
    length(versionMatch) > 0 &&
      package_version(versionMatch) >= package_version("0.50.0")
  ) {
    tmpOption <- paste("-T", file.path(buildDir, "tmp"))
  }
  cmd <- sprintf(
    "kallisto index %s -i %s.idx %s",
    tmpOption,
    buildBase,
    pathTranscripts
  )
  ezSystem(cmd)
  ## publish: the .idx file last, it marks the index as complete
  dir.create(dirname(refBase), recursive = TRUE, showWarnings = FALSE)
  if (!ezIsSpecified(param$transcriptFasta)) {
    stopifnot(file.rename(pathTranscripts, paste0(refBase, ".fa")))
  }
  stopifnot(file.rename(paste0(buildBase, ".idx"), idxFile))
  if (ezIsSpecified(param$transcriptFasta)) {
    ezLogMethods(
      "The kallisto index (k-mer length 31) was built from the supplied transcript sequences."
    )
  } else {
    ezLogMethods(
      "The kallisto index (k-mer length 31) was built from the transcript sequences of the reference annotation",
      if (ezIsSpecified(param$transcriptTypes)) {
        paste0(
          ", restricted to transcripts of type ",
          paste(sort(param$transcriptTypes), collapse = ", ")
        )
      },
      if (ezIsSpecified(param$secondRef)) ", with the secondRef sequences appended",
      "."
    )
  }

  ## Run kallisto inspect to create transcripts.info file
  # We skip this for now as it causes memory overflow
  # inspectCmd = sprintf("kallisto inspect %s.idx > %s.info 2>&1", refBase, refBase)
  # ezSystem(inspectCmd)

  ezWriteElapsed(job, "done")
  return(refBase)
}
