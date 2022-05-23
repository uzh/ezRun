### convert a fastq file into bam file, with additional tags;
### Mainly for pooled single cell samples
### fastqI1Fn: cell barcode
### fastqI2Fn: library barcode
fastq2bam <- function(fastqFn, refFn, bamFn, fastqI1Fn = NULL, fastqI2Fn = NULL) {
  require(ShortRead)
  tempSamFn <- str_c(Sys.getpid(), "temp.sam", sep = "-")
  tempBamFn <- str_c(Sys.getpid(), "temp", sep = "-")

  ## SQ header; required by asBam
  seqlengths <- fasta.seqlengths(refFn)
  seqlengthsHeader <- str_c("@SQ\tSN:", names(seqlengths), "\tLN:", seqlengths)

  ## Add headers
  cat("@HD\tVN:1.5\tSO:unknown", file = tempSamFn, sep = "\n", append = FALSE)
  cat(seqlengthsHeader, file = tempSamFn, sep = "\n", append = TRUE)

  ## parse the fastq files
  f <- FastqStreamer(fastqFn, 1e6)
  on.exit(close(f), add = TRUE)
  if (!is.null(fastqI1Fn)) {
    stopifnot(!is.null(fastqI2Fn))
    f1 <- FastqStreamer(fastqI1Fn, 1e6)
    f2 <- FastqStreamer(fastqI2Fn, 1e6)
    on.exit(close(f1), add = TRUE)
    on.exit(close(f2), add = TRUE)
  }
  while (length(fq <- yield(f))) {
    samText <- str_c(id(fq), "4", "*", "0", "0", "*", "*", "0", "0",
      sread(fq), quality(quality(fq)),
      sep = "\t"
    )
    if (!is.null(fastqI1Fn)) {
      fq1 <- yield(f1)
      fq2 <- yield(f2)
      samText <- str_c(
        samText, "\t",
        "BC:Z:", sread(fq2), "\t",
        "QT:Z:", quality(quality(fq2)), "\t",
        "CR:Z:", sread(fq1), "\t",
        "CY:Z:", quality(quality(fq1))
      )
    }
    cat(samText, file = tempSamFn, append = TRUE, sep = "\n")
  }
  tempBamFn <- asBam(
    file = tempSamFn, destination = tempBamFn,
    indexDestination = FALSE
  )
  bamFn <- sortBam(
    file = tempBamFn, destination = bamFn, byQname = TRUE,
    maxMemory = 1024
  )
  on.exit(file.remove(tempSamFn), add = TRUE)
  on.exit(file.remove(tempBamFn), add = TRUE)
  invisible(bamFn)
}

### convert fastq files into a bam file, with read group tags
fastqs2bam <- function(fastqFns, fastq2Fns = NULL, readGroupNames = NULL,
                       bamFn, platform = "illumina", mc.cores = ezThreads()) {
  paired <- FALSE
  if (!is.null(fastq2Fns)) {
    stopifnot(length(fastqFns) == length(fastq2Fns))
    paired <- TRUE
  }
  sampleBasenames <- sub(
    "\\.(fastq|fq)(\\.gz){0,1}$", "",
    basename(fastqFns)
  )
  if (is.null(readGroupNames)) {
    readGroupNames <- sampleBasenames
  } else {
    stopifnot(length(fastqFns) == length(readGroupNames))
  }

  emptyFastqs <- countReadsInFastq(fastqFns) == 0L
  ## tempty fastq files fail in the picard tool of FastqToSam and MergeSamFiles
  ## Convert and merge the non-empty fastqs first and alter the header of merged bam file

  cmd <- str_c(
    prepareJavaTools("picard"), "FastqToSam",
    str_c("F1=", fastqFns),
    sep = " "
  )
  if (isTRUE(paired)) {
    cmd <- str_c(cmd, str_c("F2=", fastq2Fns), sep = " ")
  }
  cmd <- str_c(
    cmd, str_c("O=", sampleBasenames, ".bam"),
    str_c("SAMPLE_NAME=", readGroupNames),
    str_c("LIBRARY_NAME=", readGroupNames),
    str_c("READ_GROUP_NAME=", readGroupNames),
    str_c("RUN_DATE=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S+00:00")),
    str_c("PLATFORM=", platform),
    "SEQUENCING_CENTER=FGCZ",
    sep = " "
  )

  ezMclapply(cmd[!emptyFastqs], ezSystem,
    mc.preschedule = FALSE,
    mc.cores = mc.cores
  )

  ## picard tools merge
  cmd <- str_c(
    prepareJavaTools("picard"), "MergeSamFiles",
    str_c("I=", str_c(sampleBasenames, ".bam")[!emptyFastqs], collapse = " "),
    str_c("O=", bamFn), "USE_THREADING=true", "SORT_ORDER=queryname",
    sep = " "
    )
  ezSystem(cmd)
  file.remove(str_c(sampleBasenames, ".bam"))

  if (any(emptyFastqs)) {
    tempHeaderFn <- tempfile(
      pattern = "nonEmpty", tmpdir = getwd(),
      fileext = ".header"
    )
    cmd <- str_c("samtools view -H", bamFn, ">", tempHeaderFn, sep = " ")
    ezSystem(cmd)

    extraHeaders <- str_c("@RG",
      str_c("ID:", readGroupNames[emptyFastqs]),
      str_c("SM:", readGroupNames[emptyFastqs]),
      str_c("LB:", readGroupNames[emptyFastqs]),
      str_c("PL:", platform),
      "CN:FGCZ",
      str_c("DT:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S+00:00")),
      sep = "\t"
    )
    con <- file(tempHeaderFn, open = "a")
    writeLines(extraHeaders, con = con)
    close(con)

    bamReheaderFn <- tempfile(
      pattern = "reheader", tmpdir = getwd(),
      fileext = ".bam"
    )
    cmd <- str_c("samtools reheader", tempHeaderFn, bamFn, ">", bamReheaderFn,
                 ssep = " ")
    ezSystem(cmd)

    file.remove(bamFn)
    file.rename(from = bamReheaderFn, to = bamFn)
  }

  invisible(bamFn)
}

### Convert  bam/sam files into fastqs
bam2fastq <- function(bamFn, OUTPUT_PER_RG = TRUE, OUTPUT_DIR = ".",
                      paired = FALSE,
                      fastqFns = str_replace(bamFn, "(\\.bam|\\.sam)$", "_R1.fastq"),
                      fastq2Fns = str_replace(fastqFns, "_R1", "_R2.")) {
  if (isTRUE(OUTPUT_PER_RG)) {
    ## When OUTPUT_PER_RG is TRUE, we only process one uBam each time.
    stopifnot(length(bamFn) == 1L)
    ## I don't want to parse the bam header to get RG IDs
    ## Put them in a tempdir and move to OUTPUT_DIR later
    tempDIR <- str_c("SamtoFastqTempDir", Sys.getpid(), sep = "-")
    dir.create(tempDIR)
    on.exit(unlink(tempDIR, recursive = TRUE), add = TRUE)
    cmd <- str_c(
      prepareJavaTools("picard"), "SamToFastq",
      str_c("I=", bamFn),
      str_c("OUTPUT_DIR=", tempDIR),
      "OUTPUT_PER_RG=true RG_TAG=ID",
      sep = " "
    )
    ezSystem(cmd)
    fastqFns <- list.files(path = tempDIR, pattern = "_1\\.fastq$")
    fromFns <- file.path(tempDIR, fastqFns)
    toFns <- file.path(OUTPUT_DIR, str_replace(fastqFns, "_1\\.fastq$", "_R1.fastq"))
    file.rename(from = fromFns, to = toFns)
    if (isTRUE(paired)) {
      file.rename(
        from = str_replace(fromFns, "_1\\.fastq$", "_2.fastq"),
        to = str_replace(toFns, "_R1\\.fastq$", "_R2.fastq")
      )
    }
    return(invisible(toFns))
    ## This is not much slower than splitBambyRG and SamToFastq in parallel
  } else {
    cmd <- str_c(
      prepareJavaTools("picard"), "SamToFastq",
      str_c("I=", bamFn),
      str_c("FASTQ=", fastqFns),
      sep = " "
    )
    if (isTRUE(paired)) {
      cmd <- str_c(cmd, str_c("SECOND_END_FASTQ=", fastq2Fns), sep = " ")
    }
    lapply(cmd, ezSystem)
    return(invisible(fastqFns))
  }
}

.ezMethodBam2Fastq_PerRG <- function(input = NA, output = NA, param = NA) {
  output <- EzDataset(
    file = input$getFullPaths("CellDataset"),
    dataRoot = param$dataRoot
  )
  output$setColumn("Read1", file.path(
    getwd(), str_c(output$getNames(), "_R1.fastq"))
  )
  if (param$paired) {
    output$setColumn("Read2", file.path(
      getwd(), str_c(output$getNames(), "_R2.fastq"))
    )
  } else {
    if ("Read2" %in% input$colNames) {
      output$setColumn("Read2", NULL)
    }
  }
  output$dataRoot <- NULL

  if (all(file.exists(output$getColumn("Read1")))) {
    return(output)
  }
  allFastqFns <- bam2fastq(
    bamFn = input$getFullPaths("Read1"),
    OUTPUT_PER_RG = TRUE, OUTPUT_DIR = getwd(),
    paired = param$paired
  )
  fastqsToRemove <- allFastqFns[!(basename(allFastqFns) %in%
                                    basename(output$getColumn("Read1")))]
  ## When we want only subset of the cells to keep
  file.remove(fastqsToRemove)
  file.remove(str_replace(fastqsToRemove, "_R1\\.fastq$", "_R2.fastq"))
  return(output)
}

.ezMethodBam2Fastq_NotPerRG <- function(input = NA, output = NA, param = NA) {
  if (!is(output, "EzDataset")) {
    outputMeta <- input$meta
    outputMeta[["Read1"]] <- file.path(
      getwd(), str_c(rownames(outputMeta), "_R1.fastq"))
    outputMeta[["Read1 [File]"]] <- NULL
    if (param$paired) {
      outputMeta[["Read2"]] <-  file.path(
        getwd(), str_c(rownames(outputMeta), "_R2.fastq"))
      outputMeta[["Read2 [File]"]] <- NULL
    } else {
      outputMeta[["Read2"]] <- NULL
    }
    output <- EzDataset(meta = outputMeta, dataRoot = NULL)
  }
  if (all(file.exists(output$getColumn("Read1")))) {
    return(output)
  }
  bam2fastq(
    bamFn = input$getFullPaths("Read1"),
    OUTPUT_PER_RG = FALSE, paired = param$paired,
    fastqFns = output$getColumn("Read1")
  )

  output$setColumn(
    "Read Count",
    countReadsInFastq(output$getColumn("Read1"))
  )
  return(output)
}

ezMethodBam2Fastq <- function(input = NA, output = NA, param = NA,
                              OUTPUT_PER_RG = TRUE) {
  if (isTRUE(OUTPUT_PER_RG)) {
    output <- .ezMethodBam2Fastq_PerRG(input = input, output = output,
                                       param = param)
  } else {
    output <- .ezMethodBam2Fastq_NotPerRG(input = input, output = output,
                                          param = param)
  }
  return(output)
}

countReadsInFastq <- function(fastqFiles) {
  nReads <- sapply(fastqFiles, fastq.geometry)[1, ]
  return(nReads)
}

ezMethodSubsampleFastq <- function(input = NA, output = NA, param = NA, n = 1e6) {
  require(ShortRead)
  ## if output is not an EzDataset, set it!
  if (!is(output, "EzDataset")) {
    output <- input$copy()
    output$setColumn("Read1", str_c(
      getwd(), "/", input$getNames(),
      "-subsample_R1.fastq.gz"
    ))
    if (param$paired) {
      output$setColumn("Read2", str_c(
        getwd(), "/", input$getNames(),
        "-subsample_R2.fastq.gz"
      ))
    } else {
      if ("Read2" %in% input$colNames) {
        output$setColumn("Read2", NULL)
      }
    }
    output$dataRoot <- NULL
  }
  output$setColumn("Read Count", ifelse(input$getColumn("Read Count") < n, 
                                        input$getColumn("Read Count"), 
                                        n))

  dataset <- input$meta
  samples <- rownames(dataset)
  lapply(samples, function(sm, input, output, param, n) {
    subsampleFastqFile(input$getFullPaths("Read1")[sm],
                      output$getColumn("Read1")[sm],
                      nReads=n)
    if (param$paired) {
      subsampleFastqFile(input$getFullPaths("Read2")[sm],
                        output$getColumn("Read2")[sm],
                        nReads=n)
    }
  }, input = input, output = output, param = param, n)
  return(output)
}


subsampleFastqFile <- function(inFile, outFile, nReads, seed=123L){
  
  library(ShortRead)
  fqs <- FastqSampler(inFile, n = nReads, ordered = TRUE)
  set.seed(seed)
  ssReads <- yield(fqs)
  close(fqs)
  writeFastq(ssReads, file = outFile)
  ## seqtk implementation
  # cmd = paste("seqtk sample -s 42", inFile, nReads, "| pigz --fast -p1 >", outFile)
  # ezSystem(cmd)
  return(invisible(outFile))
}
