### convert a fastq file into bam file, with additional tags;
### Mainly for pooled single cell samples
### fastqI1Fn: cell barcode
### fastqI2Fn: library barcode
fastq2bam <- function(fastqFn, refFn, bamFn, fastqI1Fn=NULL, fastqI2Fn=NULL){
  require(ShortRead)
  require(Biostrings)
  tempSamFn <- paste(Sys.getpid(), "temp.sam", sep="-")
  tempBamFn <- paste(Sys.getpid(), "temp", sep="-")
  
  ## SQ header; required by asBam
  seqlengths <- fasta.seqlengths(refFn)
  seqlengthsHeader <- paste0("@SQ\tSN:", names(seqlengths), "\tLN:", seqlengths)
  
  ## Add headers
  cat("@HD\tVN:1.5\tSO:unknown", file=tempSamFn, sep="\n", append=FALSE)
  cat(seqlengthsHeader, file=tempSamFn, sep="\n", append=TRUE)
  
  ## parse the fastq files
  f <- FastqStreamer(fastqFn, 1e6)
  on.exit(close(f), add=TRUE)
  if(!is.null(fastqI1Fn)){
    stopifnot(!is.null(fastqI2Fn))
    f1 <- FastqStreamer(fastqI1Fn, 1e6)
    f2 <- FastqStreamer(fastqI2Fn, 1e6)
    on.exit(close(f1), add=TRUE)
    on.exit(close(f2), add=TRUE)
  }
  while(length(fq <- yield(f))){
    samText <- paste(id(fq), "4", "*", "0", "0", "*", "*", "0", "0",
                     sread(fq), quality(quality(fq)), sep="\t")
    if(!is.null(fastqI1Fn)){
      fq1 <- yield(f1)
      fq2 <- yield(f2)
      samText <- paste0(samText, "\t", 
                        "BC:Z:", sread(fq2), "\t",
                        "QT:Z:", quality(quality(fq2)), "\t",
                        "CR:Z:", sread(fq1), "\t",
                        "CY:Z:", quality(quality(fq1)))
    }
    cat(samText, file=tempSamFn, append=TRUE, sep="\n")
  }
  tempBamFn <- asBam(file=tempSamFn, destination=tempBamFn, 
                     indexDestination=FALSE)
  bamFn <- sortBam(file=tempBamFn, destination=bamFn, byQname=TRUE, 
                   maxMemory=1024)
  on.exit(file.remove(tempSamFn), add=TRUE)
  on.exit(file.remove(tempBamFn), add=TRUE)
  invisible(bamFn)
}

### convert fastq files into a bam file, with read group tags
fastqs2bam <- function(fastqFns, fastq2Fns=NULL, readGroupNames=NULL,
                       bamFn, platform="illumina", mc.cores=ezThreads()){
  require(Biostrings)
  paired <- FALSE
  if(!is.null(fastq2Fns)){
    stopifnot(length(fastqFns) == length(fastq2Fns))
    paired <- TRUE
  }
  sampleBasenames <- sub("\\.(fastq|fq)(\\.gz){0,1}$", "",
                         basename(fastqFns))
  if(is.null(readGroupNames)){
    readGroupNames <- sampleBasenames
  }else{
    stopifnot(length(fastqFns) == length(readGroupNames))
  }
  
  emptyFastqs <- countReadsInFastq(fastqFns) == 0L
  ## tempty fastq files fail in the picard tool of FastqToSam and MergeSamFiles
  ## Convert and merge the non-empty fastqs first and alter the header of merged bam file
    
  cmd <- paste(preparePicard(), "FastqToSam",
               paste0("F1=", fastqFns)
               )
  if(isTRUE(paired)){
    cmd <- paste(cmd, paste0("F2=", fastq2Fns))
  }
  cmd <- paste(cmd, paste0("O=", sampleBasenames, ".bam"),
               paste0("SAMPLE_NAME=", readGroupNames),
               paste0("LIBRARY_NAME=", readGroupNames),
               paste0("READ_GROUP_NAME=", readGroupNames),
               paste0("RUN_DATE=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S+00:00")),
               paste0("PLATFORM=", platform),
               "SEQUENCING_CENTER=FGCZ")
  
  ezMclapply(cmd[!emptyFastqs], ezSystem, mc.preschedule=FALSE,
             mc.cores = mc.cores)
  
  ## picard tools merge
  cmd <- paste(preparePicard(), "MergeSamFiles",
               paste0("I=", paste0(sampleBasenames, ".bam")[!emptyFastqs],
                      collapse=" "),
               paste0("O=", bamFn),
               "USE_THREADING=true", "SORT_ORDER=queryname")
  ezSystem(cmd)
  file.remove(paste0(sampleBasenames, ".bam"))
  
  if(any(emptyFastqs)){
    tempHeaderFn <- tempfile(pattern="nonEmpty", tmpdir=getwd(),
                             fileext=".header")
    cmd <- paste("samtools view -H", bamFn, ">", tempHeaderFn)
    ezSystem(cmd)
    
    extraHeaders <- paste("@RG",
                          paste0("ID:", readGroupNames[emptyFastqs]),
                          paste0("SM:", readGroupNames[emptyFastqs]),
                          paste0("LB:", readGroupNames[emptyFastqs]),
                          paste0("PL:", platform),
                          "CN:FGCZ",
                          paste0("DT:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S+00:00")),
                          sep="\t")
    con <- file(tempHeaderFn, open="a")
    writeLines(extraHeaders, con=con)
    close(con)
    
    bamReheaderFn <- tempfile(pattern="reheader", tmpdir=getwd(),
                              fileext = ".bam")
    cmd <- paste("samtools reheader", tempHeaderFn, bamFn, ">", bamReheaderFn)
    ezSystem(cmd)
    
    file.remove(bamFn)
    file.rename(from=bamReheaderFn, to=bamFn)
  }
  
  invisible(bamFn)
}

### Convert  bam/sam files into fastqs
bam2fastq <- function(bamFn, OUTPUT_PER_RG=TRUE, OUTPUT_DIR=".",
                      paired=FALSE,
                      fastqFns=sub("(\\.bam|\\.sam)$", "_R1.fastq", bamFn),
                      fastq2Fns=sub("(\\.bam|\\.sam)$", "_R2.fastq", bamFn)){
  if(isTRUE(OUTPUT_PER_RG)){
    ## When OUTPUT_PER_RG is TRUE, we only process one uBam each time.
    stopifnot(length(bamFn) == 1L)
    ## I don't want to parse the bam header to get RG IDs
    ## Put them in a tempdir and move to OUTPUT_DIR later
    tempDIR <- paste("SamtoFastqTempDir", Sys.getpid(), sep="-")
    dir.create(tempDIR)
    on.exit(unlink(tempDIR, recursive=TRUE), add = TRUE)
    cmd <- paste(preparePicard(), "SamToFastq",
                 paste0("I=", bamFn),
                 paste0("OUTPUT_DIR=", tempDIR),
                 "OUTPUT_PER_RG=true RG_TAG=ID"
                 )
    ezSystem(cmd)
    fastqFns <- list.files(path=tempDIR, pattern="_1\\.fastq$")
    fromFns <- file.path(tempDIR, fastqFns)
    toFns <- file.path(OUTPUT_DIR, sub("_1\\.fastq$", "_R1.fastq", fastqFns))
    file.rename(from=fromFns, to=toFns)
    if(isTRUE(paired)){
      file.rename(from=sub("_1\\.fastq$", "_2.fastq", fromFns), 
                  to=sub("_R1\\.fastq$", "_R2.fastq", toFns)
                  )
    }
    return(invisible(toFns))
    ## This is not much slower than splitBambyRG and SamToFastq in parallel
  }else{
    cmd <- paste(preparePicard(), "SamToFastq",
                 paste0("I=", bamFn),
                 paste0("FASTQ=", fastqFns))
    if(isTRUE(paired))
      cmd <- paste(cmd, paste0("SECOND_END_FASTQ=", fastq2Fns))
    lapply(cmd, ezSystem)
    return(invisible(fastqFns))
  }
}

ezMethodBam2Fastq <- function(input=NA, output=NA, param=NA,
                              OUTPUT_PER_RG=TRUE){
  require(Biostrings)
  
  if(isTRUE(OUTPUT_PER_RG)){
    output = EzDataset(file=input$getFullPaths("CellDataset"),
                       dataRoot=param$dataRoot)
    output$setColumn("Read1", paste0(getwd(), "/", output$getNames(),
                                     "_R1.fastq"))
    if (param$paired){
      output$setColumn("Read2", paste0(getwd(), "/", output$getNames(), 
                                       "_R2.fastq"))
    } else {
      if ("Read2" %in% input$colNames)
        output$setColumn("Read2", NULL)
    }
    output$dataRoot = NULL
    
    allFastqFns <- bam2fastq(bamFn=input$getFullPaths("Read1"),
                             OUTPUT_PER_RG=TRUE, OUTPUT_DIR=getwd(),
                             paired=param$paired)
    fastqsToRemove <- allFastqFns[!(basename(allFastqFns) %in% 
                                      basename(output$getColumn("Read1")))]
    ## When we want only subset of the cells to keep
    file.remove(fastqsToRemove)
    file.remove(sub("_R1\\.fastq$", "_R2.fastq", fastqsToRemove))
  }else{
    if (!is(output, "EzDataset")){
      outputMeta <- input$meta
      outputMeta[['Read1']] <- paste0(getwd(), "/", rownames(outputMeta), 
                                   "_R1.fastq")
      outputMeta[['Read1 [File]']] = NULL
      if (param$paired){
        outputMeta[['Read2']] <- paste0(getwd(), "/", rownames(outputMeta), 
                                     "_R2.fastq")
        outputMeta[['Read2 [File]']] = NULL
      } else {
        outputMeta[['Read2']] <- NULL
      }
      output <- EzDataset(meta=outputMeta, dataRoot=NULL)
    }
    bam2fastq(bamFn=input$getFullPaths("Read1"),
              OUTPUT_PER_RG=FALSE,
              fastqFns=output$getColumn("Read1"),
              fastq2Fns=ifelse(isTRUE(param$paired), output$getColumn("Read1"),
                               list(NULL)),
              paired=param$paired)
    output$setColumn("Read Count",
                     countReadsInFastq(output$getColumn("Read1")))
  }
  return(output)
}

countReadsInFastq = function(fastqFiles){
  require(Biostrings)
  nReads <- sapply(fastqFiles, fastq.geometry)[1, ]
  return(nReads)
}

ezMethodSubsampleFastq <- function(input=NA, output=NA, param=NA, n=1e6){
  require(ShortRead)
  ## if output is not an EzDataset, set it!
  if (!is(output, "EzDataset")){
    output = input$copy()
    output$setColumn("Read1", paste0(getwd(), "/", input$getNames(), 
                                     "-subsample_R1.fastq.gz"))
    if (param$paired){
      output$setColumn("Read2", paste0(getwd(), "/", input$getNames(),
                                       "-subsample_R2.fastq.gz"))
    } else {
      if ("Read2" %in% input$colNames){
        output$setColumn("Read2", NULL)
      }
    }
    output$dataRoot = NULL
  }
  
  dataset = input$meta
  samples = rownames(dataset)
  mclapply(samples, function(sm, input, output, param){
    fl <- input$getFullPaths("Read1")[sm]
    f1 <- FastqSampler(fl, n=n, ordered=TRUE)
    set.seed(123L)
    p1 <- yield(f1)
    close(f1)
    writeFastq(p1, file=output$getColumn("Read1")[sm])
    if(param$paired){
      fl <- input$getFullPaths("Read2")[sm]
      f1 <- FastqSampler(fl, n=n, ordered=TRUE)
      set.seed(123L)
      p1 <- yield(f1)
      close(f1)
      writeFastq(p1, file=output$getColumn("Read2")[sm])
    }
  }, input=input, output=output, param=param, mc.cores=min(4L, param$cores))
  return(output)
}
