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
                       bamFn){
  if(!isTRUE(isValidEnvironments("picard"))){
    setEnvironments("picard")
  }
  paired <- FALSE
  if(!is.null(fastq2Fns)){
    stopifnot(length(fastqFns) == length(fastq2Fns))
    paired <- TRUE
  }
  sampleBasenames <- sub("\\.(fastq|fq)(\\.gz){0,1}$", "",
                         basename(fastqFns))
  if(is.null(readGroupNames))
    readGroupNames <- sampleBasenames
  cmd <- paste("java -jar", Sys.getenv("Picard_jar"), "FastqToSam",
               paste0("F1=", fastqFns),
               ifelse(isTRUE(paired), paste0("F2=", fastq2Fns), ""),
               paste0("O=", sampleBasenames, ".bam"),
               paste0("SM=", basename(fastqFns)),
               paste0("READ_GROUP_NAME=", readGroupNames)
               )
  lapply(cmd, ezSystem)
  ## picard tools merge
  cmd <- paste("java -jar", Sys.getenv("Picard_jar"), "MergeSamFiles",
               paste0("I=", paste0(sampleBasenames, ".bam"), collapse=" "),
               paste0("O=", bamFn),
               "SORT_ORDER=queryname")
  ezSystem(cmd)
  file.remove(paste0(sampleBasenames, ".bam"))
  
  invisible(bamFn)
}

### Convert  bam/sam files into fastqs
bam2fastq <- function(bamFns,
                      fastqFns=sub("(\\.bam|\\.sam)$", "_R1.fastq", bamFns),
                      fastq2Fns=sub("(\\.bam|\\.sam)$", "_R2.fastq", bamFns),
                      paired=FALSE){
  
  if(!isTRUE(isValidEnvironments("picard"))){
    setEnvironments("picard")
  }
  cmd <- paste("java -jar", Sys.getenv("Picard_jar"), "SamToFastq",
               paste0("I=", bamFns),
               paste0("FASTQ=", fastqFns),
               ifelse(isTRUE(paired), PASTE0("SECOND_END_FASTQ=", fastq2Fns),
                      ""))
  lapply(cmd, ezSystem)
  invisible(fastqFns)
}

ezMethodBam2Fastq <- function(input=NA, output=NA, param=NA){
  ## if output is not an EzDataset, set it!
  if (!is(output, "EzDataset")){
    output = input$copy()
    output$setColumn("Read1", paste0(getwd(), "/", input$getNames(), 
                                     "-R1.fastq"))
    if (param$paired){
      output$setColumn("Read2", paste0(getwd(), "/", input$getNames(), 
                                       "-R2.fastq"))
    } else {
      if ("Read2" %in% input$colNames){
        output$setColumn("Read2", NULL)
      }
    }
    output$dataRoot = NULL
  }
  
  bam2fastq(bamFns=input$getFullPaths("Read1"),
            fastqFns=output$getColumn("Read1"),
            fastq2Fns=ifelse(isTRUE(param$paired, output$getFullPaths("Read2"))),
            paired=param$paired)
  return(output)
}
