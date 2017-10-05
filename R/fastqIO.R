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
  on.exit(close(f))
  if(!is.null(fastqI1Fn)){
    stopifnot(!is.null(fastqI2Fn))
    f1 <- FastqStreamer(fastqI1Fn, 1e6)
    f2 <- FastqStreamer(fastqI2Fn, 1e6)
    on.exit(close(f1))
    on.exit(close(f2))
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
  on.exit(file.remove(tempSamFn))
  on.exit(file.remove(tempBamFn))
  invisible(bamFn)
}

### convert fastq files into a bam file, with read group tags
fastqs2bam <- function(fastqFns, bamFn){
  if(!isTRUE(isValidEnvironments("picard"))){
    setEnvironments("picard")
  }
  if(!isTRUE(isValidEnvironments("samtools"))){
    setEnvironments("samtools")
  }
  sampleBasenames <- sub("\\.fastq(\\.gz){0,1}$", "", basename(fastqFns))
  cmd <- paste("java -jar", Sys.getenv("Picard_jar"), "FastqToSam",
               paste0("F1=", fastqFns), 
               paste0("O=", sampleBasenames, ".bam"),
               paste0("SM=", sampleBasenames),
               paste0("READ_GROUP_NAME=", LETTERS[1:4])
               #paste0("SORT_ORDER=", "unsorted")
               )
  lapply(cmd, ezSystem)
  
  ## Rsamtolls merge only keeps the header from first bam
  #cmd <- paste("samtools merge", bamFn, 
  #             paste0(sampleBasenames, ".bam", collapse=" "))
  #ezSystem(cmd)
  
  ## picard tools merge
  cmd <- paste("java -jar", Sys.getenv("Picard_jar"), "MergeSamFiles",
               paste0("I=", paste0(sampleBasenames, ".bam"), collapse=" "),
               paste0("O=", bamFn),
               "SORT_ORDER=queryname")
  ezSystem(cmd)
}
