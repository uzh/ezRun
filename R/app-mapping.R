# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodBowtie2 <- function(input = NA, output = NA, param = NA) {
  param$fastpCompression = 9
  ref <- getBowtie2Reference(param)
  bamFile <- output$getColumn("BAM")
  sampleName <- sub(".bam", "", basename(bamFile))
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  defOpt <- paste("-p", param$cores)
  readGroupOpt <- paste0(
    "--rg-id ", sampleName, " --rg SM:", sampleName,
    " --rg LB:RGLB_", sampleName,
    " --rg PL:illumina", " --rg PU:RGPU_", sampleName
  )
  cmd <- paste(
    "bowtie2", param$cmdOptions, defOpt, readGroupOpt,
    "-x", ref, if (param$paired) "-1", trimmedInput$getColumn("Read1"),
    if (param$paired) paste("-2", trimmedInput$getColumn("Read2")),
    "2>", paste0(sampleName, "_bowtie2.log"), "|",
    "samtools", "view -S -b -", " > bowtie.bam"
  )
  ezSystem(cmd)
  file.remove(trimmedInput$getColumn("Read1"))
  if (param$paired) {
    file.remove(trimmedInput$getColumn("Read2"))
  }
  
  if (!is.null(param$markDuplicates) && param$markDuplicates) {
    ezSortIndexBam("bowtie.bam", "sorted.bam",
                   ram = param$ram, removeBam = TRUE,
                   cores = param$cores
    )
    dupBam(inBam = "sorted.bam", outBam = basename(bamFile), operation = "mark", ram = param$ram, dupDistance = param$dupDistance)
    file.remove("sorted.bam")
  } else {
    ezSortIndexBam("bowtie.bam", basename(bamFile),
                   ram = param$ram,
                   removeBam = TRUE, cores = param$cores
    )
  }
  
  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  
  if (param$generateBigWig) {
    bam2bw(file = basename(bamFile), paired = param$paired, method = "Bioconductor", cores = param$cores)
  }
  
  if (ezIsSpecified(param$secondRef) & param$secondRef!='') {
      require(Rsamtools)
      require(GenomicAlignments)
      
      bam_file <- basename(bamFile)
      bam <- readGAlignments(bam_file)
      
      chrFile <- file.path(dirname(dirname(file.path('/srv/GT/reference/',param$refBuild))), 'Sequence/WholeGenomeFasta/genome-chromsizes.txt')
      chrSizes <- ezRead.table(chrFile, header = FALSE)
      
      reads_per_chromosome <- table(seqnames(bam))
      chrLengths <- seqlengths(seqinfo(bam))[names(reads_per_chromosome)]
      chrLengths <- chrLengths[names(reads_per_chromosome)]
      readSize <- mean(qwidth(bam))
      stats <- data.frame(seqname = names(reads_per_chromosome), readsPerChr= as.vector(reads_per_chromosome), lengthPerChr = chrLengths, avgCov = as.vector(reads_per_chromosome*readSize/chrLengths))
      myChr <- levels(seqnames(bam))
      myChr <- intersect(levels(seqnames(bam)),names(reads_per_chromosome))
      
      extraChr <- setdiff(myChr, rownames(chrSizes))
      extraChrCov <- list()
      for (k in 1:length(extraChr)) {
          extraChrCov[[k]] <- coverage(bam)[[extraChr[k]]]
      }
      covExtraChr <- stats[extraChr,]$avgCov
      remove(bam)
     
      png(paste0('Coverage_', sampleName, '_secondRef.png'), width = 1200, height = 500*length(extraChr), res = 100)
      par(mfrow = c(length(extraChr),1))
      for (k in 1:length(extraChr)) {
          plot(extraChrCov[[k]], type = "l", col = "blue", lwd = 2, xlab = "Genomic position", ylab = "Coverage", main = paste("Coverage of", extraChr[k], 'in', sampleName))
          abline(h = covExtraChr[k], col = "red", lwd = 2)
          text(0.9 *chrLengths[[extraChr[k]]], covExtraChr[k], paste("Mean coverage:", round(covExtraChr[k], 2)), pos = 3, col = "black")
      }
      dev.off()
      
      write.table(
          '\n',
          file = paste0(sampleName, "_bowtie2.log"),
          append = TRUE,        
          row.names = FALSE,    
          col.names = TRUE,    
          sep = ",",           
          quote = FALSE
      )
      
      write.table(
          stats,
          file = paste0(sampleName, "_bowtie2.log"),
          append = TRUE,        
          row.names = FALSE,    
          col.names = TRUE,    
          sep = ",",           
          quote = FALSE
      )
  }
  
  return("Success")
}

##' @template getref-template
##' @templateVar methodName Bowtie2
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refBuildDir}{ a character specifying the directory of the reference build.}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getBowtie2Reference <- function(param) {
  ## the refBase specifies the common stem of the index files with suffix .bt2
  if (ezIsSpecified(param$ezRef["refIndex"])) {
    refBase <- param$ezRef["refIndex"]
    stopifnot(file.exists(dirname(refBase)))
    return(stopifnot)
  }
  ## default values
  refBase <- file.path(param$ezRef["refBuildDir"], "Sequence/BOWTIE2Index/genome")
  genomeFastaFiles <- param$ezRef["refFastaFile"]
  ## update defaults if second ref exists
  ## build a temporary index if a second reference exists
  if (ezIsSpecified(param$secondRef)) {
    stopifnot(file.exists(param$secondRef))
    genomeFastaFiles <- paste0(genomeFastaFiles, ",", param$secondRef)
    refBase <- file.path(getwd(), "BOWTIE2Index/customGenome")
  }
  
  ## check the ref
  lockFile <- file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))) {
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con = lockFile)
    wd <- getwd()
    setwd(dirname(refBase))
    
    cmd <- paste("bowtie2-build", "--seed 42 --threads", as.numeric(param$cores), "-f", genomeFastaFiles, basename(refBase))
    ezSystem(cmd)
    # ezWriteElapsed(job, "done")
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refBase)))
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles <- list.files(dirname(refBase), basename(refBase))
  if (length(refFiles) < 2) {
    ## too few files, we assume the index is incomplete
    stop(paste("index not available: ", refBase))
  }
  return(refBase)
}

##' @template app-template
##' @templateVar method ezMethodBowtie2(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getBowtie2Reference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppBowtie2 <-
  setRefClass("EzAppBowtie2",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBowtie2
                  name <<- "EzAppBowtie2"
                  appDefaults <<- rbind(
                    writeIgvSessionLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated"),
                    markDuplicates = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should duplicates be marked"),
                    generateBigWig = ezFrame(Type = "logical", DefaultValue = "FALSE", Description = "should a bigwig file be generated")
                  )
                }
              )
  )

ezMethodBowtie <- function(input = NA, output = NA, param = NA) {
  ref <- getBowtieReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  defOpt <- paste("--chunkmbs 256", "--sam", "-p", param$cores)
  cmd <- paste(
    "bowtie", param$cmdOptions, defOpt,
    ref, if (param$paired) "-1", trimmedInput$getColumn("Read1"),
    if (param$paired) paste("-2", trimmedInput$getColumn("Read2")),
    "2> bowtie.log", "|", "samtools", "view -S -b -", " > bowtie.bam"
  )
  ezSystem(cmd)
  ezSortIndexBam("bowtie.bam", basename(bamFile),
                 ram = param$ram, removeBam = TRUE,
                 cores = param$cores
  )
  
  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName Bowtie
##' @inheritParams getBowtie2Reference
getBowtieReference <- function(param) {
  refBase <- ifelse(param$ezRef["refIndex"] == "",
                    file.path(param$ezRef["refBuildDir"], "Sequence/BOWTIEIndex/genome"),
                    param$ezRef["refIndex"]
  )
  ## check the ref
  lockFile <- file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))) {
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con = lockFile)
    wd <- getwd()
    setwd(dirname(refBase))
    
    fastaFile <- param$ezRef["refFastaFile"]
    # job = ezJobStart("bowtie index")
    ezSystem(paste("ln -s", fastaFile, "."))
    buildOpts <- ""
    if (any(grepl("--threads", system("bowtie-build --help", intern = T)))) {
      buildOpts <- paste("--threads", param$cores)
    }
    cmd <- paste("bowtie-build", buildOpts, "-f", basename(fastaFile), basename(refBase))
    ezSystem(cmd)
    # ezWriteElapsed(job, "done")
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refBase)))
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles <- list.files(dirname(refBase), basename(refBase))
  if (length(refFiles) < 3) {
    ## we assume the index is built and complete
    stop(paste("index not available: ", refBase))
  }
  return(refBase)
}

##' @template app-template
##' @templateVar method ezMethodBowtie(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getBowtieReference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppBowtie <-
  setRefClass("EzAppBowtie",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBowtie
                  name <<- "EzAppBowtie"
                  appDefaults <<- rbind(writeIgvSessionLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated"))
                }
              )
  )

ezMethodSTAR <- function(input = NA, output = NA, param = NA) {
  refDir <- getSTARReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  
  if(ezIsSpecified(param$barcodePattern) && param$barcodePattern!=''){ #Extract UMI
      require(Herper)
      local_CondaEnv("gi_umi_tools", pathToMiniConda = "/usr/local/ngseq/miniforge3")
      ##Extract UMI from R2
      markedFile_R1 <- sub('R1', 'markedUMI_R1', trimmedInput$getColumn("Read1"))
      markedFile_R2 <- sub('R2', 'markedUMI_R2', trimmedInput$getColumn("Read2"))
      cmd <- paste0('umi_tools extract --stdin=',trimmedInput$getColumn("Read2"), ' --read2-in=',trimmedInput$getColumn("Read1"), 
                    ' --stdout=',markedFile_R2,' --read2-out=',markedFile_R1,' --bc-pattern=', param$barcodePattern)
      ezSystem(cmd)
      ezSystem(paste('mv', markedFile_R1, trimmedInput$getColumn("Read1")))
      ezSystem(paste('mv', markedFile_R2, trimmedInput$getColumn("Read2")))
  }
  
  if (!str_detect(param$cmdOptions, "outSAMattributes")) {
    param$cmdOptions <- str_c(param$cmdOptions, "--outSAMattributes All",
                              sep = " "
    )
  }
  
  ## we use the --genomeFastaFiles option only if we have spliced sequences; 
  ## if we have .fa and .gtf we build a separate index
  if (ezIsSpecified(param$secondRef)) {
    stopifnot(file.exists(param$secondRef))
    secondGtf <- sub(".fa", ".gtf", param$secondRef)
    if (!file.exists(secondGtf)){
      genomeFastaFileOption <- str_c("--genomeFastaFiles", param$secondRef, sep = " ")
    }
  } else {
    genomeFastaFileOption <- ""
  }
  
  
  cmd <- str_c(
    "STAR", "--genomeDir", refDir,
    "--readFilesIn", trimmedInput$getColumn("Read1"),
    if (param$paired) trimmedInput$getColumn("Read2"),
    "--twopassMode", if_else(param$twopassMode, "Basic", "None"),
    "--runThreadN", param$cores, param$cmdOptions,
    genomeFastaFileOption,
    "--outStd BAM_Unsorted --outSAMtype BAM Unsorted",
    "--outSAMattrRGline", str_c("ID:", trimmedInput$getNames()), str_c("SM:", trimmedInput$getNames()),
    if_else(str_detect(trimmedInput$getColumn("Read1"), "\\.gz$"), "--readFilesCommand zcat", ""),
    ">  Aligned.out.bam",
    sep = " "
  )
  ezSystem(cmd)
  
  nSortThreads <- min(param$cores, 8)
  ## if the index is loaded in shared memory we have to use only 10% of the scheduled RAM
  if (str_detect(param$cmdOptions, "--genomeLoad Load")) {
    sortRam <- param$ram / 10
  } else {
    sortRam <- param$ram
  }
  
  file.rename(
    from = "Log.final.out",
    to = basename(output$getColumn("STARLog"))
  )
  
  if (!is.null(param$markDuplicates) && param$markDuplicates) {
    ezSortIndexBam("Aligned.out.bam", "sorted.bam",
                   ram = sortRam, removeBam = TRUE,
                   cores = nSortThreads
    )
    dupBam(inBam = "sorted.bam", outBam = basename(bamFile), operation = "mark", ram = param$ram)
    file.remove("sorted.bam")
  } else {
    ezSortIndexBam("Aligned.out.bam", basename(bamFile),
                   ram = sortRam,
                   removeBam = TRUE, cores = nSortThreads
    )
  }
  
  if (param$getJunctions) {
    file.rename(
      from = "SJ.out.tab",
      to = basename(output$getColumn("Junctions"))
    )
    file.rename(
      from = "Chimeric.out.junction",
      to = basename(output$getColumn("Chimerics"))
    )
  }
  
  if(ezIsSpecified(param$barcodePattern) && param$barcodePattern!=''){ #Deduplicated based on UMI
      deDupBamFile <- sub('.bam', '_dedup.bam', basename(bamFile))
      cmd <- paste0('umi_tools dedup --stdin=',basename(bamFile),' --stdout=', deDupBamFile,' --log=',basename(bamFile),'.log', ' --output-stats=',basename(bamFile),'.stats')
      ezSystem(cmd)
      ezSystem(paste('mv', deDupBamFile, basename(bamFile)))
      ezSystem(paste('samtools index', basename(bamFile)))
      tryCatch({local_CondaEnv("gi_py3.11.5", pathToMiniConda = "/usr/local/ngseq/miniforge3")}, warning = function(warning_condition) {cat('warning')}, 
               error = function(error_condition) {cat('error')})
  }
  
  
  ## check the strandedness
  tryCatch({
    ezSystem(str_c(
      "infer_experiment.py", "-r", getReferenceFeaturesBed(param),
      "-i", basename(bamFile), "-s 1000000",
      sep = " "
    ), stopOnFailure = FALSE)
  }, error = function(e) {
    warning(paste0("Could not get the reference features for the ", 
                   "strandedness check. Skipping."))
  })
  
  ## write an igv link
  if (param$writeIgvLink && "IGV" %in% output$colNames) {
    writeIgvHtml(param, output)
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName STAR
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refFeatureFile}{ a character specifying the file path to the annotation feature file (.gtf).}
##'   \item{ram}{ an integer specifying how many gigabytes of RAM to use.}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getSTARReference <- function(param) {
  if (ezIsSpecified(param$ezRef["refIndex"])) {
    refDir <- param$ezRef["refIndex"]
    stopifnot(file.exists(file.path(refDir, "SAindex")))
    return(refDir)
  }
  if (!ezIsSpecified(param$ezRef["refFeatureFile"])) {
    stop("refFeatureFile not defined")
  }
  
  ## default values  
  gtfFile <- param$ezRef["refFeatureFile"]
  genomeFastaFiles <- param$ezRef["refFastaFile"]
  refDir <- sub(".gtf$", "_STARIndex", param$ezRef["refFeatureFile"])
  ## update if second .fa and .gtf exists
  if (ezIsSpecified(param$secondRef)) {
    stopifnot(file.exists(param$secondRef))
    secondGtf <- sub(".fa", ".gtf", param$secondRef)
    if (file.exists(secondGtf)) {
      gtfFile <- tempfile(
        pattern = "genes", tmpdir = getwd(),
        fileext = ".gtf"
      )
      gtf1 <- ezRead.table(param$ezRef["refFeatureFile"], quote="", sep="\t", comment.char = "#", row.names = NULL, header = FALSE)
      gtf2 <- ezRead.table(secondGtf, quote="", sep="\t", comment.char = "#", row.names = NULL, header = FALSE)
      ## if the same chromosome name shows up in both GTF files, we can't concatenate them
      stopifnot(length(intersect(gtf1$V1, gtf2$V1)) == 0) ## first column holds the seqid
      ezSystem(paste("cp", param$ezRef["refFeatureFile"], gtfFile))
      ezSystem(paste("cat", secondGtf, ">>", gtfFile))
      genomeFastaFiles <- paste(param$ezRef["refFastaFile"], param$secondRef)
      refDir <- file.path(getwd(), "Custom_STARIndex")
    }
  }
  
  ## random sleep to avoid parallel ref building
  Sys.sleep(runif(1, max = 20))
  
  lockFile <- paste0(refDir, ".lock")
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (i >= INDEX_BUILD_TIMEOUT) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  if (file.exists(file.path(refDir, "SAindex"))) {
    ## we assume the index is built and complete
    return(refDir)
  }
  
  ## no lock file and no refFiles, so we build the reference
  ezWrite(Sys.info(), con = lockFile)
  dir.create(refDir)
  
  fai <- ezRead.table(paste0(param$ezRef["refFastaFile"], ".fai"), header = FALSE)
  colnames(fai) <- c("LENGTH", "OFFSET", "LINEBASES", "LINEWDITH")
  if (nrow(fai) > 50) {
    binOpt <- "--genomeChrBinNbits 16"
  } else {
    binOpt <- ""
  }
  
  genomeLength <- sum(fai$LENGTH)
  readLength <- 150 ## assumption
  indexNBasesOpt <- paste("--genomeSAindexNbases", min(13, floor(log2(genomeLength) / 2 - 1)))
  if (binOpt == "") {
    genomeChrBinNbits <- paste("--genomeChrBinNbits", floor(min(
      18,
      log2(max(genomeLength / nrow(fai), readLength))
    )))
  } else {
    genomeChrBinNbits <- ""
  }
  
  job <- ezJobStart("STAR genome build")
  cmd <- paste(
    "STAR", "--runMode genomeGenerate --genomeDir", refDir, binOpt, indexNBasesOpt, genomeChrBinNbits,
    "--limitGenomeGenerateRAM", format(param$ram * 1e9, scientific = FALSE),
    "--genomeFastaFiles", genomeFastaFiles,
    "--sjdbGTFfile", gtfFile, "--sjdbOverhang 150", "--runThreadN", param$cores, "--genomeSAsparseD 2"
  )
  ezSystem(cmd)
  file.remove(lockFile)
  ezWriteElapsed(job, "done")
  file.remove("Log.out")
  return(refDir)
}

##' @template app-template
##' @templateVar method ezMethodSTAR(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getSTARReference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppSTAR <-
  setRefClass("EzAppSTAR",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSTAR
                  name <<- "EzAppSTAR"
                  appDefaults <<- rbind(
                    getJunctions = ezFrame(Type = "logical", DefaultValue = "FALSE", Description = "should junctions be returned"),
                    writeIgvLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated"),
                    markDuplicates = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should duplicates be marked with picard"),
                    twopassMode = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "1-pass mapping or basic 2-pass mapping")
                  )
                }
              )
  )

ezMethodBWA <- function(input = NA, output = NA, param = NA) {
  refIdx <- getBWAReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  sampleName <- sub(".bam", "", basename(bamFile))
  readGroupOpt <- paste0("\'@RG\\tID:",sampleName,"\\tSM:",sampleName,"\\tLB:",sampleName,"\\tPL:ILLUMINA\'")
  
  if (param$algorithm == "aln") {
    cmd <- paste(
      "bwa", param$algorithm, param$cmdOptions, "-R", readGroupOpt, "-t", param$cores,
      refIdx, trimmedInput$getColumn("Read1"), ">", "read1.sai", "2> bwa.log"
    )
    ezSystem(cmd)
    if (param$paired) {
      cmd <- paste(
        "bwa", param$algorithm, param$cmdOptions, "-R", readGroupOpt, "-t", param$cores,
        refIdx, trimmedInput$getColumn("Read2"), ">", "read2.sai", "2> bwa.log"
      )
      ezSystem(cmd)
      cmd <- paste(
        "bwa", "sampe", refIdx, "read1.sai", "read2.sai",
        trimmedInput$getColumn("Read1"), trimmedInput$getColumn("Read2"),
        "2> bwa.log", "|",
        "samtools", "view -S -b -", " > aligned.bam"
      )
      ezSystem(cmd)
    } else {
      cmd <- paste(
        "bwa", "samse", refIdx, "read1.sai",
        trimmedInput$getColumn("Read1"), "2> bwa.log", "|",
        "samtools", "view -S -b -", " > aligned.bam"
      )
      ezSystem(cmd)
    }
  } else {
    if (param$algorithm == "bwasw" && param$paired) {
      stop("paired is not supported for algorithm bwasw")
    }
    cmd <- paste(
      "bwa", param$algorithm, param$cmdOptions, "-R", readGroupOpt, "-t", param$cores,
      refIdx, trimmedInput$getColumn("Read1"),
      if (param$paired) trimmedInput$getColumn("Read2"),
      "2> bwa.log", "|", "samtools", "view -S -b -", " > aligned.bam"
    )
    message(cmd)
    system(cmd)
  }
  file.remove(trimmedInput$getColumn("Read1"))
  if (param$paired) {
    file.remove(trimmedInput$getColumn("Read2"))
  }
  
  if (!is.null(param$markDuplicates) && param$markDuplicates) {
    ezSortIndexBam("aligned.bam", "sorted.bam",
                   ram = param$ram, removeBam = TRUE,
                   cores = param$cores
    )
    dupBam(inBam = "sorted.bam", outBam = basename(bamFile), operation = "mark", ram = round(0.8 * param$ram))
    file.remove("sorted.bam")
  } else {
    ezSortIndexBam("aligned.bam", basename(bamFile),
                   ram = param$ram,
                   removeBam = TRUE, cores = param$cores
    )
  }
  
  
  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  return("Success")
}

ezMethodBWATrimmomatic <- function(input = NA, output = NA, param = NA) { # Perform BWA using fastp for read pre-processing
  
  refIdx <- getBWAReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodTrim(input = input, param = param)
  if (param$algorithm == "aln") {
    cmd <- paste(
      "bwa", param$algorithm, param$cmdOptions, "-t", param$cores,
      refIdx, trimmedInput$getColumn("Read1"), ">", "read1.sai", "2> bwa.log"
    )
    ezSystem(cmd)
    if (param$paired) {
      cmd <- paste(
        "bwa", param$algorithm, param$cmdOptions, "-t", param$cores,
        refIdx, trimmedInput$getColumn("Read2"), ">", "read2.sai", "2> bwa.log"
      )
      ezSystem(cmd)
      cmd <- paste(
        "bwa", "sampe", refIdx, "read1.sai", "read2.sai",
        trimmedInput$getColumn("Read1"), trimmedInput$getColumn("Read2"),
        "2> bwa.log", "|",
        "samtools", "view -S -b -", " > aligned.bam"
      )
      ezSystem(cmd)
    } else {
      cmd <- paste(
        "bwa", "samse", refIdx, "read1.sai",
        trimmedInput$getColumn("Read1"), "2> bwa.log", "|",
        "samtools", "view -S -b -", " > aligned.bam"
      )
      ezSystem(cmd)
    }
  } else {
    if (param$algorithm == "bwasw" && param$paired) {
      stop("paired is not supported for algorithm bwasw")
    }
    cmd <- paste(
      "bwa", param$algorithm, param$cmdOptions, "-t", param$cores,
      refIdx, trimmedInput$getColumn("Read1"),
      if (param$paired) trimmedInput$getColumn("Read2"),
      "2> bwa.log", "|", "samtools", "view -S -b -", " > aligned.bam"
    )
    ezSystem(cmd)
  }
  file.remove(trimmedInput$getColumn("Read1"))
  if (param$paired) {
    file.remove(trimmedInput$getColumn("Read2"))
  }
  
  ezSortIndexBam("aligned.bam", basename(bamFile),
                 ram = param$ram, removeBam = TRUE,
                 cores = param$cores
  )
  
  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName BWA
##' @inheritParams getBowtie2Reference
getBWAReference <- function(param) {
  refPath <- ifelse(param$ezRef["refIndex"] == "",
                    file.path(param$ezRef["refBuildDir"], "Sequence/BWAIndex/genome.fa"),
                    param$ezRef["refIndex"]
  )
  ## check the ref
  lockFile <- file.path(dirname(refPath), "lock")
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  if (!file.exists(dirname(refPath))) {
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refPath))
    ezWrite(Sys.info(), con = lockFile)
    wd <- getwd()
    setwd(dirname(refPath))
    
    fastaFile <- param$ezRef["refFastaFile"]
    file.symlink(from = fastaFile, to = ".")
    cmd <- paste("bwa", "index", "-a", "bwtsw", basename(fastaFile))
    ezSystem(cmd)
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refPath)))
  if (!file.exists(paste0(refPath, ".sa"))) {
    stop(paste("sa index not found for:", refPath))
  }
  return(refPath)
}

##' @template app-template
##' @templateVar method ezMethodBWA(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getBWAReference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppBWA <-
  setRefClass("EzAppBWA",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBWA
                  name <<- "EzAppBWA"
                  appDefaults <<- rbind(
                    algorithm = ezFrame(Type = "character", DefaultValue = "mem", Description = "bwa's alignment algorithm. One of aln, bwasw, mem."),
                    writeIgvSessionLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated"),
                    markDuplicates = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should duplicates be marked")
                  )
                }
              )
  )

EzAppBWATrimmomatic <-
  setRefClass("EzAppBWATrimmomatic",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBWATrimmomatic
                  name <<- "EzAppBWATrimmomatic"
                  appDefaults <<- rbind(
                    algorithm = ezFrame(Type = "character", DefaultValue = "mem", Description = "bwa's alignment algorithm. One of aln, bwasw, mem."),
                    writeIgvSessionLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated")
                  )
                }
              )
  )



ezMethodMinimap2 <- function(input = NA, output = NA, param = NA) {
  refIdx <- getMinimapReference(param)
  
  bedFile <- "genes.bed"
  cmd <- paste("paftools.js gff2bed", param$ezRef["refFeatureFile"], ">", bedFile)
  
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  sampleName <- sub(".bam", "", basename(bamFile))
  cmd <- paste("minimap2", "-a", "-t", param$cores, 
               "-2", ## use separate threads for input and output
               "-K 500M", ## that's the default batch size of bases; don't know how to adapt that to the available param$ram
               "--junc-bed", bedFile,
               param$cmdOptions, refIdx, trimmedInput$getColumn("Read1"), "> output.sam")
  ezSystem(cmd)
  ezSortIndexBam("output.sam", basename(bamFile),
                 ram = param$ram, removeBam = TRUE,
                 cores = param$cores
  )
  file.remove(trimmedInput$getColumn("Read1"))
  
  ## write an igv link
  if (param$writeIgvLink) {
    if ("IGV" %in% output$colNames) {
      writeIgvHtml(param, output)
    }
  }
  return("Success")
}


getMinimapReference <- function(param) {
  ## this does not build a reference but simply uses the genome fasta
  ## see the getBWAReference as an example on how to build an actual reference
  return(param$ezRef["refFastaFile"])
}


##' @template app-template
##' @templateVar method ezMethodMinimap2(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
##' @seealso \code{\link{getMinimap2Reference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppMinimap2 <-
  setRefClass("EzAppMinimap2",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMinimap2
                  name <<- "EzAppMinimap2"
                  appDefaults <<- rbind(
                    writeIgvSessionLink = ezFrame(Type = "logical", DefaultValue = "TRUE", Description = "should an IGV link be generated")
                  )
                }
              )
  )




ezMethodBismark <- function(input = NA, output = NA, param = NA) {
  param$fastpCompression = 9
  ref <- getBismarkReference(param)
  bamFile <- output$getColumn("BAM")
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)
  defOpt <- paste("-p", max(2, param$cores / 2))
  if (param$paired) {
    cmd <- paste(
      "bismark", param$cmdOptions,
      "--path_to_bowtie", paste0("$Bowtie2", "/bin"), defOpt, ref,
      "-1", trimmedInput$getColumn("Read1"),
      if (param$paired) paste("-2", trimmedInput$getColumn("Read2")),
      "2> bismark.log"
    )
  } else {
    cmd <- paste(
      "bismark", param$cmdOptions,
      "--path_to_bowtie", paste0("$Bowtie2", "/bin"), defOpt, ref,
      trimmedInput$getColumn("Read1"),
      "2> bismark.log"
    )
  }
  ezSystem(cmd)
  bamFileNameBismark <- list.files(".", pattern = "bam$")
  reportFileNameBismark <- list.files(".", pattern = "report.txt$")
  reportFile <- paste0(names(bamFile), ".report.txt")
  ezSystem(paste("mv ", reportFileNameBismark, reportFile))
  if (param$deduplicate) {
    cmd <- paste("deduplicate_bismark", ifelse(param$paired, "-p", "-s"), bamFileNameBismark)
    ezSystem(cmd)
    bamFileNameBismark <- list.files(".", pattern = "deduplicated.bam$")
    deduplicationReportFile <- list.files(".", pattern = "deduplication_report.txt$")
    ezSystem(paste("cat ", deduplicationReportFile, ">>", reportFile))
  }
  
  cmd <- paste("bismark_methylation_extractor", ifelse(param$paired, "-p", "-s"), "--comprehensive --gzip", bamFileNameBismark, "--parallel", max(2, param$cores / 2))
  ezSystem(cmd)
  cmd <- paste("samtools", "view -S -b ", bamFileNameBismark, " > bismark.bam")
  ezSystem(cmd)
  ezSortIndexBam("bismark.bam", basename(bamFile),
                 ram = param$ram, removeBam = TRUE,
                 cores = param$cores
  )
  
  if (param$generateBigWig) {
    destination = sub("\\.bam$", "_Cov.bw", basename(bamFile), ignore.case = TRUE)
    bam2bw(file = basename(bamFile), destination = destination, paired = param$paired, method = "Bioconductor", cores = param$cores)
  }
  
  mBiasImages <- list.files(".", pattern = "png$")
  if ((length(mBiasImages)) > 0) {
    for (i in 1:length(mBiasImages)) {
      ezSystem(paste("mv ", mBiasImages[i], paste0(names(bamFile), ".M-bias_R", i, ".png")))
    }
  } else {
    ezSystem(paste("touch", paste0(names(bamFile), ".M-bias_R1.png")))
    ezSystem(paste("touch", paste0(names(bamFile), ".M-bias_R2.png")))
  }
  
  CpGFile <- list.files(".", pattern = "^CpG.*txt.gz$")
  ezSystem(paste('gunzip',CpGFile))
  CpGFile <- sub('.gz', '', CpGFile)
  
  cmd <- paste("bismark2bedGraph --scaffolds", CpGFile, "-o", names(bamFile))
  ezSystem(cmd)
  ezSystem(paste("mv ", CpGFile, paste0(names(bamFile), ".CpG_context.txt")))
  #ezSystem(paste('pigz --best', paste0(names(bamFile), ".CpG_context.txt")))
  
  splittingReportFile <- list.files(".", pattern = "splitting_report.txt$")
  ezSystem(paste("cat ", splittingReportFile, ">>", reportFile))
  ## write an igv link
  # if (param$writeIgvSessionLink){
  #  writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
  #                  bamUrls = paste(PROJECT_BASE_URL, bamFile, sep="/") )
  #  writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", bamFile),
  #               sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
  # }
 if(grepl('Lambda',ref)){
  library(ggplot2)
  dataFile <- paste0(names(bamFile),'.gz.bismark.cov.gz')
  system(paste('gunzip', dataFile))
  dataFile <- sub('.gz$', '', dataFile)
  minCov <- 20
  data <- ezRead.table(dataFile, header = FALSE, row.names = NULL)
  colnames(data) <- c('Ref', 'Pos1', 'Pos2', 'Meth', 'MethCount', 'UnmethCount')
  data[['Cov']] = data$MethCount + data$UnmethCount
      data <- data[data$Cov >= minCov,]
      p <- ggplot(data, aes(x=Ref, y=Meth)) +  geom_boxplot(fill='#A4A4A4', color="black") + theme_classic()
      p <- p + labs(title=sub('.gz.*', '', dataFile))
      ggsave(paste0(sub('.gz.*', '_Meth.png', dataFile)), p, width = 6, height = 6)
  system(paste('gzip', dataFile))
 }
  return("Success")
}


##' @template getref-template
##' @templateVar methodName Bismark
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refBuildDir}{ a character specifying the directory of the reference build.}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getBismarkReference <- function(param) {
  refBase <- ifelse(param$ezRef["refIndex"] == "",
                    file.path(param$ezRef["refBuildDir"], "Sequence/WholeGenomeFasta/Bisulfite_Genome/CT_conversion"),
                    param$ezRef["refIndex"]
  )
  ## check the ref
  lockFile <- file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))) {
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con = lockFile)
    wd <- getwd()
    cmd <- paste(
      "bismark_genome_preparation", dirname(param$ezRef["refFastaFile"]),
      "2> bismarkGenomePrep.log"
    )
    ezSystem(cmd)
    
    # ezWriteElapsed(job, "done")
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refBase)))
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles <- list.files(dirname(refBase), basename(refBase))
  if (length(refFiles) < 1) {
    ## we assume the index is built and complete
    stop(paste("index not available: ", refBase))
  }
  return(dirname(param$ezRef@refFastaFile))
}



##' @template app-template
##' @templateVar method ezMethodBismark(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppBismark <-
  setRefClass("EzAppBismark",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBismark
                  name <<- "EzAppBismark"
                  appDefaults <<- rbind(
                    generateBigWig = ezFrame(Type = "logical", DefaultValue = "FALSE", Description = "should a bigwig coverage file be generated")
                  )
                }
              )
  )
