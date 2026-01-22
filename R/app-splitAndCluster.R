ezMethodSplitAndCluster = function(input = NA, output = NA, param = NA) {
  cdHitOpt = "-c 0.98 -d 0" ## TODO put this options in the sushi interface

  fastqFile1 = input$getFullPaths("Read1")
  fastqFile2 = input$getFullPaths("Read2")
  localFastq1 = sub(".gz", "", basename(fastqFile1))
  localFastq2 = sub(".gz", "", basename(fastqFile2))
  fastqJoined = "joined.fastq"
  ## <<- makes global vars
  FASTQ_JOIN <<-
    "/usr/local/ngseq/src/ea-utils.1.1.2-686/bin/fastq-join"
  CD_HIT <<- " /usr/local/ngseq/bin/cd-hit"
  ezSystem(paste("gunzip -c", fastqFile1, ">", localFastq1))
  ezSystem(paste("gunzip -c", fastqFile2, ">", localFastq2))
  fastqJoinCmd = paste(
    FASTQ_JOIN,
    localFastq1,
    localFastq2,
    "-o unjoined_R1.fastq",
    "-o unjoined_R2.fastq",
    "-o",
    fastqJoined,
    sep = " "
  )
  ezSystem(fastqJoinCmd)
  forwardAdapter = readDNAStringSet(param$forwardPrimerFile)
  reverseAdapter = readDNAStringSet(param$reversePrimerFile)
  require(ShortRead)
  nYield = 1e5
  max.mismatch = 0

  ## TODO make sure that directory does not yet exist
  fqs = FastqStreamer(fastqJoined, nYield)
  setwdNew(basename(output$getColumn("Clustered")))
  outputFiles = character()
  while (length(x <- yield(fqs))) {
    outputFiles = union(
      outputFiles,
      splitByAdapters(x, forwardAdapter, reverseAdapter, max.mismatch = 0)
    )
  }
  close(fqs)
  for (outputFile in outputFiles) {
    if (
      file.info(outputFile)$size > 0 & basename(outputFile) != "none.none.fasta"
    ) {
      cdHitOutForward = paste0(
        "cdHit_",
        sub(".fasta", "", basename(outputFile))
      )
      cdHitCmd = paste(
        CD_HIT,
        cdHitOpt,
        "-i",
        outputFile,
        "-o",
        cdHitOutForward,
        sep = " "
      )
      ezSystem(cdHitCmd)
    }
  }
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSplitAndCluster(input=NA, output=NA, param=NA)
##' @description Cluster stuff
EzAppSplitAndCluster <-
  setRefClass(
    "EzAppSplitAndCluster",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSplitAndCluster
        name <<- "EzAppSplitAndCluster"
        appDefaults <<- rbind(
          maxMismatch = ezFrame(
            Type = "integer",
            DefaultValue = "0",
            Description = "number of mismatches allowed in primer search"
          )
        )
      }
    )
  )


splitByAdapters <- function(
  reads,
  forwardAdapter,
  reverseAdapter,
  max.mismatch
) {
  foundFwdAdapterNames = rep("none", length(reads))
  adapterStartPos = rep(1, length(reads))
  foundRevAdapterNames = rep("none", length(reads))
  adapterEndPos = width(reads)

  readSeq = sread(reads)
  names(readSeq) = sub(" .*", "", id(reads))
  for (i in 1:(length(names(forwardAdapter)))) {
    faName = names(forwardAdapter)[i]
    vp1 = vmatchPattern(
      forwardAdapter[[faName]],
      readSeq,
      max.mismatch = max.mismatch
    )
    start = sapply(startIndex(vp1), function(x) {
      ifelse(is.null(x), -1, x[1])
    })
    hasFwdAdapt = start >= 1 & start <= 10
    if (any(hasFwdAdapt)) {
      foundFwdAdapterNames[hasFwdAdapt] = faName
      adapterStartPos[hasFwdAdapt] = start[hasFwdAdapt]
    }
  }

  for (i in 1:(length(names(reverseAdapter)))) {
    faNameRev = names(reverseAdapter)[i]
    vp1 = vmatchPattern(
      reverseAdapter[[faNameRev]],
      readSeq,
      max.mismatch = max.mismatch
    )
    end = sapply(endIndex(vp1), function(x) {
      ifelse(is.null(x), -1, x[1])
    })
    hasRevAdapt = end >= width(readSeq) - 10
    if (any(hasRevAdapt)) {
      foundRevAdapterNames[adapterEndPos] = faNameRev
      adapterEndPos[hasRevAdapt] = end[hasRevAdapt]
    }
  }
  idsByAdapter = split(
    1:length(reads),
    list(foundFwdAdapterNames, foundRevAdapterNames)
  )
  idsByAdapter = idsByAdapter[sapply(idsByAdapter, length) > 0]
  outputFiles = character()
  for (adapterCombName in names(idsByAdapter)) {
    ## adapterCombName will be the <adapterName1>.<adapterName2>
    indices = idsByAdapter[[adapterCombName]]
    trimmedReads = subseq(
      readSeq[indices],
      start = adapterStartPos[indices],
      end = adapterEndPos[indices]
    )
    nameOfFile = paste0(adapterCombName, ".fasta")
    writeFasta(trimmedReads, file = nameOfFile, mode = "a")
    outputFiles[adapterCombName] = nameOfFile
  }
  return(outputFiles)
}
