###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCrisprScreenQC <- function(input, output, param) {
  require(Herper)
  require(Biostrings)
  require(seqLogo)
  require(ShortRead)

  local_CondaEnv("gi_mageck", pathToMiniConda = "/usr/local/ngseq/miniforge3")
  sgRNADirs <- list.dirs(param$libPath, full.names = TRUE)
  reportDir <- basename(output$getColumn("Report"))
  dir.create(reportDir)
  mergedRef <- c()
  for (i in 1:length(sgRNADirs)) {
    refFile <- list.files(
      sgRNADirs[i],
      pattern = '*_MAGeCK.csv',
      full.names = TRUE
    )
    if (length(refFile) == 1) {
      myRef <- read.csv(
        refFile,
        quote = '',
        header = FALSE,
        stringsAsFactors = FALSE
      )
      refName <- basename(dirname(refFile))
      myRef$V1 <- paste(refName, myRef$V1, sep = '--')
      mergedRef <- rbind(mergedRef, myRef)
    }
  }
  refFile <- 'GEML_Mageck_mergedRefs.csv'
  write.csv(mergedRef, refFile, quote = FALSE, row.names = FALSE)

  ##Subsample/fastp
  inputRaw <- ezMethodSubsampleFastq(
    input = input,
    param = param,
    n = param$nReads
  )
  param$trimAdapter <- TRUE
  nReads <- param$nReads
  param$nReads <- 0 #prevent second round of subsampling
  inputProc <- ezMethodFastpTrim(input = inputRaw, param = param)
  param$nReads <- nReads

  sampleNames <- inputProc$getNames()
  samplesToRemove <- c()
  inputFiles <- inputProc$getFullPaths("Read1")
  countsPerLib <- list()
  topFeatureResults <- list()

  PWMs <- list()

  for (i in 1:length(sampleNames)) {
    system2(
      "mageck",
      args = c(
        "count",
        "-l",
        refFile,
        "--fastq",
        inputFiles[i],
        "-n",
        sampleNames[i]
      )
    )
    resultFile <- paste0(sampleNames[i], '.count.txt')
    if (file.exists(resultFile)) {
      counts <- ezRead.table(resultFile, row.names = NULL)
      counts[['Lib']] = sub('--.*', '', counts$sgRNA)
      countsPerLib[[i]] <- tapply(counts$sample1, INDEX = counts$Lib, FUN = sum)
      topFeatureResults[[i]] <- counts[
        order(counts$sample1, decreasing = TRUE),
      ][1:param$topFeatures, ]
    } else {
      samplesToRemove <- c(samplesToRemove, sampleNames[i])
    }
    fq <- readFastq(inputFiles[i])
    reads <- sread(fq)
    consMatrix <- consensusMatrix(reads, as.prob = TRUE)
    consMatrix <- consMatrix[rownames(consMatrix) %in% c('A', 'C', 'G', 'T'), ]
    PWMs[[i]] <- makePWM(consMatrix)
  }
  names(PWMs) <- sampleNames
  if (length(samplesToRemove) > 0) {
    sampleNames <- sampleNames[!(sampleNames %in% samplesToRemove)]
  }
  names(topFeatureResults) <- sampleNames
  names(countsPerLib) <- sampleNames
  sgRNAPerLib <- table(counts[['Lib']])

  data <- list(
    sgRNAPerLib = sgRNAPerLib,
    countsPerLib = countsPerLib,
    topFeatureResults = topFeatureResults,
    PWMs = PWMs
  )

  setwd(reportDir)
  makeRmdReport(
    output = output,
    param = param,
    input = input,
    data = data,
    rmdFile = "CrisprScreenQC.Rmd",
    reportTitle = paste("CRISPR Screen QC", param$name)
  )
  return("Success")

  # Zac mail 7.10.22: The app could do something like
  # 1. detect and report the spacer sequences (seqLogo)
  # 2. map against the 5'-sequence and report mapping stats (seqLogo)
  # 3. remove 5'- and 3'-sequences (automatically done by MAGECK)
  # 4. map guide RNA against our list of guide library files (minimum the standard libraries and all previously processed libraries, we will make a depository) (MAGECK count all installed refs for MAGECK)
  # 5. map against REFseq mRNA??? (and other typical sources of contamination) (to be done, needs PWM based identification of Spacers)
  # 6. If we switch to UMIs, should the app take care of the relevant analysis? (to be done)
}

##' @template app-template
##' @templateVar method ezMethodCrisprScreenQC(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppCrisprScreenQC <-
  setRefClass(
    "EzAppCrisprScreenQC",
    contains = "EzApp",
    methods = list(
      ## mageck count, fastp (trimAdapter forced TRUE), ShortRead/Biostrings/seqLogo
      ## (base-composition PWM, same pattern as FastqScreenApp) all unconditional.
      citation = function() {
        c(
          "Li, W. et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology 15, 554 (2014). https://doi.org/10.1186/s13059-014-0554-4",
          "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560",
          "Morgan, M., Anders, S., Lawrence, M., Aboyoun, P., Pagès, H. & Gentleman, R. ShortRead: a bioconductor package for input, quality assessment and exploration of high-throughput sequence data. Bioinformatics 25(19), 2607-2608 (2009). https://doi.org/10.1093/bioinformatics/btp450",
          "Pagès, H., Aboyoun, P., Gentleman, R. & DebRoy, S. Biostrings: Efficient manipulation of biological strings. R package version 2.80.1. https://doi.org/10.18129/B9.bioc.Biostrings",
          "Bembom, O. & Ivanek, R. seqLogo: Sequence logos for DNA sequence alignments. R package version 1.78.0. https://doi.org/10.18129/B9.bioc.seqLogo"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodCrisprScreenQC
        name <<- "CrisprScreenQC"
        appDefaults <<- rbind(
          libPath = ezFrame(
            Type = "character",
            DefaultValue = "/srv/GT/databases/GEML/sgRNA_Libs",
            Description = "sgRNA Library Path"
          ),
          topFeatures = ezFrame(
            Type = "numeric",
            DefaultValue = 10,
            Description = "list top sgRNA per sample"
          )
        )
      }
    )
  )
