###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

### EzAppSingleCellSTAR
EzAppSingleCellSTAR <-
  setRefClass("EzAppSingleCellSTAR",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodSingleCellSTAR
        name <<- "EzAppSingleCellSTAR"
        appDefaults <<- rbind(
          getJunctions = ezFrame(
            Type = "logical",
            DefaultValue = "FALSE",
            Description = "should junctions be returned"
          ),
          markDuplicates = ezFrame(
            Type = "logical",
            DefaultValue = "FALSE",
            Description = "should duplicates be marked with picard"
          ),
          twopassMode = ezFrame(
            Type = "logical",
            DefaultValue = "FALSE",
            Description = "1-pass mapping or basic 2-pass mapping"
          ),
          controlSeqs = ezFrame(
            Type = "charVector",
            DefaultValue = "",
            Description = "control sequences to add"
          )
        )
      }
    )
  )

### STAR for single cell data: reads in a unmapped bam
ezMethodSingleCellSTAR <- function(input = NA, output = NA, param = NA) {
  require(withr)
  refDir <- getSTARReference(param)
  bamFile <- output$getColumn("BAM")

  isSingleBam <- !is.na(input$readType()) && input$readType() == "bam"
  if (isSingleBam) {
    ## Read 1 is uBam
    fastqInput <- ezMethodBam2Fastq(
      input = input, param = param,
      OUTPUT_PER_RG = TRUE
    )
  } else {
    ## The read data is in CellDataset and input is fastq files
    fastqInput <- EzDataset(
      file = input$getFullPaths("CellDataset"),
      dataRoot = DEFAULT_DATA_ROOT
    )
  }

  trimmedInput <- ezMethodFastpTrim(input = fastqInput, param = param)

  ## Merge and clean prepross logs
  preprocessLogFns <- str_c(trimmedInput$getNames(), "_preprocessing.log")
  preprocessLogs <- map_chr(preprocessLogFns, read_lines)
  file.remove(preprocessLogFns)
  write_lines(preprocessLogs, file = str_c(input$getNames(), "_preprocessing.log"))


  # Clean converted fastqs
  if (isSingleBam) {
    file.remove(fastqInput$getFullPaths("Read1"))
    if (param$paired) {
      file.remove(fastqInput$getFullPaths("Read2"))
    }
  }

  ## fastq to bam
  trimmedBamFn <- tempfile(
    pattern = "trimmedBam", tmpdir = getwd(), fileext = ".bam"
  )
  if (param$paired) {
    fastqs2bam(
      fastqFns = trimmedInput$getFullPaths("Read1"),
      fastq2Fns = trimmedInput$getFullPaths("Read2"),
      readGroupNames = trimmedInput$getNames(),
      bamFn = trimmedBamFn, mc.cores = param$cores
    )
    file.remove(c(
      trimmedInput$getFullPaths("Read1"),
      trimmedInput$getFullPaths("Read2")
    ))
  } else {
    fastqs2bam(
      fastqFns = trimmedInput$getFullPaths("Read1"),
      readGroupNames = trimmedInput$getNames(),
      bamFn = trimmedBamFn, mc.cores = param$cores
    )
    file.remove(trimmedInput$getFullPaths("Read1"))
  }
  ## We can concatenate the fastqs for the aligner
  ## But it will takes more space. So we delete and convert here.

  inputTrimmed <- input$copy()
  inputTrimmed$setColumn("Read1", trimmedBamFn)
  inputTrimmed$dataRoot <- NULL
  alignerInput <- ezMethodBam2Fastq(
    input = inputTrimmed, param = param,
    OUTPUT_PER_RG = FALSE
  )

  if (param$cmdOptions == "") {
    param$cmdOptions <- "--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --outSAMattributes All"
  }

  if (!grepl("outSAMattributes", param$cmdOptions)) {
    param$cmdOptions <- paste(param$cmdOptions, "--outSAMattributes All")
  }

  genomeFn <- param$ezRef@refFastaFile

  if (ezIsSpecified(param$controlSeqs)) {
    ## control sequences
    controlSeqsLocalFn <- tempfile(
      pattern = "controlSeqs", tmpdir = getwd(), fileext = ".fa"
    )
    writeXStringSet(getControlSeqs(param$controlSeqs), filepath = controlSeqsLocalFn)
    defer(file.remove(controlSeqsLocalFn))

    genomeLocalFn <- tempfile(
      pattern = "genome", tmpdir = getwd(), fileext = ".fa"
    )
    file.copy(from = genomeFn, to = genomeLocalFn)
    writeXStringSet(getControlSeqs(param$controlSeqs),
      filepath = genomeLocalFn,
      append = TRUE
    )
    dictFile <- sub(".fa$", ".dict", genomeLocalFn)
    cmd <- str_c(
      prepareJavaTools("picard"), "CreateSequenceDictionary",
      str_c("R=", genomeLocalFn), str_c("O=", dictFile),
      sep = " "
    )
    ezSystem(cmd)

    genomeFn <- genomeLocalFn
    defer(file.remove(c(genomeLocalFn, dictFile)))
  }

  cmd <- str_c(
    "STAR", " --genomeDir", refDir, "--sjdbOverhang 150",
    "--readFilesIn", alignerInput$getColumn("Read1"),
    if_else(param$paired, alignerInput$getColumn("Read2"), ""),
    "--twopassMode", if_else(param$twopassMode, "Basic", "None"),
    if_else(ezIsSpecified(param$controlSeqs),
            str_c("--genomeFastaFiles", controlSeqsLocalFn, sep = " "), ""),
    "--runThreadN", param$cores, param$cmdOptions,
    "--outStd BAM_Unsorted --outSAMtype BAM Unsorted", ">  Aligned.out.bam",
    sep = " "
  )
  ezSystem(cmd)
  file.remove(alignerInput$getColumn("Read1"))
  if (param$paired) {
    file.remove(alignerInput$getColumn("Read2"))
  }

  ## clean star log files
  defer(file.remove(c("Log.progress.out", "Log.out", "Log.std.out")))
  defer(unlink(c("_STARgenome", "_STARpass1"), recursive = TRUE, force = TRUE))

  ## Merge unmapped and mapped bam to recover the tags
  mergeBamAlignments(
    alignedBamFn = "Aligned.out.bam",
    unmappedBamFn = inputTrimmed$getFullPaths("Read1"),
    outputBamFn = "Aligned.out.merged.bam",
    fastaFn = genomeFn
  )
  file.remove("Aligned.out.bam")
  file.rename(from = "Aligned.out.merged.bam", to = "Aligned.out.bam")
  file.remove(trimmedBamFn)

  nSortThreads <- min(param$cores, 8)

  file.rename("Log.final.out", to = basename(output$getColumn("STARLog")))

  if (!is.null(param$markDuplicates) && param$markDuplicates) {
    ezSortIndexBam("Aligned.out.bam", "sorted.bam",
      ram = param$ram, removeBam = TRUE, cores = nSortThreads
    )
    dupBam(inBam = "sorted.bam", outBam = basename(bamFile), operation = "mark")
    file.remove("sorted.bam")
  } else {
    ezSortIndexBam("Aligned.out.bam", basename(bamFile),
      ram = param$ram, removeBam = TRUE, cores = nSortThreads
    )
  }

  if (param$getJunctions) {
    file.rename(from = "SJ.out.tab", to = basename(output$getColumn("Junctions")))
    file.rename(from = "Chimeric.out.junction", to = basename(output$getColumn("Chimerics")))
  } else {
    defer(file.remove(c("SJ.out.tab", "Chimeric.out.junction", "Chimeric.out.sam")))
  }

  ## check the strandedness
  ezSystem(str_c(
    "infer_experiment.py", "-r", getReferenceFeaturesBed(param),
    "-i", basename(bamFile), "-s 1000000",
    sep = " "
  ))

  return("Success")
}
