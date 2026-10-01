###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## Human contamination is always removed before removing the (non-human)
## host, since samples can pick up human reads at the bench regardless of
## which species the sample itself comes from. If the user-selected host
## build is itself human, this step is skipped since it would just repeat
## the same mapping twice. Shared by EzAppVirDetect and EzAppDeHost.
DEFAULT_HUMAN_REFBUILD <- "Homo_sapiens/GENCODE/GRCh38.p14"
## VirDetect filters human reads against the telomere-to-telomere assembly,
## which also covers centromeres/satellites/rDNA arrays missing from GRCh38
## and thus removes more human reads before the virome mapping
VIRDETECT_HUMAN_REFBUILD <- "Homo_sapiens/Ensembl/T2T-CHM13v2.0"

##' @title Counts the reads in a (gzipped) fastq file
countFastqReads <- function(fastqFile) {
  nLines <- as.integer(ezSystem(
    paste("zcat", fastqFile, "| wc -l"),
    intern = TRUE
  ))
  nLines %/% 4L
}

##' @title Depletes reads mapping to one reference genome
##' @description Maps read1/read2 to \code{refBuild} with bowtie2 and keeps
##' only the reads/pairs for which nothing mapped, writing them out as
##' gzipped fastq. Returns the kept fastq paths and the number of kept
##' reads/pairs. bowtie2's own log is appended to \code{logFile}.
depleteAgainstReference <- function(
  read1,
  read2 = NULL,
  param,
  refBuild,
  cmdOptions,
  stageLabel,
  logFile
) {
  paramRef <- param
  paramRef$refBuild <- refBuild
  paramRef$ezRef <- EzRef(paramRef)
  ref <- getBowtie2Reference(paramRef)

  ## stream bowtie2 -> filter unmapped -> fastq, so that no intermediate
  ## BAM or collate temp files hit the scratch disk; bowtie2 emits both
  ## mates of a pair adjacently, so no collate is needed
  r1Out <- paste0(stageLabel, "_R1.fastq.gz")
  r2Out <- if (param$paired) paste0(stageLabel, "_R2.fastq.gz") else NULL
  cmd <- paste(
    "bowtie2",
    cmdOptions,
    "-p",
    param$cores,
    "-x",
    ref,
    if (param$paired) {
      paste("-1", read1, "-2", read2)
    } else {
      paste("-U", read1)
    },
    "2>>",
    logFile,
    "|",
    "samtools view -u",
    if (param$paired) "-f 12 -F 256" else "-f 4 -F 256",
    "-",
    "|",
    "samtools fastq",
    if (param$paired) {
      paste("-1", r1Out, "-2", r2Out, "-0 /dev/null -s /dev/null -")
    } else {
      paste("-0", r1Out, "-")
    }
  )
  ezSystem(cmd)
  nKept <- countFastqReads(r1Out)

  list(read1 = r1Out, read2 = r2Out, count = nKept)
}

ezMethodVirDetect <- function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  sampleName <- input$getNames()
  setwdNew(sampleName)
  ## keep temp files of child processes (e.g. samtools) on the job's scratch
  ## dir instead of the node's small /tmp
  oldTmpDir <- Sys.getenv("TMPDIR", unset = NA)
  Sys.setenv(TMPDIR = getwd())
  on.exit(
    {
      if (is.na(oldTmpDir)) {
        Sys.unsetenv("TMPDIR")
      } else {
        Sys.setenv(TMPDIR = oldTmpDir)
      }
    },
    add = TRUE
  )

  ## trim reads
  param$fastpCompression <- 9
  trimmedInput <- ezMethodFastpTrim(input = input, param = param)

  countSummary <- list(
    QCReads = data.frame(
      Stage = "QCReads",
      Count = countFastqReads(trimmedInput$getColumn("Read1"))
    )
  )

  ## build the depletion chain: always remove human contamination first,
  ## unless the user-selected host build already is human
  refChain <- c(Host = param$hostBuild)
  if (!grepl("^Homo_sapiens", param$hostBuild)) {
    refChain <- c(Human = VIRDETECT_HUMAN_REFBUILD, refChain)
  }

  curR1 <- trimmedInput$getColumn("Read1")
  curR2 <- if (param$paired) trimmedInput$getColumn("Read2") else NULL

  for (stageName in names(refChain)) {
    res <- depleteAgainstReference(
      read1 = curR1,
      read2 = curR2,
      param = param,
      refBuild = refChain[[stageName]],
      cmdOptions = param$cmdOptionsHost,
      stageLabel = stageName,
      logFile = paste0(stageName, "_bowtie2.log")
    )
    file.remove(curR1)
    if (param$paired) {
      file.remove(curR2)
    }
    curR1 <- res$read1
    curR2 <- res$read2
    countSummary[[stageName]] <- data.frame(
      Stage = paste0(stageName, "Removed"),
      Count = res$count
    )
  }

  ## align filtered reads to the viral reference database, get sorted bam
  ## file and index, output idxstats into a text file
  paramVirom <- param
  paramVirom$refBuild <- param$virBuild
  paramVirom$ezRef <- EzRef(paramVirom)
  vir <- getBowtie2Reference(paramVirom)
  readGroupOpt <- paste0(
    "--rg-id ",
    sampleName,
    " --rg SM:",
    sampleName,
    " --rg LB:RGLB_",
    sampleName,
    " --rg PL:illumina",
    " --rg PU:RGPU_",
    sampleName
  )
  viromeLog <- "Virome_bowtie2.log"
  cmd <- paste(
    "bowtie2",
    param$cmdOptions,
    "-p",
    param$cores,
    readGroupOpt,
    "-x",
    vir,
    if (param$paired) {
      paste("-1", curR1, "-2", curR2)
    } else {
      paste("-U", curR1)
    },
    "2>",
    viromeLog,
    "|",
    "samtools view -b -F 4 -",
    "> virome.bam"
  )
  ezSystem(cmd)
  file.remove(curR1)
  if (param$paired) {
    file.remove(curR2)
  }

  ezSortIndexBam(
    "virome.bam",
    "virome.sorted.bam",
    ram = param$ram * 0.7, ## put additional safety margin in
    removeBam = TRUE,
    cores = param$cores
  )
  ezSystem("samtools idxstats virome.sorted.bam > virome.idxstats.txt")
  bamFile <- "virome.sorted.bam"

  ## unique/multi mapped read counts, matched by content (not line position)
  ## so they don't depend on how many lines bowtie2 happens to print for this
  ## paired/single-end + cmdOptions combination (e.g. --no-mixed --no-discordant
  ## shortens the paired summary compared to the host-depletion runs above)
  viromeLogLines <- readLines(viromeLog)
  uniqLine <- grep(
    "aligned (concordantly )?exactly 1 time",
    viromeLogLines,
    value = TRUE
  )
  multiLine <- grep(
    "aligned (concordantly )?>1 times",
    viromeLogLines,
    value = TRUE
  )
  nUniqueMapped <- if (length(uniqLine) > 0) {
    as.integer(sub("^\\s*([0-9]+).*", "\\1", uniqLine[1]))
  } else {
    NA_integer_
  }
  nMultiMapped <- if (length(multiLine) > 0) {
    as.integer(sub("^\\s*([0-9]+).*", "\\1", multiLine[1]))
  } else {
    NA_integer_
  }
  countSummary[["UniqueMapped"]] <- data.frame(
    Stage = "UniqueMapped",
    Count = nUniqueMapped
  )
  countSummary[["MultiMapped"]] <- data.frame(
    Stage = "MultiMapped",
    Count = nMultiMapped
  )

  countSummary <- do.call(rbind, countSummary)
  rownames(countSummary) <- NULL
  ezWrite.table(countSummary, file = "read_count_summary.tsv", row.names = FALSE)

  ## collect summary statistics and save in a summary table, collect per base
  ## coverage of each mapped viral genome and save in individual csv files
  idx <- read.table(
    "virome.idxstats.txt",
    header = FALSE,
    stringsAsFactors = FALSE,
    colClasses = c("character", "integer", "integer", "integer")
  )
  sub <- idx[idx$V3 > 0, ]
  csvFile <- sub(".fa$", ".csv", paramVirom$ezRef["refFastaFile"])
  names <- read.csv(
    csvFile,
    quote = "",
    stringsAsFactors = FALSE,
    header = FALSE,
    colClasses = "character"
  )
  if (any(duplicated(names$V1))) {
    warning(
      "duplicate accessions found in the viral reference name table (",
      csvFile,
      "); keeping only the first occurrence of each to avoid double-counting hits"
    )
    names <- names[!duplicated(names$V1), ]
  }
  sub <- merge(sub, names, by = "V1")
  if (nrow(sub) != 0) {
    ## these three metrics are computed by this loop, so (unlike the merge()
    ## result above) we fully control their names -- assign by name instead
    ## of by numeric column position
    sub$mappedBases <- NA_real_
    sub$genomeCov_pect <- NA_real_
    sub$aveDepth <- NA_real_
    for (i in seq_len(nrow(sub))) {
      tryCatch(
        {
          chr <- sub[i, 1]
          len <- sub[i, 2]
          common_name <- sub[i, 6]
          temp.df <- data.frame(
            c1 = c(chr),
            c2 = c("0"),
            c3 = c(len),
            c4 = c(common_name)
          )
          bed.file <- paste0(chr, ".bed")
          csv.file <- paste0(chr, ".csv")
          write.table(
            temp.df,
            file = bed.file,
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t"
          )
          ezSystem(paste0(
            "samtools view -b ",
            bamFile,
            " ",
            chr,
            " > ",
            chr,
            ".bam"
          ))
          ezSystem(paste0("samtools index ", chr, ".bam"))
          subsam <- 1
          if (sub[i, 3] > 1000000) {
            subsam <- 1000000 / sub[i, 3]
            ezSystem(paste0(
              "samtools view -b -O BAM -o ",
              chr,
              ".subsam.bam",
              " -s ",
              subsam,
              " ",
              chr,
              ".bam"
            ))
            ezSystem(paste0("samtools index ", chr, ".subsam.bam"))
            ezSystem(paste0("mv ", chr, ".subsam.bam", " ", chr, ".bam"))
            ezSystem(paste0("mv ", chr, ".subsam.bam.bai", " ", chr, ".bam.bai"))
          }
          ezSystem(paste0(
            "bedtools coverage -sorted -a ",
            bed.file,
            " -b ",
            chr,
            ".bam",
            " -d > ",
            csv.file
          ))
          cov <- read.table(
            csv.file,
            header = FALSE,
            sep = "\t",
            quote = "",
            stringsAsFactors = FALSE
          )
          sub$mappedBases[i] <- sum(cov$V6 != 0)
          sub$genomeCov_pect[i] <- sum(cov$V6 != 0) / len * 100
          sub$aveDepth[i] <- sum(cov$V6) / len / subsam
          ezSystem(paste0("rm ", chr, ".bam"))
          ezSystem(paste0("rm ", chr, ".bam.bai"))
        },
        error = function(e) {
          ## keep going for the other detected viruses even if coverage
          ## computation fails for this one (e.g. a corrupt/truncated
          ## per-chromosome BAM, or a transient I/O error)
          message(
            "VirDetect: coverage computation failed for ",
            sub[i, 1],
            ": ",
            conditionMessage(e)
          )
        }
      )
    }
    nFailed <- sum(is.na(sub$aveDepth))
    if (nFailed > 0) {
      message(
        "VirDetect: excluding ",
        nFailed,
        " detected virus(es) from the summary table due to failed coverage computation"
      )
    }
    sub <- sub[!is.na(sub$aveDepth), ]
    sub <- sub[order(sub$aveDepth, decreasing = TRUE), ]
    out <- sub[, c(1, 6, 7, 2, 3, 8, 9, 10)]
    colnames(out) <- c(
      "ID",
      "CommonName",
      "Family",
      "Len",
      "mappedReads",
      "mappedBases",
      "genomeCov_pect",
      "aveDepth"
    )
    write.table(
      out,
      file = "summary_table.tsv",
      col.names = TRUE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
  }
  ## delete intermediate result files; the per-sample folder gets copied back
  ## to gstore, so only what the report needs (and per-stage logs, the
  ## summary table, and the per-virus coverage csv files, for provenance)
  ## is left behind
  ezSystem("rm -f *.bed")

  ## html file
  htmlFile <- output$getColumn("OutReport")
  styleFiles <- file.path(
    system.file("templates", package = "ezRun"),
    c("fgcz.css", "VirDetect.Rmd", "fgcz_header.html", "banner.png")
  )
  file.copy(from = styleFiles, to = ".", overwrite = TRUE)
  params <- list(
    sample = sampleName,
    minReadCount = param$minReadCount,
    hostBuild = param$hostBuild,
    virBuild = param$virBuild,
    paired = param$paired
  )
  rmarkdown::render(
    input = "VirDetect.Rmd",
    envir = new.env(),
    output_dir = ".",
    output_file = htmlFile,
    quiet = TRUE
  )
  cdir <- "."
  csv.files <- file.path(cdir, "*.csv")
  ezSystem(paste0("rm ", csv.files))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run
EzAppVirDetect <-
  setRefClass(
    "EzAppVirDetect",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodVirDetect
        name <<- "EzAppVirDetect"
        appDefaults <<- rbind(
          minReadCount = ezFrame(
            Type = "integer",
            DefaultValue = "9",
            Description = "use for reporting only viral genomes with mapped reads higher than"
          )
        )
      }
    )
  )
