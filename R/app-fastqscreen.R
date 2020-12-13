###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodFastqScreen <- function(input = NA, output = NA, param = NA,
                                htmlFile = "00index.html") {

  if (ezIsSpecified(param$samples)){
    input <- input$subset(param$samples)
  }  
  
  # Override the pairing info and treat as single end to be faster
  param$paired = FALSE
  
  # Override the virus check parameter for human data
  if (grepl("^Human|^Homo", dataset$Species[1])) {
    param[["virusCheck"]] <- T
  }
  
  if (input$readType() == "bam") {
    stopifnot(input$getLength() == 1L) ## We only support one uBam now.
    input <- ezMethodBam2Fastq(
      input = input, param = param,
      OUTPUT_PER_RG = TRUE
    )
  }
  
  
  inputRaw <- ezMethodSubsampleFastq(input = input, param = param, n=param$nReads)
  param$trimAdapter <- TRUE
  inputProc <- ezMethodFastpTrim(input = inputRaw, param = param)

  ## map to adapters ----
  rawScreenResult = getFastqScreenStats(param,
                                        confFile = FASTQSCREEN_ADAPTER_CONF,
                                        inputRaw$getFullPaths("Read1"), workDir="rawReads")
  procScreenResult = getFastqScreenStats(param,
                                         confFile = FASTQSCREEN_GENOMICDNA_RIBORNA_CONF,
                                         inputProc$getFullPaths("Read1"), workDir="procReads")
  if (param[["virusCheck"]]) {
    unmappedFiles = gsub(".fastq.gz$", ".tagged_filter.fastq.gz", 
                         file.path("procReads", basename(inputProc$getFullPaths("Read1"))))
    names(unmappedFiles) =inputProc$getNames()
    virusResult <- map_and_count_virus(param, unmappedFiles, workDir="virusResult")
  }  
  refseqResult = map_and_count_refseq(param, inputProc$getFullPaths("Read1"), workDir="refseqResult", 
                       readCount = inputProc$getColumn("Read Count"))
  
  rRNAstrandResult <- get_rRNA_Strandness(param, inputProc)
  krakenResult <- runKraken(param, inputProc)

  file.remove(inputProc$getFullPaths("Read1"))

  makeRmdReport(
    output = output, param = param,
    rawScreenResult=rawScreenResult, procScreenResult=procScreenResult, virusResult=virusResult,
    rRNAstrandResult=rRNAstrandResult, krakenResult=krakenResult, refseqResult=refseqResult,
    rmdFile = "FastqScreen.Rmd", reportTitle = paste("Fastq Screen", param$name)
  )
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodFastqScreen(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
##' @section Functions:
##' \itemize{
##'   \item{\code{executeFastqscreenCMD(param, files): }}
##'   {Executes the fastqscreen command with given parameters on the input files.}
##'   \item{\code{collectFastqscreenOutput(dataset, files, resultFiles): }}
##'   {Collects the fastqscreen output after the result files have been obtained by \code{executeFastqscreenCMD()}.}
##'   \item{\code{executeBowtie2CMD(param, input): }}
##'   {Executes the bowtie2 command with given parameters on the input files.}
##'   \item{\code{collectBowtie2Output(param, dataset, countFiles): }}
##'   {Collects the bowtie2 output after the count files have been obtained by \code{executeBowtie2CMD()}.}
##'   \item{\code{fastqscreenReport(dataset, param, htmlFile="00index.html", fastqData, speciesPercentageTop): }}
##'   {Generates a report with plots and other information about the outcome of the run.}
##' }
EzAppFastqScreen <-
  setRefClass("EzAppFastqScreen",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodFastqScreen
        name <<- "EzAppFastqScreen"
        appDefaults <<- rbind(
          nTopSpecies = ezFrame(Type = "integer", DefaultValue = 10,
                                Description = "number of species to show in the plots"),
          confFile = ezFrame(Type = "character", DefaultValue = "",
                             Description = "the configuration file for fastq screen"),
          virusCheck = ezFrame(Type = "logical", DefaultValue = FALSE,
                               Description = "check for viruses in unmapped data"),
          minAlignmentScore = ezFrame(Type = "integer", DefaultValue = "-20",
                                      Description = "the min alignment score for bowtie2"),
          trimAdapter = ezFrame(Type = "logical", DefaultValue = TRUE,
                                Description = "whether to search for the adapters and trim them"),
          copyReadsLocally = ezFrame(Type = "logical", DefaultValue = FALSE,
                                     Description = "copy reads to scratch first")
        )
      }
    )
  )



getFastqScreenStats <- function(param, confFile = NULL, files, workDir="fastqScreenTmp") {
  dir.create(workDir)
  cmd <- paste(
    "fastq_screen", " --threads", param$cores, " --conf ", confFile,
    paste(files, collapse = " "), "--outdir", workDir, "--nohits --aligner bowtie2",
    "> fastqscreen.out", "2> fastqscreen.err"
  )
  ezSystem(cmd)

  fastqData <- list()
  # fastqData$PercentMapped <- setNames(double(length(files)), names(files))
  # fastqData$PercentMapped <- setNames(double(length(files)), names(files))
  # fastqData$Reads <- setNames(integer(length(files)), names(files))
  # fastqData$CommonResults <- setNames(vector("list", length(files)), names(files))
  for (nm in names(files)) {
    resultFile = str_replace(basename(files[nm]), "\\.gz$", "") %>%
      str_replace("\\.fastq$", "") %>%
      str_c("_screen.txt") ## remove the suffix .fastq[.gz] with _screen.txt
    #cat("Process ", files[nm], " - ", resultFiles[nm], " :")
    emptyLinePos = match("", readLines(file.path(workDir, resultFile)))
    x <- ezRead.table(file.path(workDir, resultFile),
                      skip = 1, stringsAsFactors = F,
                      nrows=emptyLinePos-3) ## first line and
    fastqData[[nm]] = x
    # fastqData$Reads[nm] <- x$"#Reads_processed"[1]
    # if (nrow(x) > 0) {
    #   UnmappedPercent <- as.numeric(unlist(strsplit(x$Genome[nrow(x)], split = " "))[2])
    #   fastqData$PercentMapped[nm] <- round((100 - UnmappedPercent), digits = 2)
    # } else {
    #   fastqData$MappingRate[nm] <- 0
    # }
    # rownames(x) <- x$Genome
    # x <- x[-nrow(x), grep("%.*hit", colnames(x))]
    # fastqData$CommonResults[[nm]] <- x ## can be a data.frame with 0 rows
    
    # noHit_files <- gsub(".fastq.gz$", ".tagged_filter.fastq.gz", files_ppData)
    # readCount <- ezFrame(
    #   totalReadCount = integer(length(files_ppData)),
    #   unmappedReadCount = integer(length(files_ppData)),
    #   row.names = names(files_ppData)
    # )
    # readCount[, "totalReadCount"] <- countReadsInFastq(files_ppData)
    # readCount[, "unmappedReadCount"] <- countReadsInFastq(noHit_files)
    
    
  }
  return(fastqData)
}


executeFastqscreenCMD <- function(param, confFile = NULL, files) {
  cmd <- paste(
    "fastq_screen", " --threads", param$cores, " --conf ", confFile,
    paste(files, collapse = " "), "--outdir . --nohits --aligner bowtie2",
    "> fastqscreen.out", "2> fastqscreen.err"
  )
  ezSystem(cmd)
  resultFiles <- str_replace(basename(files), "\\.gz$", "") %>%
    str_replace("\\.fastq$", "") %>%
    str_c("_screen.txt") ## remove the suffix .fastq[.gz] with _screen.txt
  names(resultFiles) <- names(files)
  return(resultFiles)
}

collectFastqscreenOutput <- function(files, resultFiles) {
  fastqData <- list()
  fastqData$MappingRate <- numeric()
  fastqData$Reads <- integer()
  fastqData$CommonResults <- list()
  for (nm in names(files)) {
    cat("Process ", files[nm], " - ", resultFiles[nm], " :")
    x <- ezRead.table(resultFiles[nm],
      skip = 1, stringsAsFactors = F,
      blank.lines.skip = T, fill = T, row.names = NULL
    )
    fastqData$Reads[nm] <- x$"#Reads_processed"[1]
    if (nrow(x) > 0) {
      UnmappedReads <- as.numeric(unlist(strsplit(x$Genome[nrow(x)], split = " "))[2])
      fastqData$MappingRate[nm] <- round((100 - UnmappedReads), digits = 2)
    } else {
      fastqData$MappingRate[nm] <- 0
    }
    rownames(x) <- x$Genome
    x <- x[-nrow(x), grep("%.*hit", colnames(x))]
    fastqData$CommonResults[[nm]] <- x ## can be a data.frame with 0 rows
  }
  return(fastqData)
}

executeBowtie2CMD <- function(param, input) {
  r1Files <- input$getFullPaths("Read1")
  if (param$paired) {
    r2Files <- input$getFullPaths("Read2")
  }
  countFiles <- character()
  for (nm in names(r1Files)) {
    countFiles[nm] <- paste0(nm, "-counts.txt")
    bowtie2options <- param$cmdOptions
    writeLines("ReadID\tRefSeqID\tAlignmentScore", countFiles[nm])
    if (!param$paired) {
      cmd <- paste(
        "bowtie2", "-x", REFSEQ_mRNA_REF,
        " -U ", r1Files[nm], bowtie2options, "-p", param$cores,
        "--no-unal --no-hd --mm", "2> ", paste0(nm, "_bowtie2.err"),
        "| cut -f1,3,12", " |sed s/AS:i://g", ">>", countFiles[nm]
      )
    } else {
      cmd <- paste(
        "bowtie2", "-x", REFSEQ_mRNA_REF,
        " -1 ", r1Files[nm], " -2 ", r2Files[nm], bowtie2options, "-p", param$cores,
        "--no-discordant --no-mixed --no-unal --no-hd --mm",
        "2> ", paste0(nm, "_bowtie2.err"),
        "| cut -f1,3,12", " |sed s/AS:i://g", ">>", countFiles[nm]
      )
    }
    ezSystem(cmd)
  }
  return(countFiles)
}

executeBowtie2CMD_Virus <- function(param, files) {
  r1Files <- files
  countFiles <- character()
  for (nm in names(r1Files)) {
    countFiles[nm] <- paste0(nm, "-counts.txt")
    bowtie2options <- param$cmdOptions
    writeLines("ReadID\tRefSeqID\tAlignmentScore", countFiles[nm])
    cmd <- paste(
      "bowtie2", "-x", REFSEQ_pathogenicHumanViruses_REF,
      " -U ", r1Files[nm], bowtie2options, "-p", param$cores,
      "--no-unal --no-hd", "2> ", paste0(nm, "_bowtie2.err"),
      "| cut -f1,3,12", " |sed s/AS:i://g", ">>", countFiles[nm]
    )
    ezSystem(cmd)
  }
  return(countFiles)
}

get_rRNA_Strandness <- function(param, input) {
  rRNA_REF <- "/srv/GT/reference/Silva/silva/release_123_1/SILVA_123.1_LSU_SSU"
  r1Files <- input$getFullPaths("Read1")

  countFiles <- character()
  for (nm in names(r1Files)) {
    countFiles[nm] <- paste0(nm, "-counts.txt")
    bowtie2options <- param$cmdOptions
    writeLines("ReadID\tFlag\tID\tAlignmentScore", countFiles[nm])

    cmd <- paste(
      "bowtie2", "-x", rRNA_REF,
      " -U ", r1Files[nm], bowtie2options, "-p", param$cores,
      "--no-unal --no-hd --mm", "2> ", paste0(nm, "_bowtie2.err"),
      "| cut -f1,2,3,12", " |sed s/AS:i://g", ">>", countFiles[nm]
    )
    ezSystem(cmd)
  }

  rRNA_strandInfo = input$meta[ , "Read Count", drop=FALSE]
  rRNA_strandInfo$Sense = 0
  rRNA_strandInfo$Antisense = 0
  for (nm in names(countFiles)) {
    countData <- ezRead.table(countFiles[nm], row.names = NULL)
    bestScores <- tapply(countData$AlignmentScore, countData$ReadID, max)
    countData <- countData[countData$AlignmentScore == bestScores[countData$ReadID], , drop = FALSE]
    countData <- countData[countData$AlignmentScore >= param$minAlignmentScore, , drop = FALSE]
    rRNA_strandInfo[nm, "Sense"] <- length(which(countData$Flag == 0))
    rRNA_strandInfo[nm, "Antisense"] <- length(which(countData$Flag == 16))
  }
  return(rRNA_strandInfo)
}

runKraken <- function(param, input) {
  r1Files <- input$getFullPaths("Read1")
  krakenResult <- list()
  for (i in 1:length(r1Files)) {
    resultFile <- paste0("report_", names(r1Files)[i], ".kraken")
    cmd <- paste("kraken2 --db", KRAKEN2_DB, r1Files[i],
                 "--gzip-compressed --threads", param$cores, "--report", 
                 resultFile, ">sequences.kraken")
    ezSystem(cmd)
    ## simple result filtering
    data <- ezRead.table(resultFile, row.names = NULL)
    colnames(data) <- c("readPercentage", "nreads_clade", "nreads_taxon", "rankCode", "ncbi", "name")
    data <- data[data$rankCode %in% c("U", "S"), ]
    data <- data[order(data$readPercentage, decreasing = TRUE), ]
    data <- data[!(data$ncbi %in% c(1, 131567, 136843)), ] # remove general terms

    ## save table for report
    ezWrite.table(data, resultFile, row.names = FALSE)
    data <- data[!(data$ncbi %in% c(0)), ]
    data <- data[1:min(nrow(data), 10), ]
    krakenResult[[i]] <- data
  }
  names(krakenResult) <- names(r1Files)
  file.remove("sequences.kraken")
  return(krakenResult)
}



map_and_count_virus <- function(param, files, workDir="virusResult", readCount = countReadsInFastq(files)) {
  dir.create(workDir)
  countFiles <- character()
  for (nm in names(files)) {
    countFiles[nm] <- paste0(workDir, "/", nm, "-counts.txt")
    bowtie2options <- param$cmdOptions
    writeLines("ReadID\tRefSeqID\tAlignmentScore", countFiles[nm])
    cmd <- paste(
      "bowtie2", "-x", REFSEQ_pathogenicHumanViruses_REF,
      " -U ", files[nm], bowtie2options, "-p", param$cores,
      "--no-unal --no-hd", "2> ", paste0(workDir, "/", nm, "_bowtie2.err"),
      "| cut -f1,3,12", " |sed s/AS:i://g", ">>", countFiles[nm]
    )
    ezSystem(cmd)
  }
  countList = collectBowtie2Output(param, countFiles, readCount, virusResult = T)
  return(countList)
}


map_and_count_refseq <- function(param, files, workDir="refseqResult", readCount = countReadsInFastq(files)) {
  dir.create(workDir)
  countFiles <- character()
  for (nm in names(files)) {
    countFiles[nm] <- paste0(workDir, "/", nm, "-counts.txt")
    bowtie2options <- param$cmdOptions
    writeLines("ReadID\tRefSeqID\tAlignmentScore", countFiles[nm])
    cmd <- paste(
      "bowtie2", "-x", REFSEQ_mRNA_REF,
      " -U ", files[nm], bowtie2options, "-p", param$cores,
      "--no-unal --no-hd", "2> ", paste0(workDir, "/", nm, "_bowtie2.err"),
      "| cut -f1,3,12", " |sed s/AS:i://g", ">>", countFiles[nm]
    )
    ezSystem(cmd)
  }
  countList = collectBowtie2Output(param, countFiles, readCount, virusResult = F)
  return(countList)
}





collectBowtie2Output <- function(param, countFiles, readCount, virusResult = F) {
  tax2name <- read.table(sub("Sequence.*", "Annotation/tax2name.txt", REFSEQ_mRNA_REF),
    header = F, stringsAsFactors = F, sep = "\t",
    colClasses = "character", quote = "", comment.char = ""
  )
  tax2name <- set_names(tax2name[, 3], tax2name[, 1])
  speciesPercentageTop <- list()
  for (nm in names(countFiles)) {
    countData <- ezRead.table(countFiles[nm], row.names = NULL)
    bestScores <- tapply(countData$AlignmentScore, countData$ReadID, max)
    countData <- countData[countData$AlignmentScore == bestScores[countData$ReadID], , drop = FALSE]
    countData <- countData[countData$AlignmentScore >= param$minAlignmentScore, , drop = FALSE]
    if (nrow(countData) > 0) {
      if (virusResult) {
        countData$species <- substr(sub("_NC_[0-9].*", "", countData$RefSeqID), 1, 30)
      } else {
        countData$species <- sub("_.*", "", countData$RefSeqID)
      }
      speciesHitsPerRead <- tapply(countData$species, countData$ReadID, unique)
      uniqSpeciesHitsPerRead <- names(speciesHitsPerRead)[lengths(speciesHitsPerRead) == 1]
      ### Result UniqHits:
      uniqSpeciesHits <- sort(table(unlist(speciesHitsPerRead[uniqSpeciesHitsPerRead])),
        decreasing = T
      )
      ### Results MultipleHits:
      multipleSpeciesHitsPerRead <- countData[!(countData$ReadID %in% uniqSpeciesHitsPerRead), ]
      by <- paste(multipleSpeciesHitsPerRead$ReadID,
        multipleSpeciesHitsPerRead$species,
        sep = "_"
      )
      ## multipleSpeciesHits = sort(table(multipleSpeciesHitsPerRead$species[!duplicated(by)]), decreasing=T) # is equivalent to the row below
      multipleSpeciesHits <- sort(table(tapply(multipleSpeciesHitsPerRead$species, by, unique)),
        decreasing = T
      )

      if (length(uniqSpeciesHits) > param$nTopSpecies) {
        topSpeciesUniq <- uniqSpeciesHits[1:param$nTopSpecies]
      } else {
        topSpeciesUniq <- uniqSpeciesHits
      }
      ## Special case where all hits are multi hits --- in that case we sort according to the multi-hits
      if (length(uniqSpeciesHits) == 0) {
        topSpeciesUniq <- integer(min(param$nTopSpecies, length(multipleSpeciesHits)))
        names(topSpeciesUniq) <- names(multipleSpeciesHits)[1:min(
          param$nTopSpecies,
          length(multipleSpeciesHits)
        )]
      }

      multipleSpeciesHits[setdiff(names(topSpeciesUniq), names(multipleSpeciesHits))] <- 0
      topSpeciesMultiple <- multipleSpeciesHits[names(topSpeciesUniq)]

      taxIds <- names(topSpeciesUniq)
      taxNames <- tax2name[taxIds]
      hasNoName <- is.na(taxNames)
      taxNames[hasNoName] <- taxIds[hasNoName]
      x <- ezFrame(
        UniqueSpeciesPercent = signif(100 * as.matrix(topSpeciesUniq) / readCount[nm], digits = 4),
        MultipleSpeciesPercent = signif(100 * as.matrix(topSpeciesMultiple) / readCount[nm], digits = 4),
        row.names = taxNames
      )
      speciesPercentageTop[[nm]] <- x
    } else {
      speciesPercentageTop[[nm]] <- NULL
    }
  }
  return(speciesPercentageTop)
}

makeScatterplot <- function(dataset, colname1, colname2) {
  if (colname1 %in% colnames(dataset) && colname2 %in% colnames(dataset) && nrow(dataset) > 1) {
    if (!all(dataset[[colname1]] == 0) && !any(is.na(dataset[[colname1]])) && !all(dataset[[colname2]] == 0) && !any(is.na(dataset[[colname2]]))) {
      ## LibCon column can all be 0 or NA. Then don't plot.
      dataset <- dataset[order(dataset[[colname1]], decreasing = T), ]
      corResult <- cor.test(dataset[[colname1]], dataset[[colname2]],
                            method = "spearman"
      )
      regressionResult <- lm(dataset[[colname2]] ~ dataset[[colname1]])
      label1 <- sub(" \\[.*", "", colname1)
      label2 <- sub(" \\[.*", "", colname2)
      
      ## plotly
      require(plotly)
      # a function to calculate your abline
      xmin <- min(dataset[[colname1]]) - 5
      xmax <- max(dataset[[colname1]]) + 5
      intercept <- regressionResult$coefficients[1]
      slope <- regressionResult$coefficients[2]
      p_abline <- function(x, a, b) {
        y <- a * x + b
        return(y)
      }

      p <- plot_ly(
        x = dataset[[colname1]], y = dataset[[colname2]],
        text = rownames(dataset)
      ) %>%
        add_markers() %>%
        add_text(textposition = "top right") %>%
        plotly::layout(showlegend = FALSE)
      a <- list(
        x = max(dataset[[colname1]]),
        y = max(dataset[[colname2]]),
        text = paste0("r=", round(corResult$estimate, 2)),
        xref = "x",
        yref = "y",
        showarrow = FALSE
      )
      p <- p %>% plotly::layout(
        shapes = list(
          type = "line", line = list(dash = "dash"),
          x0 = xmin, x1 = xmax,
          y0 = p_abline(xmin, slope, intercept),
          y1 = p_abline(xmax, slope, intercept)
        ),
        annotations = a,
        title = paste(label1, "vs", label2), yaxis = list(title = label2),
        xaxis = list(title = colname1)
      )
      return(p)
    }
  }
}

# dataset = ezRun::ezRead.table('/srv/gstore/projects/p2821/NovaSeq_20180829_NOV17_o4694/dataset.tsv')
# makeScatterplot(dataset, colname1 ='RIN [Factor]', colname2='Read Count')
# makeScatterplot(dataset, colname1 ='RIN [Factor]', colname2='LibConc_100_800bp [Characteristic]')
