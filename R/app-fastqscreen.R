###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodFastqScreen <- function(input = NA, output = NA, param = NA,
                                htmlFile = "00index.html") {

  require(seqLogo)
  require(ShortRead)
  
  # Override the virus check parameter for human data
  if ("Species" %in% input$colNames && grepl("^Human|^Homo", input$getColumn("Species")[1])) {
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
  param$nReads <- 0 #prevent second round of subsampling
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
  }  else {
    virusResult = ''
  }
  refseqResult = map_and_count_refseq(param, inputProc$getFullPaths("Read1"), workDir="refseqResult", 
                       readCount = inputProc$getColumn("Read Count"))
  
  rRNAstrandResult <- get_rRNA_Strandness(param, inputProc)
  krakenResult <- runKraken(param, inputProc)
  
  sampleNames <- inputProc$getNames()
  PWMs <- list()
  PWMs[['R1']] <- list()
  inputFiles_R1  <- inputProc$getFullPaths("Read1")
  for (i in 1:length(sampleNames)){
      fq <- readFastq(inputFiles_R1[i])
      reads <- sread(fq)
      consMatrix <- consensusMatrix(reads, as.prob = TRUE)
      consMatrix <- consMatrix[rownames(consMatrix) %in% c('A', 'C','G','T'),]
      PWMs[['R1']][[i]] <- makePWM(consMatrix)
  }
  names(PWMs[['R1']]) <- sampleNames
  
if(param$paired){
  PWMs[['R2']] <- list()
  inputFiles_R2  <- inputProc$getFullPaths("Read2")
  for (i in 1:length(sampleNames)){
      fq <- readFastq(inputFiles_R2[i])
      reads <- sread(fq)
      consMatrix <- consensusMatrix(reads, as.prob = TRUE)
      consMatrix <- consMatrix[rownames(consMatrix) %in% c('A', 'C','G','T'),]
      PWMs[['R2']][[i]] <- makePWM(consMatrix)
  }
  names(PWMs[['R2']]) <- sampleNames
}
  
  file.remove(inputProc$getFullPaths("Read1"))
  file.remove(inputRaw$getFullPaths("Read1"))
  if(param$paired){  
    file.remove(inputProc$getFullPaths("Read2"))
    file.remove(inputRaw$getFullPaths("Read2"))
  }
  
  setwdNew(basename(output$getColumn("Report")))
  makeRmdReport(
    output = output, param = param,
    input=input,
    rawScreenResult=rawScreenResult, procScreenResult=procScreenResult, virusResult=virusResult,
    rRNAstrandResult=rRNAstrandResult, krakenResult=krakenResult, refseqResult=refseqResult,
    PWMs = PWMs,
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
  
  if(length(files) > 384){
      cat(cmd, file = 'fastqScreenCall.sh')
      result <- ezSystem('bash fastqScreenCall.sh')
  } else {
      result <- ezSystem(cmd)
  }
  gc()

  fastqData <- list()
  for (nm in names(files)) {
    resultFile = str_replace(basename(files[nm]), "\\.gz$", "") %>%
      str_replace("\\.fastq$", "") %>%
      str_c("_screen.txt") ## remove the suffix .fastq[.gz] with _screen.txt
    emptyLinePos = match("", readLines(file.path(workDir, resultFile)))
    if (is.na(emptyLinePos)){ ## no hits found!
      emptyLinePos = 3
    }
    x <- ezRead.table(file.path(workDir, resultFile),
                      skip = 1, stringsAsFactors = F,
                      nrows=emptyLinePos-3) ## first line and
    fastqData[[nm]] = x
  }
  return(fastqData)
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
  for (nm in names(r1Files)) {
    resultFile <- paste0("report_", nm, ".kraken")
    cmd <- paste("kraken2 --db", KRAKEN2_DB, r1Files[nm],
                 "--gzip-compressed --threads", param$cores, "--report", 
                 resultFile, ">sequences.kraken")
    ezSystem(cmd)
    if (length(readLines(resultFile)) > 0){
      data <- ezRead.table(resultFile, row.names = NULL)
      colnames(data) <- c("readPercentage", "nreads_clade", "nreads_taxon", "rankCode", "ncbi", "name")
      ## sort and filter results -- could also be done later
      data <- data[data$rankCode %in% c("U", "S"), ]
      data <- data[order(data$readPercentage, decreasing = TRUE), ]
      termsToRemove = c(0, 1, 131567, 136843)
      data <- data[!(data$ncbi %in% termsToRemove), ] # remove general terms
      data$name <- sub('^ *', '',data$name)
      data <- data[1:min(10, nrow(data)),] # Keep max top10 hits
    } else {
      data = ezFrame("readPercentage"=numeric(0), "nreads_clade"=numeric(0), "nreads_taxon"=numeric(0), 
                     "rankCode"=character(0), "ncbi"=character(0), "name"=numeric(0))
    }
    krakenResult[[nm]] <- data
  }
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
  bamFile <- file.path(workDir, "allSamples_unmapped.bam")
  fastqs2bam(fastqFns = files, readGroupNames = names(files), bamFn = bamFile)
  readToReadGroupFile <- file.path(workDir, "readID2RG.txt") 
  ezSystem(paste('samtools view', bamFile, '|cut -f 1,12| sed s/RG:Z:// >', readToReadGroupFile))
  bowtie2options <- param$cmdOptions
  inputFastq <- file.path(workDir, "allSamples.fastq")
  ezSystem(paste('samtools bam2fq', bamFile, '>', inputFastq))
  outputCountFile <-  file.path(workDir, "refSeq_Counts.txt")
  cmd = paste('bowtie2 -x', REFSEQ_mRNA_REF, '-U', inputFastq, bowtie2options, "-p", param$cores, '-t --no-unal', '2> ', paste0(workDir, '/bowtie2.err'), '| grep ^@ -v|cut -f1,3,12', '|sed s/AS:i://g >', outputCountFile)
  tryCatch(expr = {ezSystem(cmd)}, 
           error = function(e) {message('No reads aligned to RefSeq')})
  gc()
  completeOutputCountFile <-  file.path(workDir, "refSeq_Counts_allSamples.txt")
  if(file.exists(outputCountFile) && file.size(outputCountFile) > 0){
    ezSystem(paste('sort -k1', outputCountFile, '>', completeOutputCountFile))
    system(paste('join -j 1 -o 1.1,1.2,1.3,2.2', completeOutputCountFile, readToReadGroupFile,'>', outputCountFile)) #ezSystem thinks that it fails
    file.remove(completeOutputCountFile, bamFile, inputFastq, readToReadGroupFile)
  
    countFiles <- character()
    for (nm in names(files)) {
        countFiles[nm] <- paste0(workDir, "/", nm, "-counts.txt")
        writeLines("ReadID\tRefSeqID\tAlignmentScore\tName", countFiles[nm])
        system(paste('grep', paste0('\" ', nm, '$\"'), outputCountFile, '|sed s/[[:blank:]]/\\\t/g >>', countFiles[nm]))
    }
    countList = collectBowtie2Output(param, countFiles, readCount, virusResult = FALSE)
    return(countList)
  } else {
      return(NULL)
  }
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
