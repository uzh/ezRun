###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCountSpacer = function(input = NA, output = NA, param = NA) {
  require(ShortRead)
  require(htmlwidgets)
  require(stringi)

  sampleName = input$getNames()
  setwdNew(sampleName)
  csvFiles = list.files(
    file.path('/srv/GT/databases/GEML/sgRNA_Libs/', param[['dictPath']]),
    pattern = 'final.csv$',
    full.names = TRUE
  )
  if (length(csvFiles) < 1) {
    csvFiles <- list.files(
      file.path('/srv/GT/databases/GEML/sgRNA_Libs/', param[['dictPath']]),
      pattern = '.csv$',
      full.names = TRUE
    )
    param[['dictPath']] = csvFiles[grep('MAGeCK', csvFiles, invert = TRUE)]
  } else {
    param[['dictPath']] = csvFiles[1]
  }

  dict = ezRead.table(
    param[['dictPath']],
    header = FALSE,
    sep = ',',
    row.names = NULL
  )
  colnames(dict) = c('TargetID', 'Sequence', 'GeneSymbol', 'isControl')
  dict[['ID']] = paste(dict$TargetID, dict$Sequence, sep = '-')
  stats = list()
  stats[['rawReads']] = as.numeric(input$meta[['Read Count']])

  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  readFile = trimmedInput$getColumn("Read1")
  stats[['filteredReads']] = as.numeric(ezSystem(
    paste('zcat', readFile, '|wc -l'),
    intern = TRUE
  )) /
    4
  ## Spacer length: use the library's own sgRNA length when the dict is uniform
  ## (the norm). Anchoring the extraction to this fixed length -- rather than
  ## taking everything between the read start and the right pattern -- removes the
  ## constant 5' base (the U6 'G') and any staggered filler that otherwise leaves
  ## the spacer one or more bases too long and unalignable.
  spacerLength <- as.integer(param$spacerLength)
  if (is.na(spacerLength) || spacerLength <= 0) {
    seqLengths <- unique(nchar(dict$Sequence))
    spacerLength <- if (length(seqLengths) == 1) seqLengths else 0L
  }

  ## Read base-composition PWM (for the read-structure logo) and, when enabled,
  ## inferred flanking patterns. Computed on the trimmed reads before extraction.
  patternGuess <- guessFlankingPatterns(readFile)
  leftUsed <- param$leftPattern
  rightUsed <- param$rightPattern
  guessPatterns <- isTRUE(as.logical(param$guessPatterns))
  if (guessPatterns) {
    if (leftUsed == '') leftUsed <- patternGuess$leftGuess
    if (rightUsed == '') rightUsed <- patternGuess$rightGuess
  }
  patternInfo <- list(
    leftUsed = leftUsed,
    rightUsed = rightUsed,
    leftGuess = patternGuess$leftGuess,
    rightGuess = patternGuess$rightGuess,
    spacerLength = spacerLength,
    guessPatterns = guessPatterns
  )

  reads <- twoPatternReadFilter(
    readFile,
    leftUsed,
    rightUsed,
    param$maxMismatch,
    spacerLength = spacerLength
  )

  ###Export as fasta file
  readFile = paste0(sampleName, '.fa')
  reads = reads[width(reads) >= param$minReadLength]
  stats[['validSpacerReads']] = length(reads)
  writeXStringSet(reads, readFile, format = 'fasta')
  remove(reads)
  gc()

  resultFile = paste0(sampleName, '_bowtie.txt')
  bowtieIndex <- sub(
    '\\.[0-9]\\.ebwt',
    '',
    list.files(
      dirname(param[['dictPath']]),
      pattern = 'ebwt$',
      full.names = TRUE
    )[1]
  )
  cmd = paste(
    'bowtie',
    bowtieIndex,
    readFile,
    '-f -p',
    param$cores,
    '|cut -f3,8|sort >',
    resultFile
  )
  ezSystem(cmd)
  if (file.size(resultFile) == 0) {
    stop(
      'No read were aligned to your reference. Please check the correctness of your flanking patterns and the choice of the Crispr reference.'
    )
  }
  result_bowtie = ezRead.table(resultFile, row.names = NULL, header = FALSE)
  colnames(result_bowtie) = c('target', 'mismatches')
  result_bowtie$mismatches[result_bowtie$mismatches == ''] = 0
  result_bowtie$mismatches[which(
    result_bowtie$mismatches != 0
  )] = stri_count_regex(
    result_bowtie$mismatches[which(result_bowtie$mismatches != 0)],
    '>'
  )
  result_bowtie$mismatches = as.numeric(result_bowtie$mismatches)
  stats[['sgRNA_hits_0MM']] = length(which(result_bowtie$mismatches == 0))
  stats[['sgRNA_hits_1MM']] = length(which(result_bowtie$mismatches == 1))
  stats[['sgRNA_hits_2MM']] = length(which(result_bowtie$mismatches == 2))
  stats[['sgRNA_hits_3MM']] = length(which(result_bowtie$mismatches == 3))

  result = data.frame(table(result_bowtie$target), stringsAsFactors = FALSE)
  colnames(result) = c('ID', 'Count')
  result$Count = as.numeric(result$Count)
  dict = merge(dict, result, by.x = 'ID', by.y = 'ID', all.x = TRUE)
  dict[is.na(dict$Count), 'Count'] = 0

  ###Export Tables
  resultFile = paste0(sampleName, '-result.xlsx')
  writexl::write_xlsx(dict, resultFile)

  countFile_sgRNA = paste0(sampleName, '-sgRNA_counts.txt')
  sgRNA_counts = data.frame(
    Identifier = dict$ID,
    matchCounts = dict$Count,
    stringsAsFactors = FALSE
  )
  ezWrite.table(sgRNA_counts, countFile_sgRNA, row.names = FALSE)

  if (exists('annotationFile', where = param)) {
    countFile_gene = paste0(sampleName, '-gene_counts.xlsx')
    annot = ezRead.table(param$annotationFile, row.names = NULL)[, c(
      'gene_id',
      'gene_name'
    )]
    res = dict[!dict$isControl, ]
    ## TODO: why does the written excel file not include the controls but the plots seem to use it???
    res = res[order(res$GeneSymbol), ]
    countsPerGene = tapply(res$Count, INDEX = res$GeneSymbol, FUN = sum)
    countsPerGene = data.frame(
      ID = names(countsPerGene),
      matchCounts = countsPerGene,
      stringsAsFactors = FALSE
    )
    countsPerGene = merge(annot, countsPerGene, by.x = 'gene_name', by.y = 'ID')
    countsPerGene = countsPerGene[, c(2:3)]
    colnames(countsPerGene)[1] = 'Identifier'
    writexl::write_xlsx(countsPerGene, countFile_gene)
  }

  ## Count-distribution plots are now built in the Quarto report (log-scale,
  ## control-vs-targeting), so no static PNGs are written here.
  sortedCounts = log2(1 + sort(dict$Count[!dict$isControl]))
  meanCounts <- mean(sortedCounts)
  upperCutOff = meanCounts + param$diffToLogMeanThreshold
  lowerCutOff = meanCounts - param$diffToLogMeanThreshold

  dict2 = dict[order(dict$TargetID), ]
  dict2 = dict2[!dict2$isControl, ]
  targets = unique(dict2$TargetID)
  targetView = data.frame(TargetID = targets, stringsAsFactors = FALSE)
  targetView[['GeneSymbol']] = tapply(dict2$GeneSymbol, dict2$TargetID, unique)
  targetView[['Count']] = tapply(
    dict2$Count,
    dict2$TargetID,
    paste,
    collapse = ','
  )
  targetView[['Count_Sum']] = tapply(dict2$Count, dict2$TargetID, sum)
  targetView[['#sgRNAs > 0']] = tapply(dict2$Count > 0, dict2$TargetID, sum)
  targetView[['#sgRNAs > lowerCutOff']] = tapply(
    dict2$Count > 2^lowerCutOff,
    dict2$TargetID,
    sum
  )
  targetView = targetView[order(targetView[['#sgRNAs > lowerCutOff']]), ]

  underrepTargets = targetView[targetView[['#sgRNAs > lowerCutOff']] < 2, ]
  if (nrow(underrepTargets) > 0 & nrow(underrepTargets) < 1000) {
    underrepTargets = targetView[targetView[['#sgRNAs > lowerCutOff']] < 2, ]
    myDT = DT::datatable(
      underrepTargets,
      escape = F,
      rownames = FALSE,
      filter = 'bottom',
      caption = paste(sampleName, '- UnderrepresentedTargets', sep = ''),
      extensions = c('Buttons'),
      options = list(
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#0000A0', 'color': '#fff'});",
          "}"
        ),
        dom = c('Bfrtip'),
        buttons = c('colvis', 'copy', 'csv', 'excel', 'pdf', 'print'),
        pageLength = 100,
        autoWidth = TRUE
      )
    )
    DT::saveWidget(myDT, 'underrepresentedTargets.html')
  } else {
    myMessage = 'Too many underrepresented target for HTML output. Please check txt-files.'
    write.table(
      myMessage,
      'underrepresentedTargets.html',
      col.names = FALSE,
      row.names = FALSE
    )
  }
  writexl::write_xlsx(targetView, paste0(sampleName, '-targetBasedResult.xlsx'))

  makeQuartoReport(
    param = param,
    output = output,
    dict = dict,
    stats = stats,
    targetView = targetView,
    patternInfo = patternInfo,
    pwm = patternGuess$pwm,
    htmlFile = "00index.html",
    qmdFile = "CountSpacer.qmd",
    reportTitle = paste0("CountSpacer: ", sampleName),
    number = TRUE,
    buttons = TRUE,
    colour = TRUE
  )
  ezWrite.table(
    unlist(stats),
    paste0(sampleName, '-stats.txt'),
    row.names = TRUE
  )
  ezSystem('rm *.fastq.gz')
  ezSystem('pigz --best *.fa')
  return("Success")
}


selectFirst <- function(x) {
  ## if null we return NA
  if (is.null(x)) {
    as.integer(NA)
  } else {
    x[1]
  }
}


##' Infer the constant flanking patterns of a CRISPR read from per-position base
##' composition. The spacer region is variable (max base frequency ~0.25) while
##' the constant leading base(s) and the downstream scaffold are near-invariant.
##' Returns the read-composition PWM (for a seqLogo) plus the guessed left/right
##' flanking sequences ("" when none is evident).
guessFlankingPatterns <- function(
  readFile,
  nSample = 1e5,
  minCoverageFrac = 0.5,
  constFreq = 0.9,
  varFreq = 0.5,
  minSpacerRun = 10L,
  maxFlank = 12L
) {
  require(ShortRead)
  require(Biostrings)
  require(seqLogo)
  emptyGuess <- list(pwm = NULL, leftGuess = "", rightGuess = "", consensus = "")
  strm <- FastqStreamer(readFile, n = nSample)
  on.exit(close(strm))
  fq <- yield(strm)
  reads <- sread(fq)
  if (length(reads) == 0) {
    return(emptyGuess)
  }
  cm <- consensusMatrix(reads, baseOnly = TRUE)
  cm <- cm[c("A", "C", "G", "T"), , drop = FALSE]
  coverage <- colSums(cm)
  keep <- coverage >= (minCoverageFrac * length(reads))
  if (!any(keep)) {
    return(emptyGuess)
  }
  lastPos <- max(which(keep))
  cm <- cm[, seq_len(lastPos), drop = FALSE]
  coverage <- coverage[seq_len(lastPos)]
  probs <- sweep(cm, 2, pmax(coverage, 1), "/")
  probs[, coverage == 0] <- 0.25
  maxFreq <- apply(probs, 2, max)
  consBase <- rownames(cm)[apply(probs, 2, which.max)]
  pwm <- seqLogo::makePWM(probs)
  isConst <- maxFreq >= constFreq
  isVar <- maxFreq < varFreq

  ## Leading constant run -> left flanking pattern (e.g. the U6 'G').
  leftGuess <- ""
  if (isConst[1]) {
    run <- 1L
    while (run < length(isConst) && isConst[run + 1L]) {
      run <- run + 1L
    }
    leftGuess <- paste(consBase[seq_len(min(run, maxFlank))], collapse = "")
  }

  ## First long variable run = spacer; the first constant run that starts after
  ## it = scaffold. A transition base can sit between them (e.g. a 60%-conserved
  ## position), so scan for the next constant RUN rather than the single position
  ## immediately after the spacer.
  rightGuess <- ""
  rv <- rle(isVar)
  vEnds <- cumsum(rv$lengths)
  spacerRuns <- which(rv$values & rv$lengths >= minSpacerRun)
  if (length(spacerRuns) > 0) {
    spacerEnd <- vEnds[spacerRuns[1]]
    rc <- rle(isConst)
    cEnds <- cumsum(rc$lengths)
    cStarts <- cEnds - rc$lengths + 1L
    constRuns <- which(rc$values & cStarts > spacerEnd)
    if (length(constRuns) > 0) {
      s <- cStarts[constRuns[1]]
      e <- min(cEnds[constRuns[1]], s + maxFlank - 1L)
      rightGuess <- paste(consBase[s:e], collapse = "")
    }
  }

  list(
    pwm = pwm,
    leftGuess = leftGuess,
    rightGuess = rightGuess,
    consensus = paste(consBase, collapse = "")
  )
}


twoPatternReadFilter <- function(
  readFile,
  leftPattern,
  rightPattern,
  maxMismatch,
  spacerLength = 0L
) {
  allReads = DNAStringSet()
  processedReads = 0
  dataChunks = 5 * 10^6
  strm <- FastqStreamer(readFile, n = 5 * 10^6)
  repeat {
    currentReads <- yield(strm)
    if (length(currentReads) == 0) {
      break
    }
    reads <- sread(currentReads)
    if (leftPattern != '') {
      vp <- vmatchPattern(leftPattern, reads, max.mismatch = maxMismatch)
      leftEnd <- vp %>% endIndex() %>% vapply(selectFirst, integer(1))
    } else {
      leftEnd <- rep(0L, length(reads))
    }

    if (rightPattern != '') {
      vp <- vmatchPattern(rightPattern, reads, max.mismatch = maxMismatch)
      rightStart <- vp %>% startIndex() %>% vapply(selectFirst, integer(1))
    } else {
      rightStart <- width(reads)
    }

    if (spacerLength > 0 && rightPattern != '') {
      ## Right-anchored, fixed-length spacer: the spacerLength bases immediately
      ## 5' of the right pattern. This drops the constant leading base (U6 'G')
      ## and any staggered filler that a start-of-read extraction keeps -- the
      ## cause of over-long spacers that fail to align.
      spStart <- rightStart - spacerLength
      spEnd <- rightStart - 1L
      ok <- !is.na(rightStart) & spStart >= (leftEnd + 1L) & spStart >= 1L
      reads <- DNAStringSet(substr(reads[ok], spStart[ok], spEnd[ok]))
    } else if (spacerLength > 0 && leftPattern != '') {
      ## Left-anchored fixed-length spacer (no right pattern available).
      spStart <- leftEnd + 1L
      spEnd <- leftEnd + spacerLength
      ok <- !is.na(leftEnd) & spEnd <= width(reads)
      reads <- DNAStringSet(substr(reads[ok], spStart[ok], spEnd[ok]))
    } else {
      ## Legacy behaviour: everything between the two patterns.
      toNA <- which(rightStart < leftEnd)
      rightStart[toNA] <- NA
      patternPositions <- cbind(leftEnd = leftEnd, rightStart = rightStart)
      patternInRead <- !apply(is.na(patternPositions), 1, any)
      patternPositions <- as.data.frame(patternPositions[patternInRead, , drop = FALSE])
      if (rightPattern != '' || leftPattern != '') {
        reads <- reads[patternInRead]
        reads <- DNAStringSet(substr(
          reads,
          patternPositions$leftEnd + 1,
          patternPositions$rightStart - 1
        ))
      }
    }
    processedReads = processedReads + length(currentReads)
    allReads <- c(allReads, reads)
    print(paste0(processedReads / 10^6, 'M reads processed \n'))
  }
  return(allReads)
}


##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodCountSpacer(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppCountSpacer <-
  setRefClass(
    "EzAppCountSpacer",
    contains = "EzApp",
    methods = list(
      ## fastp/ShortRead/Biostrings/bowtie all unconditional. Note: classic Bowtie
      ## (v1), not Bowtie2 -- a different tool with its own citation.
      citation = function() {
        c(
          "Chen, S., Zhou, Y., Chen, Y. & Gu, J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34(17), i884-i890 (2018). https://doi.org/10.1093/bioinformatics/bty560",
          "Morgan, M., Anders, S., Lawrence, M., Aboyoun, P., Pagès, H. & Gentleman, R. ShortRead: a bioconductor package for input, quality assessment and exploration of high-throughput sequence data. Bioinformatics 25(19), 2607-2608 (2009). https://doi.org/10.1093/bioinformatics/btp450",
          "Pagès, H., Aboyoun, P., Gentleman, R. & DebRoy, S. Biostrings: Efficient manipulation of biological strings. R package version 2.80.1. https://doi.org/10.18129/B9.bioc.Biostrings",
          "Langmead, B., Trapnell, C., Pop, M. & Salzberg, S.L. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology 10, R25 (2009). https://doi.org/10.1186/gb-2009-10-3-r25"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodCountSpacer
        name <<- "EzAppCountSpacer"
        appDefaults <<- rbind(
          minReadLength = ezFrame(
            Type = "integer",
            DefaultValue = 18,
            Description = "minimum length of sgRNA"
          ),
          diffToLogMeanThreshold = ezFrame(
            Type = "numeric",
            DefaultValue = 2,
            Description = "log 2 difference relative to the mean log2 counts above/below which counts are called significant"
          ),
          guessPatterns = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "infer the left/right flanking patterns from read base composition when they are not supplied"
          ),
          spacerLength = ezFrame(
            Type = "integer",
            DefaultValue = 0,
            Description = "spacer length to extract, anchored on the flanking pattern; 0 = derive from the library (recommended)"
          )
        )
      }
    )
  )
