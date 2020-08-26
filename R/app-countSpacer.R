###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCountSpacer = function(input=NA, output=NA, param=NA){
  require(Biostrings)
  require(ShortRead)
  require(ggplot2)
  require(htmlwidgets)
  require(stringi)
  
  setwdNew(param[['name']])
  sampleName = input$getNames()
  param[['dictPath']] = list.files(file.path('/srv/GT/databases/GEML/sgRNA_Libs/',param[['dictPath']]), pattern = 'csv', full.names = TRUE)
  dict = ezRead.table(param[['dictPath']], header = FALSE, sep = ',', row.names = NULL)
  colnames(dict) = c('TargetID', 'Sequence', 'GeneSymbol', 'isControl')
  dict[['ID']] = paste(dict$TargetID, dict$Sequence, sep = '-')
  stats = list()
  stats[['rawReads']] = as.numeric(input$meta[['Read Count']])
  
  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  readFile = trimmedInput$getColumn("Read1")
  stats[['filteredReads']] = as.numeric(ezSystem(paste('cat', readFile, '|wc -l'), intern = TRUE))/4
  
  reads <- .getReadsFromFastq(readFile)
  reads <- twoPatternReadFilter(reads, param$leftPattern, param$rightPattern, param$maxMismatch)
  
  ###Export as fasta file
  readFile = paste0(sampleName,'.fa')
  reads = reads[width(reads) >= param$minReadLength]
  stats[['validSpacerReads']] = length(reads)
  writeXStringSet(reads, readFile, format = 'fasta')
  
  ###Run Bowtie1: ###TODO: runBowtie1 only 1 time (combined BAM) -> filter reads for MM
  resultFile = paste0(sampleName, '_bowtie.txt')
  cmd = paste('bowtie', sub('.csv', '', param[['dictPath']]), readFile, '-f -p', param$cores, '|cut -f3,8|sort >', resultFile)
  ezSystem(cmd)
  result_bowtie = ezRead.table(resultFile, row.names = NULL, header = FALSE)
  colnames(result_bowtie) = c('target', 'mismatches')
  result_bowtie$mismatches[result_bowtie$mismatches == ''] = 0
  result_bowtie$mismatches[which(result_bowtie$mismatches != 0)] = stri_count_regex(result_bowtie$mismatches[which(result_bowtie$mismatches != 0)], '>')
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
  
  ###Export Table
  ezWrite.table(dict, paste0(sampleName,'-result.txt') ,row.names = FALSE)
  data = data.frame(group = rep('Count', nrow(dict)), counts = c(dict$Count))
  # Basic violin plot
  p <- ggplot(data, aes(x=group, y=counts))
  p <- p + geom_violin(fill="royalblue", alpha= 0.5, trim = FALSE, adjust = 0.5) 
  p <- p + geom_boxplot(width = 0.1)
  p <- p +  ggtitle(paste0(sampleName, '-ReadCount Distribution')) + ylab('ReadCount per sgRNA')
  p <- p + theme(plot.title = element_text(size = 15, face = "bold"), axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  png(paste0(sampleName, '_violin.png'), 600, 1000, res = 100)
    print(p)
  dev.off()
  
  h <- ggplot(dict, aes(x=log2(1+Count))) + geom_histogram(binwidth=0.1)
  h <- h +  ggtitle(paste0(sampleName, '-Histogram')) + ylab('Number of sgRNAs') + xlab('Log2 count per sgRNA')
  png(paste0(sampleName, '_histogram.png'), 800, 500, res = 100)
    print(h)
  dev.off()
  
  sortedCounts = log2(1+sort(dict$Count))
  upperCutOff = 2+mean(sortedCounts)
  lowerCutOff = mean(sortedCounts)-2
  
  up_sgRNAs = length(which(sortedCounts>upperCutOff))
  relUp_sgRNA = round(100 * (up_sgRNAs/length(sortedCounts)), digits = 2)
  down_sgRNAs = length(which(sortedCounts<lowerCutOff))
  relDown_sgRNA = round(100 * (down_sgRNAs/length(sortedCounts)), digits = 2)
  mu = round(mean(sortedCounts), 2)
  png(paste0(sampleName, '_overview.png'), 700, 700, res = 100)
  plot(sortedCounts, pch = c(15), cex = 0.7, main = paste(sampleName, '- sgCount Overview'), ylab = 'log2(sgRNA Count)')
  abline(h = mean(sortedCounts))
  abline(h = upperCutOff, lty = 2)
  abline(h = lowerCutOff, lty = 2)
  text(length(sortedCounts)*0.05, 1.05*mean(sortedCounts), bquote(mu==.(mu)), cex = 0.8)
  text(length(sortedCounts)*0.8, 1.05*upperCutOff, paste0('#',up_sgRNAs, ' (',relUp_sgRNA,'%)' ), cex = 0.8)
  text(length(sortedCounts)*0.15, 0.95*lowerCutOff, paste0('#',down_sgRNAs, ' (',relDown_sgRNA,'%)' ), cex = 0.8)
  dev.off()
  
  dict2 = dict[order(dict$TargetID), ]
  dict2 = dict2[!dict2$isControl,]
  targets = unique(dict2$TargetID)
  targetView = data.frame(TargetID = targets, stringsAsFactors = FALSE)
  targetView[['GeneSymbol']] = tapply(dict2$GeneSymbol, dict2$TargetID, unique)
  targetView[['Count']] = tapply(dict2$Count, dict2$TargetID, paste, collapse = ',')
  targetView[['Count_Sum']] = tapply(dict2$Count, dict2$TargetID, sum)
  targetView[['#sgRNAs > 0']] = tapply(dict2$Count, dict2$TargetID, .greaterMin, 0)
  targetView[['#sgRNAs > lowerCutOff']] = tapply(dict2$Count, dict2$TargetID, .greaterMin, 2^lowerCutOff)
  targetView = targetView[order(targetView[['#sgRNAs > lowerCutOff']]),]
  
  underrepTargets = targetView[targetView[['#sgRNAs > lowerCutOff']] < 2,]
  if(nrow(underrepTargets) > 0 & nrow(underrepTargets) < 1000){
    underrepTargets = targetView[targetView[['#sgRNAs > lowerCutOff']] < 2, ]
    myDT = DT::datatable(underrepTargets, escape = F, rownames = FALSE, filter = 'bottom',
                         caption = paste(sampleName, '- UnderrepresentedTargets',sep=''),extensions = c('Buttons'),
                         options = list(initComplete = JS(
                           "function(settings, json) {",
                           "$(this.api().table().header()).css({'background-color': '#0000A0', 'color': '#fff'});",
                           "}"),
                           dom = c('Bfrtip'),buttons = c('colvis','copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 100, autoWidth=TRUE))
    saveWidget(myDT, 'underrepresentedTargets.html')
  }
  ezWrite.table(targetView, paste0(sampleName,'-targetBasedResult.txt') ,row.names = FALSE)
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "CountSpacer.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  htmlFile = '00index.html'
  rmarkdown::render(input="CountSpacer.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  ezWrite.table(unlist(stats), paste0(sampleName,'-stats.txt'), row.names = TRUE)
  remove(reads)
  ezSystem('rm *.fastq')
  ezSystem('pigz --best *.fa')
  return("Success")
}

.dummyFunction <- function(x) {
  if (is.null(x))
    NA
  else
    x[1]
}

twoPatternReadFilter <- function(reads, leftPattern, rightPattern, maxMismatch) {
  vp <- vmatchPattern(leftPattern, reads, max.mismatch = maxMismatch)
  leftEnd <- vapply(endIndex(vp),
                    .dummyFunction, c('endIndex' = 0))
  vp <- vmatchPattern(rightPattern, reads, max.mismatch = maxMismatch)
  rightStart <- vapply(startIndex(vp),
                       .dummyFunction, c('startIndex' = 0))
  toNA <- which(rightStart < leftEnd)
  rightStart[toNA] <- NA
  patternPositions <- cbind(leftEnd = leftEnd,
                            rightStart = rightStart)
  patternInRead <- !apply(is.na(patternPositions), 1, any)
  patternPositions <- as.data.frame(patternPositions[patternInRead, ])
  reads <- reads[patternInRead]
  reads <- DNAStringSet(substr(reads, patternPositions$leftEnd+1, patternPositions$rightStart-1))
  return(reads)
}

.getReadsFromFastq <- function(file) {
  f <- FastqFile(file)
  myReads <- sread(readFastq(f, ))
  close(f)
  return(myReads)
}

.greaterMin <- function(n, minVal){
  return(length(which(n > minVal)))
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodCountSpacer(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppCountSpacer <-
  setRefClass("EzAppCountSpacer",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCountSpacer
                  name <<- "EzAppCountSpacer"
                }
              )
  )
