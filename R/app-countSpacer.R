###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCountSpacer = function(input=NA, output=NA, param=NA){
  require(ShortRead)
  require(htmlwidgets)
  require(stringi)
  
  sampleName = input$getNames()
  setwdNew(sampleName)
  csvFiles = list.files(file.path('/srv/GT/databases/GEML/sgRNA_Libs/',param[['dictPath']]), pattern = 'final.csv$', full.names = TRUE)
  if(length(csvFiles) < 1){
      csvFiles <- list.files(file.path('/srv/GT/databases/GEML/sgRNA_Libs/',param[['dictPath']]), pattern = '.csv$', full.names = TRUE)
      param[['dictPath']] = csvFiles[grep('MAGeCK', csvFiles, invert = TRUE)]
  } else {
      param[['dictPath']] = csvFiles[1]
  }
  
  dict = ezRead.table(param[['dictPath']], header = FALSE, sep = ',', row.names = NULL)
  colnames(dict) = c('TargetID', 'Sequence', 'GeneSymbol', 'isControl')
  dict[['ID']] = paste(dict$TargetID, dict$Sequence, sep = '-')
  stats = list()
  stats[['rawReads']] = as.numeric(input$meta[['Read Count']])
  
  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  readFile = trimmedInput$getColumn("Read1")
  stats[['filteredReads']] = as.numeric(ezSystem(paste('zcat', readFile, '|wc -l'), intern = TRUE))/4
  reads <- twoPatternReadFilter(readFile, param$leftPattern, param$rightPattern, param$maxMismatch)
  
  ###Export as fasta file
  readFile = paste0(sampleName,'.fa')
  reads = reads[width(reads) >= param$minReadLength]
  stats[['validSpacerReads']] = length(reads)
  writeXStringSet(reads, readFile, format = 'fasta')
  remove(reads)
  gc()
  
  resultFile = paste0(sampleName, '_bowtie.txt')
  bowtieIndex <- sub('\\.[0-9]\\.ebwt', '', list.files(dirname(param[['dictPath']]), pattern = 'ebwt$', full.names = TRUE)[1])
  cmd = paste('bowtie', bowtieIndex, readFile, '-f -p', param$cores, '|cut -f3,8|sort >', resultFile)
  ezSystem(cmd)
  if(file.size(resultFile)==0){
    stop('No read were aligned to your reference. Please check the correctness of your flanking patterns and the choice of the Crispr reference.')}
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
  
  ###Export Tables
  resultFile = paste0(sampleName,'-result.xlsx')
  writexl::write_xlsx(dict, resultFile)
  
  countFile_sgRNA = paste0(sampleName,'-sgRNA_counts.txt')
  sgRNA_counts = data.frame(Identifier = dict$ID, matchCounts = dict$Count, stringsAsFactors = FALSE)
  ezWrite.table(sgRNA_counts, countFile_sgRNA, row.names=FALSE)
  
  if(exists('annotationFile', where = param)){
    countFile_gene = paste0(sampleName,'-gene_counts.xlsx')
    annot = ezRead.table(param$annotationFile, row.names = NULL)[ ,c('gene_id', 'gene_name')]
    res = dict[!dict$isControl,]
    ## TODO: why does the written excel file not include the controls but the plots seem to use it???
    res = res[order(res$GeneSymbol), ]
    countsPerGene = tapply(res$Count, INDEX = res$GeneSymbol, FUN = sum)
    countsPerGene = data.frame(ID = names(countsPerGene), matchCounts = countsPerGene, stringsAsFactors = FALSE)
    countsPerGene = merge(annot, countsPerGene, by.x = 'gene_name', by.y = 'ID')
    countsPerGene = countsPerGene[,c(2:3)]
    colnames(countsPerGene)[1] = 'Identifier'
    writexl::write_xlsx(countsPerGene, countFile_gene)
  }

  
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
  
  ## TODO: the equivalent code is also in the RMD but there the controls are not included
  sortedCounts = log2(1+sort(dict$Count))
  meanCounts <- mean(sortedCounts)
  upperCutOff = meanCounts + param$diffToLogMeanThreshold
  lowerCutOff = meanCounts - param$diffToLogMeanThreshold
  

  dict2 = dict[order(dict$TargetID), ]
  dict2 = dict2[!dict2$isControl,]
  targets = unique(dict2$TargetID)
  targetView = data.frame(TargetID = targets, stringsAsFactors = FALSE)
  targetView[['GeneSymbol']] = tapply(dict2$GeneSymbol, dict2$TargetID, unique)
  targetView[['Count']] = tapply(dict2$Count, dict2$TargetID, paste, collapse = ',')
  targetView[['Count_Sum']] = tapply(dict2$Count, dict2$TargetID, sum)
  targetView[['#sgRNAs > 0']] = tapply(dict2$Count > 0, dict2$TargetID, sum)
  targetView[['#sgRNAs > lowerCutOff']] = tapply(dict2$Count > 2^lowerCutOff, dict2$TargetID, sum)
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
  } else {
      myMessage ='Too many underrepresented target for HTML output. Please check txt-files.'
      write.table(myMessage, 'underrepresentedTargets.html', col.names = FALSE, row.names = FALSE)
  }
  writexl::write_xlsx(targetView, paste0(sampleName,'-targetBasedResult.xlsx'))
  
  makeRmdReport(param=param, output=output, dict=dict, stats=stats, htmlFile = "00index.html", rmdFile = "CountSpacer.Rmd", reportTitle = paste0("CountSpacer: ", sampleName))  
  ezWrite.table(unlist(stats), paste0(sampleName,'-stats.txt'), row.names = TRUE)
  ezSystem('rm *.fastq.gz')
  ezSystem('pigz --best *.fa')
  return("Success")
}


selectFirst <- function(x) {
  ## if null we return NA
  if (is.null(x))
    as.integer(NA)
  else
    x[1]
}


twoPatternReadFilter <- function(readFile, leftPattern, rightPattern, maxMismatch) {
  allReads = DNAStringSet()
  processedReads = 0
  dataChunks = 5*10^6
  strm <- FastqStreamer(readFile, n = 5*10^6)
  repeat {
    currentReads <- yield(strm)
    if (length(currentReads) == 0)
      break
    reads <- sread(currentReads)
    if(leftPattern != ''){
      vp <- vmatchPattern(leftPattern, reads, max.mismatch = maxMismatch)
      leftEnd <- vp %>% endIndex() %>% vapply(selectFirst, integer(1))
    } else {
      leftEnd <- rep(0, length(reads))
    }
    
    if(rightPattern != ''){
        vp <- vmatchPattern(rightPattern, reads, max.mismatch = maxMismatch)
        rightStart <- vp %>% startIndex() %>% vapply(selectFirst, integer(1))
    } else {
        rightStart <- width(reads)
    }
    
    toNA <- which(rightStart < leftEnd)
    rightStart[toNA] <- NA
    patternPositions <- cbind(leftEnd = leftEnd,
                              rightStart = rightStart)
    patternInRead <- !apply(is.na(patternPositions), 1, any)
    patternPositions <- as.data.frame(patternPositions[patternInRead, ])
    if(rightPattern != '' || leftPattern != ''){
        reads <- reads[patternInRead]
        reads <- DNAStringSet(substr(reads, patternPositions$leftEnd+1, patternPositions$rightStart-1))
    }
    processedReads = processedReads + length(currentReads)
    allReads <- c(allReads, reads)
    print(paste0(processedReads/10^6, 'M reads processed \n'))
  }
  return(allReads)
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
                  appDefaults <<- rbind(
                    minReadLength = ezFrame(Type = "integer", DefaultValue = 18,
                                          Description = "minimum length of sgRNA"),
                    diffToLogMeanThreshold = ezFrame(
                        Type = "numeric",
                        DefaultValue = 2,
                        Description = "log 2 difference relative to the mean log2 counts above/below which counts are called significant"
                    )
                    )
                    }
              )
  )
