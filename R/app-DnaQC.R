###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodDnaQC <- function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  setwdNew(basename(output$getColumn("Report")))
  param$projectId <- sub("\\/.*", "", input$getColumn("BAM")[1]) ## project id is needed for the IGV link
  result <- computeDnaBamStats(input, htmlFile, param)
  return(result)
}

##' @title Computes the DnaQC statistics
##' @description Computes the DnaQC statistics.
##' @template input-template
##' @template htmlFile-template
##' @param resultList a list of results.
##' @template roxygen-template
computeDnaBamStats <- function(input, htmlFile, param, resultList=NULL){
  #TODO: plots per avg coverage from qualimap stats, picard gc bias plot, samstats reports, qualimap as module
  require("GenomicAlignments")
  require("S4Vectors")
  require("ATACseqQC")
  require("XML")
  environment(ezLibComplexity) <- asNamespace('ATACseqQC')
  samples <- input$getNames()
  files <- input$getFullPaths("BAM")
  dataset <- input$meta
  
  resultList = list()
  for (sm in samples){
    message(sm)
    resultList[[sm]] <- list()
    resultList[[sm]]$multiMatchInFileTable <- getBamMultiMatching(param, files[sm], nReads = dataset[sm, "Read Count"])
    if(param$paired){
      png(paste0('fragmentSize_',sm,'.png'),width = 600, height = 500, res = 90)
        fragSizeDist(files[sm], sm)
      dev.off()
    
      png(paste0('libComplexity_',sm,'.png'),width = 600, height = 500, res = 90)
      ezLibComplexity(readsDupFreq(files[sm]), main = paste(sm, 'Library Complexity'))
      dev.off()
    }
    
    #run Qualimap SingleSample
    cmd = paste('unset DISPLAY; qualimap bamqc -bam', files[sm], '-c -nt', param[['cores']], paste0('--java-mem-size=', param[['ram']],'G'), '-outdir', sm)
    ezSystem(cmd)
    
    ##Extract Duplication Rate, InDel Rate:
    qmFile = file.path(sm, 'genome_results.txt')
    #qmFile = '/srv/gstore/projects/p1001/DNAQC_31169_2018-11-06--11-58-53/DNA_QC_Statistics/OBV_35/genome_results.txt'
    all_data = readLines(qmFile)
    resultList[[sm]]$dupRate = as.numeric(sub('\\%', '', sub('^.*duplication rate = ', '', all_data[grep('duplication rate', all_data)])))
    if(!param$paired){
        numReads = as.numeric(gsub(',','',sub('\\%', '', sub('^.*number of reads = ', '', all_data[grep('number of reads', all_data)]))))
    } else {
        numReads = as.numeric(gsub(',','',sub('\\%', '', sub('^.*number of reads = ', '', all_data[grep('number of reads', all_data)]))))/2 
    }
    if(length(resultList[[sm]]$dupRate) == 0){
        dupReads = as.numeric(gsub(',','',sub('\\%', '', sub('^.*number of duplicated reads \\(flagged\\) = ', '', all_data[grep('number of duplicated reads', all_data)]))))
        if(param$paired){
            dupReads = dupReads/2
        }
        resultList[[sm]]$dupRate = dupReads/numReads
        resultList[[sm]]$dupRate = resultList[[sm]]$dupRate*100
    }
    resultList[[sm]]$errorRate = as.numeric(gsub('^.*general error rate = ', '', all_data[grep('general error rate', all_data)]))
    resultList[[sm]]$insertRate = as.numeric(sub('\\%', '', sub('mapped reads with insertion percentage = ', '', all_data[grep('mapped reads with insertion percentage', all_data)])))
    resultList[[sm]]$delRate = as.numeric(sub('\\%', '',sub('mapped reads with deletion percentage = ', '', all_data[grep('mapped reads with deletion percentage', all_data)])))
    resultList[[sm]]$avgCoverage = as.numeric(sub('X', '', sub('mean coverageData = ', '', all_data[grep('mean coverageData', all_data)])))
    resultList[[sm]]$mappingRate = 100* numReads/dataset[sm, "Read Count"]
    ###TODO: add mappingQuality, GC content, insert size
    #####add duplicate rate plot to lib complexity (calc. optical duplicates with picard)
    metricFn = paste0(sm,'_picardDupReport.txt')
    cmd <- paste(preparePicard(), "MarkDuplicates",
                 paste0("I=", files[sm]),
                 paste0("O=", "toDelete.bam"),
                 paste0("M=", metricFn),
                 paste0("REMOVE_DUPLICATES=mark"),
                 paste0("OPTICAL_DUPLICATE_PIXEL_DISTANCE=",param$pixelDist),
                 "> /dev/null")
    ezSystem(cmd)
    ezSystem('rm toDelete.bam')
    metricFn = 'OBV_35_picardDupReport.txt'
    duplicateStats = read.table(metricFn, skip = 6, nrows = 1, header = TRUE)
    resultList[[sm]]$allDuplicates = duplicateStats$READ_PAIR_DUPLICATES
    resultList[[sm]]$optDuplicates = duplicateStats$READ_PAIR_OPTICAL_DUPLICATES
    resultList[[sm]]$readPairs = duplicateStats$READ_PAIRS_EXAMINED
  }
  #run Qualimap Multisample
  ezWrite.table(data.frame(SampleName = samples, Folder = samples), 'sampleKey.txt', row.names = FALSE, col.names = FALSE)
  cmd = paste('unset DISPLAY; qualimap multi-bamqc -d sampleKey.txt -outdir qualimap_MultiSample')
  ezSystem(cmd)
  
  qmFile = 'qualimap_MultiSample/multisampleBamQcReport.html'
  res <- xmlParse(qmFile, isHTML = TRUE)
  sampleSummary <- readHTMLTable(res, which = 3, as.data.frame = FALSE)
  for (j in 2:length(sampleSummary)){
    png(paste0('QualiMapStats_', names(sampleSummary)[j], '.png'))
    barplot(as.numeric(sampleSummary[[j]]), names.arg = sampleSummary[[1]], 
            col = 'royalblue', main = names(sampleSummary)[j])
    dev.off()
  }
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "DnaQC.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="DnaQC.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
}

ezLibComplexity <- function(histFile, times = 100, interpolate.sample.sizes = seq(0.1, 1, by = 0.1), 
                            extrapolate.sample.sizes = seq(5, 20, by = 5), main = c()) {
  total <- histFile[, 1] %*% histFile[, 2]
  suppressWarnings({
    result = ds.rSAC.bootstrap(histFile, r = 1, times = times)
  })
  sequences <- c(interpolate.sample.sizes, extrapolate.sample.sizes)
  estimates <- data.frame(relative.size = sequences, values = rep(NA, 
                                                                  length(sequences)))
  for (i in seq_along(sequences)) {
    estimates$values[i] <- result$f(sequences[i])
  }
  suppressWarnings(estimates$reads <- estimates$relative.size * 
                     total)
  plot(x = estimates$reads/10^6, y = estimates$values/10^6, 
       type = "o", xlab = expression(Putative ~ sequenced ~ 
                                       fragments ~ x ~ 10^6), ylab = expression(Distinct ~ 
                                                                                  fragments ~ x ~ 10^6), main = main)
  return(invisible(estimates))
}

##' @template app-template
##' @templateVar method ezMethodDnaQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
EzAppDnaQC <-
  setRefClass("EzAppDnaQC",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodDnaQC
                  name <<- "EzAppDnaQC"
                  appDefaults <<- rbind(fragSizeMax=ezFrame(Type="integer",  DefaultValue=800,	Description="maximum fragment size to plot in fragment size distribution"),
                                        writeIgvSessionLink=ezFrame(Type="logical", DefaultValue="TRUE", Description="should an IGV link be generated"))
                }
              )
)

