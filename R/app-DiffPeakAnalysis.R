###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodDiffPeakAnalysis <- function(input = NA, output = NA, param = NA){
    grouping <- input$getColumn(param$grouping)
    stopifnot(param$sampleGroup != param$refGroup)
    grouping <- grouping[grouping %in% c(param$refGroup, param$sampleGroup)]
    commonCols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

    countFiles <- input$getFullPathsList("Count")
    featureCounts <- loadCountFiles(countFiles, grouping, commonCols)

    dds <- generateDESeqDS(featureCounts, commonCols, grouping)
    peakAnno <- annotateConsensusPeaks(gtfFile = param$ezRef@refFeatureFile, fastaFile = param$ezRef@refFastaFile, 
                                       peakFile = countFiles[[1]], tool = param$annotationMethod, cores = param$cores)
    outDir <- file.path(basename(output$getColumn('ResultFolder')))
    cd = getwd()
    setwdNew(outDir)
    makeRmdReport(
      output = output, param = param, peakAnno = peakAnno, dds=dds, selfContained = TRUE,
      rmdFile = "DiffPeakAnalysis.Rmd", htmlFile = "00index.html",
      reportTitle = 'Differential Peak Analysis', use.qs2 = TRUE
    )
    setwd(cd)
}


#' @template app-template
##' @templateVar method ezMethodDiffPeakAnalysis(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run a differential expression analysis with the application edgeR on two groups.
EzAppDiffPeakAnalysis <-
  setRefClass("EzAppDiffPeakAnalysis",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodDiffPeakAnalysis
        name <<- "EzAppDiffPeakAnalysis"
        appDefaults <<- rbind(annotationMethod=ezFrame(Type="character", DefaultValue="chippeakanno", Description="peak annotation method"))  
      }
      )
    )

##' @description generate file with all counts for sampleGroup and refGroup
loadCountFiles <- function(countFiles, grouping, commonCols){
  countFilesSubset <- countFiles[names(grouping)]

  loadAllTables <- imap(countFilesSubset, function(file_i, listName) {
    group <- grouping[[listName]]
    rep_id <- match(listName, names(grouping)[grouping == group])
    new_col <- paste0(group, "_REP", rep_id)
    data.table::fread(file_i, data.table=FALSE) %>%
      rename(!!new_col := matchCounts)
  })

  purrr::reduce(loadAllTables, full_join, by = commonCols)
}

##' @description generate DESeqDataSet from the counts table
generateDESeqDS <- function(featureCounts, commonCols, grouping){
  library(DESeq2)
  countCols <- setdiff(colnames(featureCounts), commonCols)

  countData <- featureCounts %>%
    select(all_of(countCols)) %>%
    as.data.frame()
  countData <- round(countData)

  colData <- ezFrame(
    sample = countCols,
    group  = sub("_REP.*", "", countCols),
    replicate = sub(".*_REP", "", countCols),
    row.names = countCols
  )

  dds <- DESeqDataSetFromMatrix(
    countData = countData,
    colData   = colData,
    design    = ~ group
  )
  rowData(dds) <- featureCounts[, commonCols]
  rownames(dds) <- featureCounts$Geneid
  dds$Condition <- dds$group
  dds
}


##' @description generate peaks annotation file using the method specified with tool
annotateConsensusPeaks <- function(gtfFile, peakFile, fastaFile, tool, cores){
  switch (as.character(tool),
          chippeakanno = {
            require(ChIPpeakAnno)
            require(GenomicRanges)
            require(rtracklayer)
            gtf <- rtracklayer::import(gtfFile)
            if('gene' %in% unique(gtf$type)){
              idx = gtf$type == 'gene'
            } else if('transcript' %in% unique(gtf$type)) {
              idx = gtf$type == 'transcript'
            } else if('start_codon' %in% unique(gtf$type)){
              idx = gtf$type =='start_codon'
            } else {
              message('gtf is incompatabible. Peak annotation skipped!')
              return(NULL)
            }
            gtf = gtf[idx]
            if(grepl('gtf$',gtfFile)){
              names_gtf = make.unique(gtf$'gene_id')
            } else {
              names_gtf = make.unique(gtf$'ID')
            }
            names(gtf) = names_gtf
            myPeaks = ezRead.table(peakFile)
            peaksRD = makeGRangesFromDataFrame(myPeaks, keep.extra.columns = TRUE, start.field = "Start", end.field = "End", seqnames.field="Chr")
            annotChIPpeak <- annotatePeakInBatch(peaksRD, AnnotationData = gtf, output='nearestStart', multiple=FALSE, FeatureLocForDistance='TSS')
            annotChIPpeak <- as.data.frame(annotChIPpeak)
            annotChIPpeak <- annotChIPpeak %>% rename("feature_start" = "start_position" , "feature_end" = "end_position")
            annotChIPpeak <- annotChIPpeak %>% rename("peakId" = "peak", "geneName" = "feature")
            return(annotChIPpeak)
          },
          chipseeker = {
            require(ChIPseeker)
            require(GenomicRanges)
            require(rtracklayer)
            myTxDB <- txdbmaker::makeTxDbFromGFF(file=gtfFile, format='gtf')
            myPeaks = ezRead.table(peakFile)
            myPeaks$peakId = rownames(myPeaks)
            peaksRD = makeGRangesFromDataFrame(myPeaks, keep.extra.columns = TRUE, start.field = "Start", end.field = "End", seqnames.field="Chr")
            annotChIPseeker <- annotatePeak(peaksRD, TxDb=myTxDB, tssRegion=c(-1000, 1000), verbose=FALSE)
            annotChIPseeker <- data.frame(annotChIPseeker@anno)
            keepColChIPSeeker <- c("seqnames", "peakId", "annotation", "geneId", "transcriptId", "distanceToTSS", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand")
            annotChIPseeker <- annotChIPseeker[,keepColChIPSeeker]
            annotChIPseeker <- annotChIPseeker %>% rename("geneName" = "geneId")
            return(annotChIPseeker)
          },
          homer = {
            myPeaks <- ezRead.table(peakFile)
            bedFileCols <- c("Chr", "Start", "End")
            bedFile <- myPeaks[,bedFileCols]
            bedFile$Names <- rownames(myPeaks)
            bedFile$Score <- 0
            bedFile$Strand <- myPeaks[["Strand"]]
            bedFileName <- "peaks.bed"
            write.table(bedFile, bedFileName, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
            cmd = paste(
              "annotatePeaks.pl",
              bedFileName,
              fastaFile,
              "-gtf", gtfFile,
              "-cpu", cores,
              "> annotatedPeaks.txt"
            )
            system(cmd)
            annotHomer <- ezRead.table('annotatedPeaks.txt', row.names = NULL)
            colnames(annotHomer)[1] <- 'peakId'
            annotHomer <- annotHomer %>% rename("geneName" = "Gene Name")
            return(annotHomer)
          }
  )
}
