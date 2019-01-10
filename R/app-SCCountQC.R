###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSCCountQC <-
  setRefClass("EzAppSCCountQC",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSCCountQC
                  name <<- "EzAppSCCountQC"
                  appDefaults <<- rbind(minReadsPerCell=ezFrame(Type="numeric", DefaultValue=1500, Description="Filter cells with less reads counted on genes"),
                                        minReadsPerGene=ezFrame(Type="numeric", DefaultValue=3, Description="Minimal number of reads per gene to be expressed"),
                                        minGenesPerCell=ezFrame(Type="numeric", DefaultValue=500, Description="Filter cells with less genes expressed"))
                }
              )
  )

ezMethodSCCountQC = function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  require(S4Vectors)
  require(SummarizedExperiment)
  require(Rsamtools)
  require(GenomeInfoDb)
  
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  reportCwd <- getwd()
  
  ## STAR log: available for smart-seq2, not for 10x
  param$scProtocol <- ifelse("STARLog" %in% input$colNames, "smart-Seq2", "10x")
  
  sce <- loadSCCountDataset(input, param)
  
  ## SCCountQC always process one plate/sample once; remove the plate name
  colnames(sce) <- sub("^.*___", "", colnames(sce))
  
  metadata(sce)$param$name <- paste(metadata(sce)$param$name,
                                    paste(input$getNames(), collapse=", "),
                                    sep=": ")
  param <- metadata(sce)$param
  if(param$scProtocol == "smart-Seq2"){
    mlog <- read.table(input$getFullPaths("STARLog"), sep="|", 
                       as.is = TRUE, quote = "\"", fill=T)
    rownames(mlog) <- trimws(mlog[,1])
    metadata(sce)$mlog <- mlog
    
    inBam <- getBamLocally(input$getFullPaths("BAM"))
    on.exit(file.remove(file.path(reportCwd, c(inBam, paste0(inBam, ".bai")))), 
            add = TRUE)
  }
  if(all(levels(seqnames(rowRanges(sce))) %in%
         seqnames(seqinfo(FaFile(param$ezRef['refFastaFile']))))){
    param$hasControlSeqs <- FALSE
  }else{
    param$hasControlSeqs <- TRUE
  }
  
  ## 3' and 5' bias: not in use currently
  # minCount <- 20
  # minTxLength <- 1e3
  # maxGenes <- 1e3
  # maxTxs <- 1e3
  # ### Choose the most abundant genes/txs
  # if(metadata(sce)$featureLevel == "gene"){
  #   useGeneIDs <- names(head(sort(rowSums(assays(sce)$counts), decreasing=TRUE),
  #                            maxGenes))
  #   #useGeneIDs <- rownames(sce)[(rowSums(assays(sce)$counts > minCount) >= 1)]
  #   useTxIDs <- strsplit(rowData(sce)$transcript_id[rownames(sce) %in% useGeneIDs], "; ")
  #   useTxIDs <- sapply(useTxIDs, "[", 1)
  #   #useTxIDs <- unlist(useTxIDs)
  # }else{
  #   #useTxIDs <- rownames(sce)[(rowSums(assays(sce)$counts > minCount) >= 1)]
  #   useTxIDs <- names(head(sort(rowSums(assays(sce)$counts), decreasing=TRUE),
  #                          maxTxs))
  # }
  #primeBias <- txEndBias(param, inBam=inBam, minTxLength=minTxLength,
  #                       useTxIDs=useTxIDs, maxTxs=maxTxs)
  #colData(sce)$bias5 <- primeBias$bias5[colnames(sce)]
  #colData(sce)$bias3 <- primeBias$bias3[colnames(sce)]
  
  ## Picard metrics
  if(param$scProtocol == "smart-Seq2" && !param$hasControlSeqs){
    bamRGFns <- splitBamByRG(inBam, mc.cores=param$cores)
    on.exit(file.remove(file.path(reportCwd, bamRGFns)), add = TRUE)
    
    ### CollectAlignmentSummaryMetrics
    message("Start CollectAlignmentSummaryMetrics", date())
    alnMetrics <- CollectAlignmentSummaryMetrics(inBams=bamRGFns,
                                                 fastaFn=param$ezRef['refFastaFile'],
                                                 paired=param$paired,
                                                 metricLevel="SAMPLE",
                                                 mc.cores=param$cores)
    colData(sce) <- as(data.frame(colData(sce), alnMetrics[colnames(sce), ],
                                  check.names = FALSE, stringsAsFactors = FALSE),
                       "DataFrame")
    
    ### CollectRnaSeqMetrics
    message("Start CollectRnaSeqMetrics", date())
    rnaSeqMetrics <- CollectRnaSeqMetrics(inBams=bamRGFns,
                                          gtfFn=param$ezRef['refFeatureFile'],
                                          featAnnoFn=param$ezRef['refAnnotationFile'],
                                          strandMode=param$strandMode,
                                          metricLevel="SAMPLE",
                                          mc.cores=param$cores)
    colData(sce) <- as(data.frame(colData(sce), rnaSeqMetrics[colnames(sce), ],
                                  check.names = FALSE, stringsAsFactors = FALSE),
                       "DataFrame")
    
    ### DuplicationMetrics
    message("Start DuplicationMetrics", date())
    dupMetrics <- DuplicationMetrics(inBams=bamRGFns, mc.cores=param$cores)
    colData(sce) <- as(data.frame(colData(sce), dupMetrics[colnames(sce), ],
                                  check.names = FALSE, stringsAsFactors = FALSE),
                       "DataFrame")
  }
  metadata(sce)$param <- param
  saveRDS(sce, file="sce.rds")
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCCountQC.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="SCCountQC.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)

  return("Success")
}

txEndBias <- function(param, inBam, width=100L, minTxLength=NULL,
                      useTxIDs=NULL, maxTxs=NULL){
  require(matrixStats)
  
  ## 5' 100bp and 3' 100bp
  gtf5TempFn <- tempfile(pattern="trimGTF-5-", tmpdir=".", fileext=".gtf")
  gtf3TempFn <- tempfile(pattern="trimGTF-3-", tmpdir=".", fileext=".gtf")
  on.exit(file.remove(c(gtf5TempFn, gtf3TempFn)), add=TRUE)
  
  trimTxGtf(param, outGTF=c(gtf5TempFn, gtf3TempFn), 
            width=width, fix=c("start", "end"),
            minTxLength=minTxLength, useTxIDs=useTxIDs, maxTxs=maxTxs)
  counts5 <- txFeatureCounts(param, inBam, gtf5TempFn)
  counts3 <- txFeatureCounts(param, inBam, gtf3TempFn)
  stopifnot(identical(rownames(counts5), rownames(counts3)))
  
  ## All txs
  countsAll <- txFeatureCounts(param, inBam, param$ezRef@refFeatureFile)
  countsAll <- countsAll[rownames(counts5), ]
  
  widthsTx <- getTranscriptGcAndWidth(param)
  widthsTx <- widthsTx$featWidth[rownames(countsAll)]
  bias5 <- (counts5 / width) / (countsAll / widthsTx)
  bias3 <- (counts3 / width) / (countsAll / widthsTx)
  
  #bias5 <- setNames(colMedians(bias5, na.rm=TRUE), colnames(bias5))
  #bias3 <- setNames(colMedians(bias3, na.rm=TRUE), colnames(bias3))
  
  bias5 <- colMeans(bias5, na.rm=TRUE)
  bias3 <- colMeans(bias3, na.rm=TRUE)
  
  return(list(bias5=bias5, bias3=bias3))
}


fixFeatureCountsRGMatrix <- function(counts, inBam){
  require(Rsamtools)
  countsFixed <- counts
  bamHeaders <- scanBamHeader(inBam)
  
  colnames(countsFixed) <- sub(paste0(make.names(inBam), "."), "", 
                               colnames(countsFixed))
  tagsRG <- sub("ID:", "", 
                sapply(bamHeaders[[1]]$text[names(bamHeaders[[1]]$text) == "@RG"], 
                       "[", 1))
  ## RG starts with numbers will have X after make.names
  ## But featureCounts doesn't have this X.
  fixNameMapping <- setNames(tagsRG, make.names(tagsRG))
  indexStartNumber <- grep("^\\d", fixNameMapping)
  names(fixNameMapping)[indexStartNumber] <- sub("^X", "", 
                                                 names(fixNameMapping)[indexStartNumber])
  colnames(countsFixed) <- fixNameMapping[colnames(countsFixed)]
  return(countsFixed)
}

txFeatureCounts <- function(param, inBam, gtfFn){
  require(Rsubread)
  require(Rsamtools)
  bamHeaders <- scanBamHeader(inBam)
  hasRG <- "@RG" %in% names(bamHeaders[[1]]$text)
  
  countResult <- featureCounts(inBam, annot.inbuilt=NULL,
                               annot.ext=gtfFn, isGTFAnnotationFile=TRUE,
                               GTF.featureType="exon",
                               GTF.attrType="transcript_id",
                               allowMultiOverlap=TRUE,
                               isPairedEnd=param$paired,
                               nthreads=param$cores,
                               strandSpecific=switch(param$strandMode, "both"=0, "sense"=1, "antisense"=2, stop("unsupported strand mode: ", param$strandMode)),
                               minMQS=param$minMapQuality,
                               minOverlap=param$minFeatureOverlap,
                               countMultiMappingReads=param$keepMultiHits,
                               primaryOnly=TRUE,
                               byReadGroup=ifelse(hasRG, TRUE, FALSE))
  ans <- fixFeatureCountsRGMatrix(countResult$counts, inBam)
  return(ans)
}