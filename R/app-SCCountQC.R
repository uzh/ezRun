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
                  #appDefaults <<- rbind(runGO=ezFrame(Type="logical", DefaultValue=TRUE, Description="whether to run the GO analysis"))
                }
              )
  )

ezMethodSCCountQC = function(input=NA, output=NA, param=NA, 
                             htmlFile="00index.html"){
  if(length(input$getNames()) > 1L)
    stop("Currently we support one pooled bam file!")
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd), add=TRUE)
  
  sce <- loadSCCountDataset(input, param)
  
  inBam <- getBamLocally(input$getFullPaths("BAM"))
  on.exit(file.remove(inBam), add = TRUE)
  
  ## STAR log
  mlog <- read.table(input$getFullPaths("STARLog"), sep="|", 
                     as.is = TRUE, quote = "\"", fill=T)
  rownames(mlog) <- trimws(mlog[,1])
  metadata(sce)$mlog <- mlog
  

  ## 3' and 5' bias
  minCount <- 20
  minTxLength <- 1e3
  maxTxs <- 2e3
  if(metadata(sce)$featureLevel == "gene"){
    useGeneIDs <- rownames(sce)[(rowSums(assays(sce)$counts > minCount) >= 1)]
    useTxIDs <- strsplit(rowData(sce)$transcript_id[rownames(sce) %in% useGeneIDs], "; ")
    useTxIDs <- unlist(useTxIDs)
  }else{
    useTxIDs <- rownames(sce)[(rowSums(assays(sce)$counts > minCount) >= 1)]
  }
 
  #txEndBias(param, inBam=inBam, minTxLength=minTxLength, useTxIDs=useTxIDs, maxTxs=maxTxs)
  
  
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "SCCountQC.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  ## debug
  saveRDS(sce, file="sce.rds")
  rmarkdown::render(input="SCCountQC.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  # Picard metrics
  # inBam <- input$getFullPaths("BAM")
  # bamRGFns <- splitBamByRG(inBam, mc.cores=param$cores)
  # on.exit(file.remove(bamRGFns), add = TRUE)
  # 
  # ## CollectAlignmentSummaryMetrics
  # message("Start CollectAlignmentSummaryMetrics", date())
  # alnMetrics <- CollectAlignmentSummaryMetrics(inBams=bamRGFns,
  #                                              fastaFn=param$ezRef['refFastaFile'],
  #                                              metricLevel="SAMPLE",
  #                                              mc.cores=param$cores)
  # save(alnMetrics, file="alnMetrics.rda")
  # message("End CollectAlignmentSummaryMetrics", date())
  # ## CollectRnaSeqMetrics
  # rnaSeqMetrics <- CollectRnaSeqMetrics(inBams=bamRGFns,
  #                                       gtfFn=param$ezRef['refFeatureFile'],
  #                                       featAnnoFn=param$ezRef['refAnnotationFile'],
  #                                       strandMode=param$strandMode,
  #                                       metricLevel="SAMPLE",
  #                                       mc.cores=param$cores)
  # save(rnaSeqMetrics, file="rnaSeqMetrics.rda")
  # message("End CollectRnaSeqMetrics", date())
  # ## DuplicationMetrics
  # dupMetrics <- DuplicationMetrics(inBams=bamRGFns, mc.cores=param$cores)
  # save(dupMetrics, file="dupMetrics.rda")
  # message("End DuplicationMetrics", date())

}

txEndBias <- function(param, inBam, width=100L, minTxLength=NULL,
                      useTxIDs=NULL, maxTxs=NULL){
  require(matrixStats)
  
  ## 5' 100bp and 3' 100bp
  gtf5TempFn <- tempfile(pattern="trimGTF-5-", tmpdir=".", fileext=".gtf")
  gtf3TempFn <- tempfile(pattern="trimGTF-3-", tmpdir=".", fileext=".gtf")
  
  trimTxGtf(param, outGTF=c(gtf5TempFn, gtf3TempFn), 
            width=width, fix=c("start", "end"),
            minTxLength=minTxLength, useTxIDs=useTxIDs, maxTxs=maxTxs)
  counts5 <- txFeatureCounts(param, inBam, gtf5TempFn)
  counts3 <- txFeatureCounts(param, inBam, gtf3TempFn)
  stopifnor(identical(rownames(counts5), rownames(counts3)))
  
  ## All txs
  countsAll <- txFeatureCounts(param, inBam, param$ezRef@refFeatureFile)
  countsAll <- countsAll[rownames(counts5), ]
  
  widthsTx <- getTranscriptGcAndWidth(param)
  widthsTx <- widthsTx$width[rownames(counts5)]
  bias5 <- (counts5 / 100) / (countsAll / widthsTx)
  bias3 <- (counts3 / 100) / (countsAll / widthsTx)
  
  bias5 <- setNames(colMedians(bias5, na.rm=TRUE), colnames(bias5))
  bias3 <- setNames(colMedians(bias3, na.rm=TRUE), colnames(bias3))
  
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