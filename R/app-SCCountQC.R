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
  useTxIDs <- sample(useTxIDs, size=min(length(useTxIDs), maxTxs))
  
  txEndBias(param, inBam=input$getFullPaths("BAM"), minTxLength=minTxLength, maxTxs=maxTxs, useTxIDs=useTxIDs)
  
  
  ## debug
  #save(sce, file="sce.rdata")

  
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

txEndBias <- function(param, inBam, maxTxs=1e3, width=100L, minTxLength=NULL,
                      useTxIDs=NULL){
  require(Rsubread)
  ## 5' 100bp
  gtfTempFn <- tempfile(pattern="trimGTF-5-", fileext=".gtf")
  trimTxGtf(param, outGTF=gtfTempFn, width=width, fix="start",
            minTxLength=minTxLength, useTxIDs=useTxIDs)
  
  ## 3' 100 bp
  gtfTempFn <- tempfile(pattern="trimGTF-3-", fileext=".gtf")
  trimTxGtf(param, outGTF=gtfTempFn, width=width, fix="end",
            minTxLength=minTxLength, useTxIDs=useTxIDs)
}
