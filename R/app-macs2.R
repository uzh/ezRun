###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMacs2 = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  if(param$paired)
    opt = paste(opt,'-f BAMPE')
  ## With BAMPE file, --shift cannot be set.
  dataset = input$meta

  ## -g option: mappable genome size
  if(!grepl("-g", opt)){
    gsizes <- c("Homo sapiens (human)"="hs",
                "Mus musculus (mouse)"="mm",
                "Caenorhabditis elegans (worm)"="ce",
                "Drosophila melanogaster (fruitfly)"="dm")
    isGsize <- grepl(input$getColumn("Species"), names(gsizes),
                     ignore.case = TRUE)
    if(any(isGsize)){
      message("Use predefined gsize: ", gsizes[isGsize])
      gsize <- gsizes[isGsize]
    }else{
      gsize <- sum(as.numeric(fasta.seqlengths(param$ezRef["refFastaFile"])))
      gsize <- round(gsize * 0.8)
      message("Use calculated gsize: ", gsize)
    }
    opt <- paste(opt, "-g", gsize)
  }
  ## --keep-dup: behavior towards duplicate tags at the exact same location
  if(!grepl("--keep-dup", opt)){
    opt <- paste(opt, "--keep-dup all")
  }

  if(param$mode == "ChIP-seq"){
    ## --extsize: extend reads in 5'->3' direction to fix-sized fragments when model building is deactivated.
    ## This size is taken from sushi app.
    if(grepl("--nomodel", opt) && !grepl("--extsize", opt)){
      opt <- paste(opt, "--extsize 147")
    }
    bamFile <- input$getFullPaths("BAM")
    outBam <- basename(output$getColumn("BAM"))
    if(param$removeDuplicates){
      dupBam(inBam=bamFile, outBam=outBam, operation="remove", ram = param$ram)
    } else {
      file.copy(from=bamFile, to=outBam, overwrite=TRUE)
      Rsamtools::indexBam(outBam)
    }

    if (isTRUE(param$useControl)){
      if(!any(grepl("Control", input$colNames, ignore.case = TRUE)))
        stop("The parameter 'useControl' is 'true' but no column named 'Control [File]' is available.")

      cmd = paste("macs2", "callpeak -t", outBam,
                  "-c", input$getFullPaths("Control"),
                  "-B", opt,"-n", output$getNames())
      ezSystem(cmd)
      bedgraphFileTreat = paste0(output$getNames(), '_treat_pileup.bdg')
      bedgraphFileControl = paste0(output$getNames(), '_control_lambda.bdg')
      cmd = paste("macs2", " bdgcmp -t", bedgraphFileTreat,
                  "-c", bedgraphFileControl, "-o",
                  paste0(output$getNames(),"_FE.bdg"), "-m FE")
      ezSystem(cmd)
      bdgSorted = "sorted.bdg"
      ezSystem(paste("bedSort", paste0(output$getNames(), "_FE.bdg"), bdgSorted))
      cmd = paste("bedGraphToBigWig", bdgSorted, param$ezRef@refChromSizesFile,
                  paste0(output$getNames(), "_processed.bw"))
      ezSystem(cmd)
      ezSystem("rm *.bdg")
    } else {
      cmd = paste("macs2", "callpeak -t", outBam, opt,
                  "-n", output$getNames())
      ezSystem(cmd)
      bam2bw(file=outBam, destination=basename(output$getColumn("BigWigFile")),
             paired=param$paired,
             method="deepTools", cores=param$cores)
    }
  }else if(param$mode == "ATAC-seq"){
    if(!param$paired)
      stop("For ATAC-seq, we only support paired-end data.")

    ## --extsize: extend reads in 5'->3' direction to fix-sized fragments.
    if(!grepl("--extsize", opt)){
      ## https://github.com/taoliu/MACS/issues/145
      opt <- paste(opt, "--extsize 200")
    }

    ## Preprocess ATAC-seq bam file
    atacBamProcess(input=input, output=output, param=param)

    cmd = paste("macs2", "callpeak -t", basename(output$getColumn("BAM")),
                opt, "-n", output$getNames())
    ezSystem(cmd)
    bam2bw(file=basename(output$getColumn("BAM")),
           destination=basename(output$getColumn("BigWigFile")),
           paired=param$paired,
           method="deepTools", cores=param$cores)
  }else{
    stop("MACS2 only supports ChIP-seq or ATAC-seq data.")
  }

  peakBedFile = basename(output$getColumn("BED"))
  if (grepl('broad', opt)){
    file.rename(from=paste0(output$getNames(),"_peaks.broadPeak"),
                to=peakBedFile)
  } else {
    file.rename(from=paste0(output$getNames(),"_peaks.narrowPeak"),
                to=peakBedFile)
  }
  peakSeqFile = basename(output$getColumn("PeakSequences"))
  cmd = paste("bedtools", " getfasta -fi", param$ezRef["refFastaFile"],
              " -bed ", peakBedFile, " -name -fo ",peakSeqFile)
  ezSystem(cmd)
  peakXlsFile <- basename(output$getColumn("CalledPeaks"))
  annotatePeaks(peakXlsFile, peakSeqFile, param)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMacs2(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppMacs2 <-
  setRefClass("EzAppMacs2",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMacs2
                  name <<- "EzAppMacs2"
                  appDefaults <<- rbind(useControl=ezFrame(Type="logical", DefaultValue="TRUE",	Description="should control samples be used"),
                                        shiftATAC=ezFrame(Type="logical", DefaultValue="FALSE",	Description="should all reads aligning to + strand were offset by +4bp, all reads aligning to the - strand are offset -5 bp"),
                                        annotatePeaks=ezFrame(Type="logical", DefaultValue="TRUE",	Description="use gtf to annotate peaks"))
                }
              )
  )

##' @title Annotates peaks
##' @description Annotates peaks and writes them into a separate table.
##' @template input-template
##' @template output-template
##' @param param a list of parameters to extract the \code{ezRef@@refFeatureFile} and the \code{ezRef@@refAnnotationFile} from.
##' @template roxygen-template
annotatePeaks = function(peakFile, peakSeqFile, param) {
  require(rtracklayer)
  require(ChIPpeakAnno)
  require(ShortRead)
  require(ChIPseeker)
  require(GenomicFeatures)
  
  data <- c()
  tryCatch(expr = {data = ezRead.table(peakFile, comment.char = "#", row.names = NULL)}, 
           error = function(e){message(paste("No peaks detected. Skip peak annotation"))})
  if (is.null(data)){
    return(NULL)
  }
  data = data[order(data$chr,data$start), ]

  if(!param$annotatePeaks){
      ezWrite.table(data, peakFile, row.names = F)
      return('done')
  }
  
  gtfFile = param$ezRef@refFeatureFile
  myTxDB <- makeTxDbFromGFF(file=gtfFile, format='gtf')
  
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
  peaksRD = makeGRangesFromDataFrame(data, keep.extra.columns = TRUE)
  names(peaksRD) = mcols(peaksRD)$name
  annot_ChIPseeker <- annotatePeak(peaksRD, TxDb=myTxDB, tssRegion=c(-1000, 1000), verbose=FALSE)
  annotatedPeaks <- annotatePeakInBatch(peaksRD,
                                        AnnotationData = gtf,
                                        output='nearestStart',
                                        multiple=FALSE,
                                        FeatureLocForDistance='TSS')
  annotatedPeaks = as.data.frame(annotatedPeaks)
  annotatedPeaks = annotatedPeaks[ , c("peak", "feature", "feature_strand",
                                       "start_position", "end_position",
                                       "insideFeature", "distancetoFeature")]
  colnames(annotatedPeaks) = c("peak", "feature", "feature_strand",
                               "feature_start", "feature_end",
                               "insideFeature", "distancetoFeature")

  annotatedPeaks = merge(data, annotatedPeaks,by.x='name', by.y='peak', all.x=T)
  localAnnotation <- ezFeatureAnnotation(param, dataFeatureType="gene")
  localAnnotation = unique(localAnnotation[, grep('^gene_id$|^description$|name$|symbol$|^type$',colnames(localAnnotation),ignore.case=TRUE)])
  if(!is.null(ncol(localAnnotation))){
    annotatedPeaks = merge(annotatedPeaks, localAnnotation, by.x='feature',
                           by.y='gene_id',all.x=T)
  }
  colnames(annotatedPeaks) = gsub('-log10','_-log10', colnames(annotatedPeaks))
  seqs = readDNAStringSet(peakSeqFile)
  dustyScores = data.frame(ID = names(seqs), dustyScore_peakSequence = dustyScore(seqs), stringsAsFactors = F)
  annotatedPeaks = merge(annotatedPeaks, dustyScores, by.x = 'name', by.y = 'ID', all.x = TRUE)
  
  annot_ChIPseeker <- data.frame(annot_ChIPseeker@anno)
  keepCol_ChIPSeeker <- c("name", "annotation", "geneId", "transcriptId", "distanceToTSS", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand")
  annot_ChIPseeker <- annot_ChIPseeker[,keepCol_ChIPSeeker]
  colnames(annot_ChIPseeker)[2:ncol(annot_ChIPseeker)] <- paste0('ChIPSeeker_', colnames(annot_ChIPseeker)[2:ncol(annot_ChIPseeker)])
  annotatedPeaks <- merge(annotatedPeaks, annot_ChIPseeker, by.x = 'name', by.y = 'name', all.x = TRUE)
  annotatedPeaks <- annotatedPeaks[!duplicated(annotatedPeaks$name),]
  annotatedPeaks <- annotatedPeaks[order(annotatedPeaks[['_-log10(pvalue)']],decreasing=TRUE),]
  
  writexl::write_xlsx(annotatedPeaks, peakFile)
  return('done')
}

### import Macs2's BED6+4 file: narrowPeak or broadPeak
import.Macs2Peaks <- function(con){
  bed <- ezRead.table(con, header=FALSE, row.names=NULL)
  colnames(bed) <- c("chr", "start", "end", "name", "score", "strand",
                     "fold-change", "p-value", "q-value", "summit")
  bed <- transform(bed, start=start+1)
  bed <- makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE)
  return(bed)
}
