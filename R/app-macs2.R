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
      require(Biostrings)
      gsize <- sum(as.numeric(fasta.seqlengths(param$ezRef["refFastaFile"])))
      gsize <- round(gsize * 0.8)
      message("Use calculated gsize: ", gsize)
    }
    opt <- paste(opt, "-g", gsize)
  }
  
  if(param$mode == "ChIP-seq"){
    ## --extsize: extend reads in 5'->3' direction to fix-sized fragments.
    if(!grepl("--extsize", opt)){
      opt <- paste(opt, "--extsize 147")
    }
    bamFile <- input$getFullPaths("BAM")
    outBam <- basename(output$getFullPaths("BAM"))
    dupBam(inBam=bamFile, outBam=outBam, operation="remove",
           cores=param$cores)
    if (isTRUE(param$useControl)){
      if(!grepl("Control", input$colNames))
        stop("Control is not available when paramter useControl is true.")
      
      cmd = paste("macs2", "callpeak -t", input$getFullPaths("BAM"), 
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
                  paste0(output$getNames(), ".bw"))
      ezSystem(cmd)
      ezSystem("rm *.bdg")
    } else {
      cmd = paste("macs2", "callpeak -t", input$getFullPaths("BAM"), opt,
                  "-n", output$getNames())
      ezSystem(cmd)
      createBigWig(input, output, param)
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
    
    cmd = paste("macs2", "callpeak -t", basename(output$getFullPaths("BAM")),
                opt, "-n", output$getNames())
    ezSystem(cmd)
    createBigWig(input, output, param)
  }else{
    stop("MACS2 only supports ChIP-seq or ATAC-seq data.")
  }
  if (grepl('broad', opt)){
    file.rename(from=paste0(output$getNames(),"_peaks.broadPeak"),
                to=paste0(output$getNames(),"_peaks.bed"))
  } else {
    file.rename(from=paste0(output$getNames(),"_peaks.narrowPeak"),
                to=paste0(output$getNames(),"_peaks.bed"))
  }
  peakBedFile = paste0(output$getNames(),"_peaks.bed")
  cmd = paste("bedtools", " getfasta -fi", param$ezRef["refFastaFile"],
              " -bed ", peakBedFile, " -name -fo ",
              paste0(output$getNames(), "_peaks.fa"))
  ezSystem(cmd)
  annotatePeaks(input, output, param)
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
                  appDefaults <<- rbind(useControl=ezFrame(Type="logical", DefaultValue="TRUE",	Description="should control samples be used"))
                }
              )
  )

##' @title Annotates peaks
##' @description Annotates peaks and writes them into a separate table.
##' @template input-template
##' @template output-template
##' @param param a list of parameters to extract the \code{ezRef@@refFeatureFile} and the \code{ezRef@@refAnnotationFile} from.
##' @template roxygen-template
annotatePeaks = function(input=NA, output=NA, param=NA) {
  peakFile = paste0(output$getNames(), '_peaks.xls')
  data = ezRead.table(peakFile, comment.char = "#", row.names = NULL)
  if (nrow(data) == 0){
    return(NULL)
  }
  data = data[order(data$chr,data$start),]
  require("rtracklayer")
  require("GenomicRanges")
  
  gtfFile = param$ezRef@refFeatureFile
  gtf = rtracklayer::import(gtfFile)
  idx = gtf$type == 'gene'
  if(!any(idx)){
    idx = gtf$type =='start_codon'
  }
  gtf = gtf[idx]
  if(grepl('gtf$',gtfFile)){
    names_gtf = make.unique(gtf$'gene_id')
  } else {
    names_gtf = make.unique(gtf$'ID')
  }
  names(gtf) = names_gtf
  annoRD = as(gtf, "RangedData")
  peaksRD = RangedData(space=data$chr, IRanges(data$start, data$end), strand=rep('*',nrow(data)))
  rownames(peaksRD) = data$name
  annotatedPeaks = as.data.frame(ChIPpeakAnno::annotatePeakInBatch(peaksRD,AnnotationData = annoRD,output='nearestStart',multiple=FALSE,FeatureLocForDistance='TSS'))
  annotatedPeaks = annotatedPeaks[,c("peak","strand","feature","start_position","end_position","insideFeature","distancetoFeature")]
  annotatedPeaks = merge(data,annotatedPeaks,by.x='name',by.y='peak',all.x=T)
  localAnnotation = ezRead.table(param$ezRef@refAnnotationFile)
  localAnnotation = unique(localAnnotation[,grep('^gene_id$|^description$|name$|symbol$|^type$',colnames(localAnnotation),ignore.case=TRUE)])
  if(!is.null(ncol(localAnnotation))){
    annotatedPeaks = merge(annotatedPeaks,localAnnotation,by.x='feature',by.y='gene_id',all.x=T)
  }
  annotatedPeaks = annotatedPeaks[order(annotatedPeaks$"-log10(pvalue)",decreasing=T),]
  colnames(annotatedPeaks) = gsub('-log10','_-log10',colnames(annotatedPeaks)) ## TODO: explain why this is done? and why gsub and not sub?
  ezWrite.table(annotatedPeaks,peakFile, row.names=F)
}

##' @title Creates a bigwig file
##' @description Creates and exports a bigwig file.
##' @template input-template
##' @template output-template
##' @param param a list of parameters to extract \code{paired} from.
##' @template roxygen-template
createBigWig = function(input=NA, output=NA, param=NA){
  require("rtracklayer")
  require("GenomicRanges")
  require("GenomicAlignments")
  if (param$paired){
    aligns = readGAlignmentPairs(file=input$getFullPaths("BAM"))
  } else {
    aligns = readGAlignments(file=input$getFullPaths("BAM"))
  }
  cov = coverage(aligns)
  export.bw(cov, paste0(output$getNames(), ".bw"))
}
