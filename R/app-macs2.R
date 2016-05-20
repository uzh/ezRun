###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMacs2 = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  dataset = input$meta
  #setwdNew(paste(basename(param$resultDir),output$getNames(),sep='_'))
  
  if (param$useControl){
    cmd = paste(MACS2, "callpeak -t", input$getFullPaths(param, "BAM"), 
                "-c", input$getFullPaths(param, "Control"), 
                "-B", opt,"-n", output$getNames())
    ezSystem(cmd)
    bedgraphFileTreat = paste0(output$getNames(), '_treat_pileup.bdg')
    bedgraphFileControl = paste0(output$getNames(), '_control_lambda.bdg')
    cmd = paste(MACS2, " bdgcmp -t", bedgraphFileTreat,
                "-c", bedgraphFileControl, "-o", paste0(output$getNames(),"_FE.bdg"), "-m FE")
    ezSystem(cmd)
    cmd = paste(BEDGRAPHBIGWIG, paste0(output$getNames(), "_FE.bdg"), param$ezRef@refChromSizesFile, paste0(output$getNames(), ".bw"))
    ezSystem(cmd)
    ezSystem("rm *.bdg")
  } else {
    cmd = paste(MACS2, "callpeak -t", input$getFullPaths(param, "BAM"), opt,"-n", output$getNames())
    ezSystem(cmd)
    createBigWig(input, output, param)
  }
  if (grepl('broad', opt)){
    ezSystem(paste("mv ",paste0(output$getNames(),"_peaks.broadPeak")," ",paste0(output$getNames(),"_peaks.bed")))
  } else {
    ezSystem(paste("mv ",paste0(output$getNames(),"_peaks.narrowPeak")," ",paste0(output$getNames(),"_peaks.bed")))
  }
  peakBedFile = paste0(output$getNames(),"_peaks.bed")
  cmd = paste(BEDTOOLS2, " getfasta -fi", param$ezRef["refFastaFile"],
            " -bed ", peakBedFile, " -name -fo ", paste0(output$getNames(), "_peaks.fa")) 
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
                  appDefaults <<- rbind(useControl=ezFrame(Type="logical",  DefaultValue="TRUE",	Description="should control samples be used"))
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
  requireNamespace("rtracklayer")
  requireNamespace("GenomicRanges")
  
  gtfFile = param$ezRef@refFeatureFile
  gtf = rtracklayer::import(gtfFile)
  idx = gtf$type =='gene'
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
  requireNamespace("rtracklayer")
  requireNamespace("GenomicRanges")
  requireNamespace("GenomicAlignments")
  if (param$paired){
    aligns = readGAlignmentPairs(file=input$getFullPaths(param, "BAM"))
  } else {
    aligns = readGAlignments(file=input$getFullPaths(param, "BAM"))
  }
  cov = coverage(aligns)
  export(cov,paste0(output$getNames(), ".bw"), format="bigWig")
}
