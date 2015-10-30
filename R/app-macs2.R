###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Macs2
##' @seealso \code{\link{EzAppMacs2}}
ezMethodMacs2 = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  dataset = input$meta
  #setwdNew(paste(basename(param$resultDir),output$getNames(),sep='_'))
  
  if(param$useControl){
    cmd = paste(MACS2, "callpeak -t", input$getFullPaths(param, "BAM"), 
                "-c", input$getFullPaths(param, "Control"), 
                "-B", opt,"-n", output$getNames())
    ezSystem(cmd)
    bedgraphFileTreat = paste0(output$getNames(), '_treat_pileup.bdg')
    bedgraphFileControl = paste0(output$getNames(), '_control_lambda.bdg')
    cmd = paste0(MACS2, "bdgcmp -t ", bedgraphFileTreat,
                "-c", bedgraphFileControl, "-o", paste(output$getNames(),"_FE.bdg"), "-m FE")
    ezSystem(cmd)
    cmd = paste0(BEDGRAPHBIGWIG, paste(output$getNames(), "_FE.bdg"), getRefChromSizesFile(param), paste(output$getNames(), ".bw"))
    ezSystem(cmd)
    ezSystem("rm *.bdg")
  } else {
    cmd = paste(MACS2, "callpeak -t", input$getFullPaths(param, "BAM"), opt,"-n", output$getNames())
    ezSystem(cmd)
    createBigWig(input,output,param)
  }
  if(grepl('broad', opt)){
    ezSystem(paste("mv ",paste0(output$getNames(),"_peaks.broadPeak")," ",paste0(output$getNames(),"_peaks.bed")))
  } else {
    ezSystem(paste("mv ",paste0(output$getNames(),"_peaks.narrowPeak")," ",paste0(output$getNames(),"_peaks.bed")))
  }
  peakBedFile = paste0(output$getNames(),"_peaks.bed")
  cmd = paste(BEDTOOLS2, " getfasta -fi", file.path(GENOMES_ROOT, dirname(dirname(input$getColumn("refBuild"))), "Sequence/WholeGenomeFasta/genome.fa"),
            " -bed ", peakBedFile, " -name -fo ", paste0(output$getNames(), "_peaks.fa")) 
  ezSystem(cmd)
  annotatePeaks(input,output,param)  
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMacs2()
##' @seealso \code{\link{ezMethodMacs2}}
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

annotatePeaks = function(input=NA,output=NA,param=NA) {
  peakFile = paste0(output$getNames(), '_peaks.xls')
  data = read.table(peakFile,sep='\t',check.names=F,header=T)
  if (nrow(data) == 0){
    return(NULL)
  }
  data = data[order(data$chr,data$start),]
  
  require(ChIPpeakAnno);require(rtracklayer);require(GenomicRanges);require(Rsamtools)
  gtfFile = param$ezRef@refFeatureFile
  gtf = import(gtfFile, asRangedData=FALSE)
  idx = gtf$type =='gene'
  if(!any(idx)){
    idx = gtf$type =='start_codon'
  }
  gtf = gtf[idx]
  if(grepl('gtf',gtfFile)){
    names_gtf = make.unique(gtf$'gene_id')
  } else {
    names_gtf = make.unique(gtf$'ID')
  }
  names(gtf) = names_gtf
  annoRD = as(gtf, "RangedData")
  peaksRD = RangedData(space=data$chr, IRanges(data$start, data$end), strand=rep('*',nrow(data)))
  rownames(peaksRD) = data$name
  annotatedPeaks = as.data.frame.RangedData(annotatePeakInBatch(peaksRD,AnnotationData = annoRD,output='nearestStart',multiple=FALSE,FeatureLocForDistance='TSS'))
  annotatedPeaks = annotatedPeaks[,c("peak","strand","feature","start_position","end_position","insideFeature","distancetoFeature")]
  annotatedPeaks = merge(data,annotatedPeaks,by.x='name',by.y='peak',all.x=T)
  localAnnotation = ezRead.table(param$ezRef@refAnnotationFile)
  localAnnotation = unique(localAnnotation[,grep('^gene_id$|^description$|name$|symbol$|^type$',colnames(localAnnotation),ignore.case=TRUE)])
  if(!is.null(ncol(localAnnotation))){
    annotatedPeaks = merge(annotatedPeaks,localAnnotation,by.x='feature',by.y='gene_id',all.x=T)
  } 
  annotatedPeaks = annotatedPeaks[order(annotatedPeaks$"-log10(pvalue)",decreasing=T),]
  colnames(annotatedPeaks) = gsub('-log10','_-log10',colnames(annotatedPeaks))
  write.table(annotatedPeaks,peakFile,sep='\t',row.names=F,quote=F)
}

createBigWig = function(input=NA,output=NA,param=NA){
  require(IRanges); require(rtracklayer);require(GenomicRanges);require(GenomicAlignments)
  if(param[['paired']]){
    aligns = readGAlignmentPairs(file=input$getFullPaths(param, "BAM"))
  } else {
    aligns = readGAlignments(file=input$getFullPaths(param, "BAM"))
  }
  cov = coverage(aligns)
  export(cov,paste0(output$getNames(), ".bw"), format="bigWig")
}
