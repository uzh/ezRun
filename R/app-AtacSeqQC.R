###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodATACSeqQC <- function(input, output, param, htmlFile="00index.html"){
  require(genomation)
  require(GenomicRanges)
  require(chromstaR)
  require(TEQC)
  require(ATACseqQC)
  
  setwdNew(basename(output$getColumn("Report")))
  ###TODO: process in parallel, add more plots from ATACseqQC package and remove hardcoding parameters, create multisample TSS LinePlot
  drawHeatmaps = T
  
  #1. Get BAM-Files:
  bamFiles = input$getFullPaths("BAM")
  bamFileList = as.list(bamFiles)
  samples = input$getNames()
  dataset = input$meta
  
  #2. Get gtf:
  gtfFile = param$ezRef@refFeatureFile
  gtf = gffToGRanges(gtfFile, filter = 'gene')
  df_gtf = as.data.frame(gtf)
  df_gtf = df_gtf[df_gtf$strand=='+',] #Use only +-strand features for QC
  tss.regions <- GRanges(seqnames = Rle(df_gtf$seqnames),
                         strand = Rle(df_gtf$strand),
                         ranges = IRanges(df_gtf$start-2000,df_gtf$start+2000))
  
  if(param$ctcfPeakFile != ''){
    ctcfPeakFile = param$ctcfPeakFile
    ctcf =  read.table(ctcfPeakFile, sep = '\t', stringsAsFactors = F, skip=1, header = F)
    ctcf  = ctcf[ctcf$V5=='+', ]
    ctcf$V2  = gsub('chr','',ctcf$V2)
    ctcf = ctcf[nchar(ctcf$V2) <3,]
    ctcf = ctcf[1:10000,]
    ctcf = data.frame(chrom=ctcf$V2,start=ctcf$V3,end=ctcf$V4, strand=ctcf$V5)
    ctcf.regions <- GRanges(seqnames = Rle(ctcf$chrom),
                          strand = Rle(ctcf$strand),
                          ranges = IRanges(ctcf$start-2000,ctcf$end+2000))
}
  
  multiScoreMatrixList = list()
  multiScoreMatrixListCTCF = list()
  for (i in 1:length(bamFiles)){
    #get granges from BAM-File:
    allReads <- readBamFileAsGRanges(bamFiles[i], chromosomes=NULL, pairedEndReads = param$paired, 
                                     max.fragment.width = 1000, min.mapq = param$minMapq, remove.duplicate.reads = FALSE)
    readsByFragmentLength = list()
    readsByFragmentLength[[1]] <- allReads[width(allReads)< param$minInsert,]
    readsByFragmentLength[[2]] <- allReads[width(allReads)>= param$minInsert & width(allReads)<= param$multiNuclInsert,]
    readsByFragmentLength[[3]] <- allReads[width(allReads)> param$multiNuclInsert,]
    names(readsByFragmentLength) = c('subNucleosomal','monoNucleosomal','multiNucleosomal')
    multiScoreMatrixList[[i]] = list() 
    multiScoreMatrixList[[i]][[1]] = ScoreMatrixBin(readsByFragmentLength[[1]], windows = tss.regions, bin.num = 50, type = 'bam')
    multiScoreMatrixList[[i]][[2]] = ScoreMatrixBin(readsByFragmentLength[[2]], windows = tss.regions, bin.num = 50, type = 'bam')
    multiScoreMatrixList[[i]][[3]] = ScoreMatrixBin(readsByFragmentLength[[3]], windows = tss.regions, bin.num = 50, type = 'bam')
    names(multiScoreMatrixList[[i]]) = c('subNucl','monoNucl','multiNucl')
    
    if(drawHeatmaps){
      for (j in 1:length(multiScoreMatrixList[[i]])){
        png(paste0('heatmap_TSS_',basename(bamFiles[i]),'_',names(multiScoreMatrixList[[i]])[j],'.png'),width = 600, height = 500, res  = 90)
        heatMatrix(multiScoreMatrixList[[i]][[j]], xcoords = c(-1*param$flankingRegion, param$flankingRegion),winsorize = c(0,99), order=T, xlab='TSS',
                   main = paste(basename(bamFiles[i]), names(multiScoreMatrixList[[i]])[j]))
        dev.off()
      }
    }
    
      x=new("ScoreMatrixList",list(multiScoreMatrixList[[i]][[1]],multiScoreMatrixList[[i]][[2]], multiScoreMatrixList[[i]][[3]]))
      png(paste0('lineplot_TSS_',basename(bamFiles[i]),'.png'),width = 700, height = 500, res = 90)
        par(mar=c(5.1,4.1,4.1,2.1))
        plotMeta(x, xcoords = c(-2000, 2000),xlab='TSS',main=samples[i], profile.names=names(multiScoreMatrixList[[i]]), lwd=3)
      dev.off()
    if(param$ctcfPeakFile != ''){  
    multiScoreMatrixListCTCF[[i]] = list()
    multiScoreMatrixListCTCF[[i]][[1]] = ScoreMatrixBin(readsByFragmentLength[[1]], windows = ctcf.regions, bin.num = 50, type = 'bam')
    multiScoreMatrixListCTCF[[i]][[2]] = ScoreMatrixBin(readsByFragmentLength[[2]], windows = ctcf.regions, bin.num = 50, type = 'bam')
    multiScoreMatrixListCTCF[[i]][[3]] = ScoreMatrixBin(readsByFragmentLength[[3]], windows = ctcf.regions, bin.num = 50, type = 'bam')
    names(multiScoreMatrixListCTCF[[i]]) = c('subNucl','monoNucl','multiNucl')
    
    if(drawHeatmaps){
      for (j in 1:length(multiScoreMatrixListCTCF[[i]])){
        png(paste0('heatmap_CTCF_',basename(bamFiles[i]),'_',names(multiScoreMatrixListCTCF[[i]])[j],'.png'), width = 600, height = 500, res  = 90)
        heatMatrix(multiScoreMatrixListCTCF[[i]][[j]], xcoords = c(-1*param$flankingRegion, param$flankingRegion), winsorize = c(0,99), order=T,
                   xlab='CTCF-Motif', main=paste(samples[i], names(multiScoreMatrixListCTCF[[i]])[j]),  cex.legend = 1.2)
        dev.off()
      }
    }
   
    x=new("ScoreMatrixList",list(multiScoreMatrixListCTCF[[i]][[1]],multiScoreMatrixListCTCF[[i]][[2]], multiScoreMatrixListCTCF[[i]][[3]]))
    png(paste0('lineplot_CTCF_',basename(bamFiles[i]),'.png'),width = 700, height = 500, res = 90)
    par(mar=c(5.1,4.1,4.1,2.1))
      plotMeta(x, xcoords = c(-2000, 2000),xlab='CTCF-Motif',main=samples[i], profile.names=names(multiScoreMatrixListCTCF[[i]]), lwd=3)
    dev.off()
    }
  }
  
  if(length(samples)>1){
    for (j in 1:3){
        if(param$ctcfPeakFile != ''){
      myList = list()
      for (k in 1:length(samples)){
        myList[[k]] = multiScoreMatrixListCTCF[[k]][[j]]
      }
      x=new("ScoreMatrixList",myList)
      png(paste0('lineplot_CTCF_multiSample_',names(multiScoreMatrixListCTCF[[1]])[j],'.png'),width = 700, height = 500, res = 90)
      plotMeta(x, xcoords = c(-2000, 2000),xlab='CTCF-Motif',
               main=paste('MultiSample CTCF',names(multiScoreMatrixListCTCF[[1]])[j]), 
               profile.names=samples, lwd=3)
      dev.off()
        }
      myList = list()
      for (k in 1:length(samples)){
        myList[[k]] = multiScoreMatrixList[[k]][[j]]
      }
      x=new("ScoreMatrixList",myList)
      png(paste0('lineplot_TSS_multiSample_',names(multiScoreMatrixList[[1]])[j],'.png'),width = 700, height = 500, res = 90)
      plotMeta(x, xcoords = c(-2000, 2000),xlab='TSS',
               main=paste('MultiSample TSS',names(multiScoreMatrixList[[1]])[j]), 
               profile.names=samples, lwd=3)
      dev.off()
    }
  }
  
  chrM_content = vector(mode = 'numeric', length = length(bamFiles))
  names(chrM_content) = gsub('.bam', '', basename(bamFiles))
  shortFragments = vector(mode = 'numeric', length = length(bamFiles))
  names(shortFragments) = names(chrM_content)
  
  
  for (i in 1:length(bamFiles)){
    res = makeFragmentSizePlot(bamFiles[i])
    shortFragments[i] = res[[1]]
    chrM_content[i] = res[[2]]
  }
  
  
  png(paste0('ChrM_content.png'),400,400, res = 80)
  par(mar=c(14.1,7.1,4.1,2.1))
  barplot(chrM_content, las = 2, ylab='ChrM content in %', main='ChrM Reads')
  dev.off()
  
  png(paste0('ShortFragmentFraction.png'),400,400, res = 80)
  par(mar=c(14.1,7.1,4.1,2.1))
  barplot(shortFragments, las = 2, ylab='ShortFragments below 100nt in %', main='Short Fragments')
  dev.off()
  
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "AtacSeqQC.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="AtacSeqQC.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  return("Success")
}

makeFragmentSizePlot = function(bamFile){
  require(GenomicAlignments)
  reads <- granges(readGAlignmentPairs(bamFile))
  sampleName = sub('.bam', '', basename(bamFile))
  
  png(paste0('fragmentSize_',sampleName,'.png'),width = 600, height = 500, res = 90)
    fragSizeDist(bamFile, sampleName)
  dev.off()
  shortFragmentFraction = 100 * table(width(reads) < as.numeric(param$minInsert))['TRUE']/length(reads)
  res = table(seqnames(reads))
  chrM_Fraction = 100 * res[grep('^M',names(res))]/sum(res)
  return(list(shortFragmentFraction, chrM_Fraction))
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodATACSeqQC(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run QC on bamFiles from ATAC-Seq
EzAppATACSeqQC <-
  setRefClass(Class = "EzAppATACSeqQC",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodATACSeqQC
                  name <<- "EzAppATACSeqQC"
                  appDefaults <<- rbind(minInsert = ezFrame(Type="numeric", DefaultValue=100, Description="define sub-nucleosomal fragments threshold"),
                                        minMapq = ezFrame(Type="numeric", DefaultValue=10, Description="define min. MappingQuality"),
                                        multiNuclInsert = ezFrame(Type="numeric", DefaultValue=250, Description="define multi-nucleosomal fragments threshold"),
                                        flankingRegion = ezFrame(Type="numeric", DefaultValue=1000, Description="define upstream/downstream region of target for heatmap"))
                }
              ))
