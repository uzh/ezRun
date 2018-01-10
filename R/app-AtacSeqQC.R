###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodATACSeqQC <- function(input, output, param){
  require(genomation)
  require(GenomicRanges)
  require(chromstaR)
  require(TEQC)
  require(ATACseqQC)
  ###TODO: process in parallel, add more plots from ATACseqQC package and remove hardcoding parameters
  
  drawHeatmaps = T
  ctcfPeakFile = param$ctcfPeakFile
  ctcf =  read.table(ctcfPeakFile, sep = '\t', stringsAsFactors = F, skip=1, header = F)
  #1. Get BAM-Files:
  bamFiles = input$getFullPaths("BAM")
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
  
  for (i in 1:length(bamFiles)){
    #get granges from BAM-File:
    allReads <- readBamFileAsGRanges(bamFiles[i], chromosomes=NULL, pairedEndReads = param$paired, max.fragment.width = 1000, min.mapq = 10, remove.duplicate.reads = FALSE)
    readsByFragmentLength = list()
    readsByFragmentLength[[1]] <- allReads[width(allReads)< 100,]
    readsByFragmentLength[[2]] <- allReads[width(allReads)>= 100 & width(allReads)<= 300,]
    readsByFragmentLength[[3]] <- allReads[width(allReads)> 300,]
    names(readsByFragmentLength) = c('subNucleosomal','monoNucleosomal','multiNucleosomal')
    sm = list()
    sm[[1]] = ScoreMatrixBin(readsByFragmentLength[[1]], windows = tss.regions, bin.num = 50, type = 'bam')
    sm[[2]] = ScoreMatrixBin(readsByFragmentLength[[2]], windows = tss.regions, bin.num = 50, type = 'bam')
    sm[[3]] = ScoreMatrixBin(readsByFragmentLength[[3]], windows = tss.regions, bin.num = 50, type = 'bam')
    names(sm) = c('subNucleosomal','monoNucleosomal','multiNucleosomal')
    if(drawHeatmaps){
      for (j in 1:length(sm)){
        png(paste0('heatmap_TSS_',basename(bamFiles[i]),'_',names(sm)[j],'.png'),width = 600, height = 500)
        heatMatrix(sm[[j]], xcoords = c(-1000, 1000),winsorize = c(0,99), order=T, xlab='TSS', main=paste(basename(bamFiles[i]), names(sm)[j]))
        dev.off()
      }
    }
    for (j in 1:length(sm)){
      png(paste0('lineplot_TSS_',basename(bamFiles[i]),'_',names(sm)[j],'.png'),width = 600, height = 500)
      plotMeta(sm[[j]], xcoords = c(-2000, 2000),xlab='TSS',main=paste(basename(bamFiles[i]), names(sm)[j]))
      dev.off()
    }
  }
  
  #Get ctcf regions
  ctcf  = ctcf[ctcf$V5=='+', ]
  ctcf$V2  = gsub('chr','',ctcf$V2)
  ctcf = ctcf[nchar(ctcf$V2) <3,]
  ctcf = ctcf[1:10000,]
  ctcf = data.frame(chrom=ctcf$V2,start=ctcf$V3,end=ctcf$V4, strand=ctcf$V5)
  
  ctcf.regions <- GRanges(seqnames = Rle(ctcf$chrom),
                          strand = Rle(ctcf$strand),
                          ranges = IRanges(ctcf$start-2000,ctcf$end+2000))
  
  for (i in 1:length(bamFiles)){
    sm = list()
    sm[[1]] = ScoreMatrixBin(readsByFragmentLength[[1]], windows = ctcf.regions, bin.num = 50, type = 'bam')
    sm[[2]] = ScoreMatrixBin(readsByFragmentLength[[2]], windows = ctcf.regions, bin.num = 50, type = 'bam')
    sm[[3]] = ScoreMatrixBin(readsByFragmentLength[[3]], windows = ctcf.regions, bin.num = 50, type = 'bam')
    names(sm) = c('subNucleosomal','monoNucleosomal','multiNucleosomal')
    
    if(drawHeatmaps){
      for (j in 1:length(sm)){
        png(paste0('heatmap_CTCF_',basename(bamFiles[i]),'_',names(sm)[j],'.png'),width = 600, height = 500)
        heatMatrix(sm[[j]], xcoords = c(-1000, 1000),winsorize = c(0,99), order=T, xlab='CTCF-Motif', main=paste(basename(bamFiles[i]), names(sm)[j]))
        dev.off()
      }
    }
    for (j in 1:length(sm)){
      png(paste0('lineplot_CTCF_',basename(bamFiles[i]),'_',names(sm)[j],'.png'),width = 600, height = 500)
      plotMeta(sm[[j]], xcoords = c(-2000, 2000),xlab='CTCF-Motif',main=paste(basename(bamFiles[i]), names(sm)[j]))
      dev.off()
    }
  }
  
  chrM_content = vector(mode = 'numeric', length = length(bamFiles))
  names(chrM_content) = gsub('.bam', '', basename(bamFiles))
  shortFragments = vector(mode = 'numeric', length = length(bamFiles))
  names(shortFragments) = names(chrM_content)
  
  for (i in 1:length(bamFiles)){
    reads <- get.reads(bamFiles[i], filetype="bam")
    readpairs <- reads2pairs(reads)
    shortFragments[i] = 100 * table(width(readpairs) < param$minInsert)['TRUE']/nrow(readpairs)
    png(paste0('InsertSize',basename(bamFiles[i]),'.png'),width = 600, height = 500)
    insert.size.hist(readpairs, breaks=50, legendpos = 'topright', main = basename(bamFiles[i]))
    dev.off()
    
    
    png(paste0('fragmentSize_',basename(bamFiles[i]),'.png'),width = 600, height = 500)
      fragSize <- fragSizeDist(bamFiles[i], samples[i])
    dev.off()
    
    
    res = table(space(readpairs))
    chrM_content[i] = 100 * res[grep('^M',names(res))]/sum(res)
  }
  
  png(paste0('ChrM_content.png'),400,400)
  par(mar=c(14.1,7.1,4.1,2.1))
  bp = barplot(chrM_content, las = 2, ylab='ChrM content in %', main='ChrM Reads')
  dev.off()
  
  png(paste0('ShortFragmentFraction.png'),400,400)
  par(mar=c(14.1,7.1,4.1,2.1))
  bp = barplot(shortFragments, las = 2, ylab='ShortFragments below 100nt in %', main='Short Fragments')
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

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodATACSeqQC(input=NA, output=NA, param=NA)
##' @templateVar htmlArg, htmlFile="00index.html"
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
                  #appDefaults <<- rbind(disp_plot = ezFrame(Type="character", DefaultValue="dispersion_estimate_plot", Description="which test method in DEXSeq to use: deseq2"))
                }
              ))
