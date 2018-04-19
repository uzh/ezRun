###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMEME = function(input=NA, output=NA, param=NA){
  opt = param$cmdOptions
  db = param$motifDB
  sampleName = input$getNames()
  if(param$filterPeaks){
    param$relPos = c('upstream','overlapStart', 'inside', 'includeFeature')
    peaks = ezRead.table(input$getFullPaths("CalledPeaks"), row.names = NULL)
    peaks = peaks[abs(peaks$distancetoFeature) <=  param$distance, ]
    peaks = peaks[peaks$insideFeature %in%  param$relPos,]
    peaks = peaks[abs(peaks$fold_enrichment) >=  param$minFold, ]
    toRemove = which(duplicated(peaks$name))
    
    if (length(toRemove)>0){
      peaks = peaks[-toRemove,] }
    ezWrite.table(peaks, basename(input$getFullPaths("CalledPeaks")),  row.names = F)
    
    #CreateBedFile:
    bed = ezRead.table(input$getFullPaths("BED"), row.names = NULL, header = F)
    bed = bed[bed$V4 %in% peaks$name,]
    bed = bed[order(bed$V4),]
    peaks = peaks[order(peaks$name),]
    bed$V4 = paste0(bed$V4,'__',paste(peaks$chr,peaks$start,peaks$end,peaks$gene_name,peaks$distancetoFeature,sep='_'))
    bedFileName = basename(input$getFullPaths("BED"))
    ezWrite.table(bed, bedFileName, row.names = F, col.names = F)
    
    refFile = param$ezRef["refFastaFile"]
    peakSeqFile = paste0(bedFileName, "_peaks.fa")
    cmd = paste("/usr/local/ngseq/bin/bedtools", " getfasta -fi", refFile,
                  " -bed ", bedFileName, " -name -fo ", peakSeqFile)
    system(cmd)
    cmd = paste("meme-chip -oc",sampleName,"-time 300 -order 1", db, opt, peakSeqFile)
  } else{
    cmd = paste("meme-chip -oc",sampleName,"-time 300 -order 1", db, opt, 
              input$getFullPaths("PeakSequences"))
  }
  
  ezSystem(cmd)
  setwd(sampleName)
  ezSystem(paste('mv meme-chip.html ',paste0(sampleName,"_meme-chip.html")))
  setwd('..')
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodMEME(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppMEME <-
  setRefClass("EzAppMEME",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMEME
                  name <<- "EzAppMEME"
                }
              )
  )
