###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCNVnator = function(input=NA, output=NA, param=NA){
  bamFile <- input$getFullPaths("BAM")
  cmd = paste(CNVNATOR, '-root out.root -unique -tree', bamFile)
  ezSystem(cmd)
  chromosomeFolder <- sub('WholeGenomeFasta/genome.fa', 'Chromosomes', param$ezRef@refFastaFile)
  cmd = paste(CNVNATOR, '-root out.root -his', param$binSize, '-d', chromosomeFolder)
  ezSystem(cmd)
  cmd = paste(CNVNATOR, '-root out.root -stat', param$binSize)
  ezSystem(cmd)
  cmd = paste(CNVNATOR, '-root out.root -partition', param$binSize)
  ezSystem(cmd)
  outFile = sub('.bam', '.txt', basename(bamFile))
  cmd = paste(CNVNATOR, '-root out.root -call', param$binSize, '>', outFile)
  ezSystem(cmd)
  annotateCNVs(outFile, param)
  }


annotateCNVs = function(file, param){
  require(GenomicRanges)
  gtfFile = param$ezRef@refFeatureFile
  gtf = rtracklayer::import(gtfFile)
  idx = gtf$type == 'gene' & gtf$source == 'protein_coding'
  gtf = gtf[idx]
  
  data = ezRead.table(file, header = F, row.names = NULL)
  colnames(data) = c('Event','Region','Width','normalized_RD', 'e-val1', 'e-val2','e-val3','e-val4', 'q0_fraction')
  data = data[data[['e-val1']] <= param$maxEVal, ]
  data = data[data[['e-val2']] <= param$maxEVal, ]
  data = data[data[['e-val3']] <= param$maxEVal, ]
  data = data[data[['e-val4']] <= param$maxEVal, ]
  data = data[data[['q0_fraction']] <= param$maxQ0, ]
  
  region = limma::strsplit2(data$Region, ':')
  coord = limma::strsplit2(region[,2],'-')
  pos = data.frame(CHR = region[,1], start = coord[,1], end = coord[,2], stringsAsFactors = F)
  segmentRanges = makeGRangesFromDataFrame(pos)
  names(segmentRanges) = rownames(data)
  
  olaps = findOverlaps(gtf, segmentRanges)
  f1 <- factor(subjectHits(olaps),
               levels=seq_len(subjectLength(olaps)))
  featureInCNVs  <- splitAsList(mcols(gtf)[['gene_id']][queryHits(olaps)], f1)
  genesInCNVs  <- splitAsList(mcols(gtf)[['gene_name']][queryHits(olaps)], f1)
  
  data = cbind(data, pos)
  data[['FeaturesInCnv']] = unlist(sapply(as.vector(featureInCNVs),paste, collapse=','))
  data[['GenesInCnv']] = unlist(sapply(as.vector(genesInCNVs),paste, collapse=','))
  
  data = data[data$FeaturesInCnv != '', ]
  ezWrite.table(data, sub('.txt','_CNV.txt',file), row.names = F)
}

##' @template app-template
##' @templateVar method ezMethodCNVNator(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppCNVnator <-
  setRefClass("EzAppCNVnator",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCNVnator
                  name <<- "ezMethodCNVnator"
                  appDefaults <<- rbind(binSize=ezFrame(Type="numeric", DefaultValue="1000",	Description="set binSize"),
                                        maxEVal=ezFrame(Type="numeric", DefaultValue="0.01",	Description="set maximum E-Value for CNV filtering"),
                                        maxQ0=ezFrame(Type="numeric", DefaultValue="0.1",	Description="set maximum fraction of low quality read mapping in CNV"))
                }
              )
  )

