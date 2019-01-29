###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template app-template
##' @templateVar method ezMethodTeqc(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppTeqc <-
  setRefClass("EzAppTeqc",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodTeqc
                  name <<- "EzAppTeqc"
                  appDefaults <<- rbind(designFile=ezFrame(Type="character",  DefaultValue="",  Description="file describing the regions selected by the enrichment kit"),
                                    covUniformityPlot=ezFrame(Type="logical", DefaultValue="TRUE", Description="generate plots for coverage uniformity?"),
                                    covTargetLengthPlot=ezFrame(Type="logical", DefaultValue="TRUE", Description="generate plots for coverage vs. target length"),
                                    duplicatesPlot=ezFrame(Type="logical", DefaultValue="TRUE", Description="generate plots for duplicates"),
                                    removeDuplicates=ezFrame(Type="logical", DefaultValue="TRUE", Description="remove Duplicates for CovCalculations"))
                }
              )
  )

ezMethodTeqc = function(input=NA, output=NA, param=NA){
  require(TEQC)
  library(chromstaR)
  require(GenomicAlignments)
  param[['build']] = unique(input$meta[['build']])
  setwdNew(basename(output$getColumn("Report")))
  if(basename(param$designFile) == param$designFile){
    param$designFile = list.files(file.path(TEQC_DESIGN_DIR, param$designFile), pattern='Covered\\.bed$', full.names = T)[1]
  }
  samples = input$getNames()
  jobList = input$getFullPaths("BAM")
  
  sGtfFile <- param$ezRef@refFeatureFile
  myGTF <- rtracklayer::import(sGtfFile)
  myGTF <- myGTF[mcols(myGTF)$type=='exon']
  myGTF <- myGTF[myGTF$gene_biotype=='protein_coding' & myGTF$source =='protein_coding',]
  
  keepCols = c('seqnames','start','end','strand','type','gene_id','gene_name')
  gtf_df = data.frame(myGTF,stringsAsFactors = FALSE)
  gtf_df = unique(gtf_df[,keepCols])
  ir <- IRanges(start = gtf_df$start, end = gtf_df$end)
  allExons <- RangedData(ranges = ir, space = gtf_df$seqnames, gene_id = gtf_df$gene_id, gene_name = gtf_df$gene_name, 
                         orientation = as.character(gtf_df$strand),typ = gtf_df$type)
  allExons <- as(allExons, "GRanges")
  
  #Create one Report per Sample:
  destDirs = ezMclapply(jobList, runTEQC, allExons, param, mc.cores = param$cores)
  
  #Create MultiSampleReport:
  destDirs = unlist(destDirs)
  TEQC::multiTEQCreport(singleReportDirs = destDirs,
                    samplenames = samples,
                    projectName = param$name,
                    targetsName = basename(dirname(param$designFile)),
                    referenceName = param[['refBuild']],
                    destDir = "multiTEQCreport",
                    k = c(1,5,10,20,30,50),
                    figureFormat = c("png"))
  
  capture.output(print(sessionInfo()), file = "sessionInfo.txt")
  htmlFile="00index.html"
  
  titles = list()
  titles[["TEQC-Report"]] = paste("TEQC-Report:", param$name)
  doc = openBsdocReport(title = titles[[length(titles)]])
  
  addDataset(doc, input$meta, param)
  
  titles[["MultiSample-Report"]] = "MultiSample-Report"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  addTxtLinksToReport(doc, "multiTEQCreport/index.html")
  titles[["Individual Reports"]] = "Individual Reports"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  addTxtLinksToReport(doc, paste0(destDirs, '/index.html'))
  
  
  titles[["Gene Coverage"]] = "GeneCoverage"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  minCov = 20
  coverageLinks = paste0(samples,'_coverage_allExons.txt')
  covData = list()
  minCovData = matrix(data = 0,nrow = 2,ncol = length(coverageLinks))
  rownames(minCovData) = c(paste0('>=',minCov,'x'),paste0('<',minCov,'x'))
  colnames(minCovData) = samples
  for(i in 1:length(coverageLinks)){
    covData[[i]] = ezRead.table(coverageLinks[i],row.names=NULL)
    avgCovPerGene = tapply(covData[[i]]$avgCoverage, INDEX=covData[[i]]$gene_id, mean)
    minCovData[,i] = c(sum(avgCovPerGene>=minCov),sum(avgCovPerGene<=minCov))
    #nReadPerGene = tapply(covData[[i]]$nReads, INDEX=covData[[i]]$gene_id, sum)
    #minExonCov = tapply(covData[[i]]$avgCoverage, INDEX=covData[[i]]$gene_id, min)
    #CovSDPerGene = tapply(covData[[i]]$avgCoverage, INDEX=covData[[i]]$gene_id, sd)
  }
  
  #### Simple Barplot
  addParagraph(doc,ezImageFileLink(plotCmd = expression({bp = barplot(minCovData,legend.text = T,
                                                                     names.arg = rep('',length(samples)),ylab='#Genes',
                                                                     main='Gene Coverage above minCov') 
                                                        text(x = bp, y = par("usr")[3] - 1, srt = 45,
                                                        adj = 1, labels = colnames(minCovData), xpd = TRUE)}), 
                                 file='genesAboveMinCov.png', 
                                 name="Genes above minCov",
                                 mouseOverText = "Genes above minCov"))
  addTxtLinksToReport(doc, coverageLinks)

  
  titles[["Misc"]] = "Misc"
  addTitle(doc, titles[[length(titles)]], 2, id=titles[[length(titles)]])
  closeBsdocReport(doc, htmlFile, titles)
  return("Success")
}

##' @title Runs the target enrichment quality control
##' @description Performs a target enrichment quality control and creates reports of the outcome.
##' @param param a list of parameters:
##' \itemize{
##'  \item{designFile}{ a file describing the regions selected by the enrichment kit.}
##'  \item{covUniformityPlot}{ a logical indicating whether to generate plots for coverage uniformity.}
##'  \item{covTargetLengthPlot}{ a logical indicating whether to generate plots for coverage vs. target length.}
##'  \item{duplicatesPlot}{ a logical indicating whether to generate plots for duplicates.}
##'  \item{paired}{ a logical indicating whether the samples are paired.}
##' }
##' @param file a character representing the path to the file containing the reads.
##' @template roxygen-template
runTEQC = function(file, allExons, param){
  sampleName = gsub('\\.bam', '', basename(file))
  destDir = paste0("report_", sampleName)
  targetsfile = param$designFile
  genomeSize = sum(as.numeric(system(paste("samtools","view -H",file,"|grep @SQ|cut -f 3|sed 's/LN://'"),intern = T)))
  
  
  reads <- readBamFileAsGRanges(file, chromosomes=NULL, pairedEndReads = param$paired, 
                                   max.fragment.width = 1000, min.mapq = 10, remove.duplicate.reads = param$removeDuplicates)
  
#  reads=TEQC::get.reads(file, filetype = "bam")
#  if(param$paired){
#    reads <- reads2pairs(reads)$readpairs
#    if(param$removeDuplicates){
#        ID.nondups <- names(unique(reads))
#        reads <- reads[names(reads) %in% ID.nondups,,drop = TRUE]
#    }
#  } else {
#      if(param$removeDuplicates){
#          reads <- unique(reads)  
#      } 
#  }
  skip = grep("^track", readLines(targetsfile, n=200))
  if (length(skip) == 0) skip = 0
  targets=TEQC::get.targets(targetsfile, skip=skip)
  #clean reads from mappings to unsupported chromosomes
  reads <- keepSeqlevels(reads, levels(seqnames(targets)), pruning.mode="coarse")
  seqlevels(reads) = as.character(unique(seqnames(reads)))
  strand(reads) = '*'
  
  TEQC::TEQCreport(sampleName=sampleName,
                   CovUniformityPlot = param$covUniformityPlot, CovTargetLengthPlot = param$covTargetLengthPlot, duplicatesPlot=param$duplicatesPlot,#CovGCPlot = T,
                   k = c(1,5,10,20,30,50),
                   targetsName=basename(dirname(targetsfile)),
                   referenceName=param[['refBuild']],
                   pairedend=FALSE,
                   destDir=destDir,
                   reads=reads,
                   targets=targets,
                   genomesize =genomeSize,figureFormat = c("png"), saveWorkspace = F)
  
  exonCoverage <- TEQC::coverage.target(reads, allExons, perBase = F, Offset = 0)$targetCoverages
  exonCoverage <- as.data.frame(TEQC::readsPerTarget(reads, exonCoverage))
  write.table(exonCoverage, file = paste0(sampleName,"_coverage_allExons.txt"), 
              sep = "\t", row.names = F, quote = F)
  
  return(destDir)
}
