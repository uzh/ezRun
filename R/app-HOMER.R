###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppHomerDiffPeaks <-
  setRefClass("EzAppHomerDiffPeaks",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodHomerDiffPeaks
                  name <<- "EzAppHomerDiffPeaks"
                  appDefaults <<- rbind(refBuildHOMER=ezFrame(Type="character", DefaultValue="hg38", Description="Genome version to use from HOMER: hg38, mm10, danRer10, etc."),
                                        repFoldChange=ezFrame(Type="numeric", DefaultValue=2, Description="Replicate fold change cutoff for peak identification (calculated by DESeq2)"),
                                        repFDR=ezFrame(Type="numeric", DefaultValue=0.05, Description="Replicate FDR cutoff for peak identification (calculated by DESeq2)"),
                                        balanced=ezFrame(Type="logical", DefaultValue=TRUE, Description="Do not force the use of normalization factors to match total mapped reads.  This can be useful when analyzing differential peaks between similar data (for example H3K27ac) where we expect similar levels in all experiments. Applying this allows the data to essentially be quantile normalized during the differential calculation."),
                                        style=ezFrame(Type="character", DefaultValue="histone", Description="Style of peaks found by findPeaks during features selection (factor, histone, super, groseq, tss, dnase, mC)")
                  )
                }
              )
  )

ezMethodHomerDiffPeaks = function(input=NA, output=NA, param=NA, 
                                  htmlFile="00index.html"){
  cwd <- getwd()
  setwdNew(basename(output$getColumn("Report")))
  on.exit(setwd(cwd))
  
  stopifnot(param$sampleGroup != param$refGroup)
  
  ## I don't have better names in mind for now. Just use first and second
  firstSamples <- input$getNames()[input$getColumn(param$grouping) %in% 
                                     param$sampleGroup]
  secondSamples <- input$getNames()[input$getColumn(param$grouping) %in% 
                                      param$refGroup]
  
  useSamples = c(firstSamples, secondSamples)
  input <- input$subset(useSamples)
  
  bamFiles <- input$getFullPaths("BAM")
  localBamFiles <- sapply(bamFiles, getBamLocally)
  localSamFiles <- sub('.bam$', '.sam', localBamFiles)

    for (i in 1:length(localBamFiles)){
      cmd <- paste('samtools view -h', localBamFiles[i], '>', localSamFiles[i])
     ezSystem(cmd)
  }
  
  mcmapply(makeTagDirectory, inBam=localSamFiles, 
           outputDir=names(localBamFiles),
           MoreArgs=list(genome=param$refBuildHOMER),
           mc.cores=param$cores)
  
  if(all(localBamFiles != bamFiles))
    file.remove(c(localBamFiles, paste0(localBamFiles, ".bai")))
  
  if(length(firstSamples) >= 2L || length(secondSamples) >= 2L){
    ## The experiments with replicates
    cmd <- paste("getDifferentialPeaksReplicates.pl -DESeq2", 
                 "-genome", param$refBuildHOMER, 
                 "-f", param$repFoldChange,
                 "-q", param$repFDR,
                 ifelse(param$balanced, "-balanced", ""),
                 "-style", param$style)
    cmd <- paste(cmd, "-t", paste(firstSamples, collapse=" "),
                 "-b", paste(secondSamples, collapse=" "),
                 ">", basename(output$getColumn("DiffPeak")))
    ezSystem(cmd)
  }else{
    ## The experiments without replicates;
    ## focus on tss regions
    #cmd <- paste("annotatePeaks.pl tss", param$ezRef["refFastaFile"], 
    #             "-gtf", param$ezRef["refFeatureFile"], "> tss.txt")
    cmd <- paste("annotatePeaks.pl tss", param$refBuildHOMER,
                 "> tss.txt")
    ezSystem(cmd)
    cmd <- paste("getDifferentialPeaks", "tss.txt",
                 firstSamples, secondSamples,
                 ">", basename(output$getColumn("DiffPeak")))
    ezSystem(cmd)
    file.remove("tss.txt")
  }
  
  file.remove(localSamFiles)
  unlink(names(localBamFiles), recursive=TRUE) ## clean the tag directory
  
  return("Success")
}

makeTagDirectory <- function(inBam, outputDir, genome=NULL, checkGC=FALSE,
                             isAntisense=FALSE, strandedPaired=FALSE){
  setEnvironments("HOMER")
  cmd <- paste("makeTagDirectory", outputDir, paste(inBam, collapse=" "),
               "-format sam")
  if(!is.null(genome)){
    cmd <- paste(cmd, "-genome", genome)
  }
  if(isTRUE(checkGC)){
    cmd <- paste(cmd, "-checkGC")
  }
  if(isTRUE(isAntisense))
    cmd <- paste(cmd, "-flip")
  if(isTRUE(strandedPaired))
    cmd <- paste(cmd, "-sspe")
  ezSystem(cmd)
}
