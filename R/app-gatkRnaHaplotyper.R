###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodGatkRnaHaplotyper = function(input=NA, output=NA, param=NA){
  library(rtracklayer)
  standardCallConfidence <- 20
  bamFile <- getBamLocally(input$getFullPaths("BAM"))
  #knownSites = list.files(param$ezRef["refVariantsDir"],pattern='vcf.gz$',full.names = T)
  #dbsnpFile = knownSites[grep('dbsnp.*vcf.gz$', knownSites)]
  dbsnpFile = param$dbsnpFile
  stopifnot(file.exists(dbsnpFile))
  stopifnot(param$markDuplicates == FALSE)
  #javaCall = paste0("java", " -Djava.io.tmpdir=. -Xmx", param$ram, "g")
  gatkCall = paste0('gatk --java-options "-Djava.io.tmpdir=. -Xmx', param$ram, 'g -Xms4g" ')
  
  ## define the exome intervals
  exomeInterVals <- "exome.bed"
  gtf <- import(param$ezRef["refFeatureFile"])
  ## TODO: we only work on the main chromosomes and use only lncRNA and protein coding
  gtf <- gtf[grepl("chr", seqnames(gtf)) & gtf$type =="exon" & gtf$gene_biotype %in% c("lncRNA", "protein_coding")]
  write.table(data.frame(seqnames(gtf), start(gtf)-1, end(gtf)), 
               file=exomeInterVals, quote = F, sep="\t", col.names = F, row.names = F)
  
  
  
  genomeSeq = param$ezRef["refFastaFile"]
  sampleName = names(bamFile)
  if(param$addReadGroup){
    cmd = paste(gatkCall, " AddOrReplaceReadGroups",
                 "-I", bamFile,
                 "-O", "withRg.bam",
                 "-LB", sampleName,
                 "-PL", "ILLUMINA",
                 "-PU", sampleName,
                 "-SM", sampleName)
    ezSystem(cmd)
  } else {
    ezSystem('mv local.bam withRg.bam')
  }

  cmd = paste(gatkCall, "SplitNCigarReads", 
              "-R", genomeSeq,
              "-I", "withRg.bam",
              "-L", exomeInterVals,
              ## this is the default! "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS",
              "-o splitNtrim.bam")
  ezSystem(cmd)

  #BaseRecalibration is done only if known sites are available
  if(param$knownSitesAvailable){
    cmd = paste(gatkCall, "BaseRecalibrator", 
                "-R", genomeSeq,
                "-I splitNtrim.bam",
                "-L", exomeInterVals,
                "--known-sites", dbsnpFile,
                "--use-original-qualities", ## as specified in the best practices
                "--output recal.table")
    ezSystem(cmd)
    
    cmd = paste(gatkCall,'ApplyBQSR',
                "-R", genomeSeq,
                "-I splitNtrim.bam",
                "-L", exomeInterVals,
                "--add-output-sam-program-record",
                "--bqsr-recal-file recal.table",
                "--output recal.bam")
  } else {
    ezSystem('mv splitNtrim.bam recal.bam')
  }
  
  ########### haplotyping
  outputFile = basename(output$getColumn("GVCF"))#  paste0(sampleName, ".g.vcf.gz")
  cmd = paste(gatkCal,'HaplotypeCaller',
              "-R", genomeSeq,
              "-I recal.bam",
              "-L", exomeInterVals,
              "-ERC GVCF",
              "-dont-use-soft-clipped-bases",
              # TODO: conflicts with -ERC GVCF "--standard-min-confidence-threshold-for-calling", standardCallConfidence,
              "--dbsnp", dbsnpFile,
              "-O", outputFile)
  ezSystem(cmd)
  
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodGatkDnaHaplotyper(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppGatkRnaHaplotyper <-
  setRefClass("EzAppGatkRnaHaplotyper",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodGatkRnaHaplotyper
                  name <<- "EzAppGatkRnaHaplotyper"
                  appDefaults <<- rbind(addReadGroup = ezFrame(Type="logical",  DefaultValue=FALSE, Description="add ReadGroup to BAM"),
                                        #getRealignedBam = ezFrame(Type="logical",  DefaultValue=FALSE, Description="for IGV check"),
                                        #splitNtrim = ezFrame(Type="logical",  DefaultValue=FALSE, Description="for RNA-Seq data"),
                                        #targetFile = ezFrame(Type="character",  DefaultValue="", Description="restrict to targeted genomic regions"),
                                        markDuplicates = ezFrame(Type="logical",  DefaultValue=FALSE, Description="not recommended for gene panels, exomes"))
                }
              )
  )
