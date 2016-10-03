#TODO: input bamFile -> output vcf with gvcf support
#Add readGroup, remove duplicates with picard yes/no, run haplotyper full or with intervall for exome kit 
#If neccessary: resort and index genome.fa file locally, use GATK 3.6

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodGatkDnaHaplotyper = function(input=NA, output=NA, param=NA){
  knownSites = list.files(param$ezRef["refVariantsDir"],pattern='vcf$',full.names = T) 
  GATK_JAR='/usr/local/ngseq/src/GenomeAnalysisTK_3.6/GenomeAnalysisTK.jar'
  javaCall = paste0(JAVA, " -Djava.io.tmpdir=. -Xmx", param$ram, "g")
  bamFile = '/srv/GT/analysis/lopitz/GATK-tutorial_data/bams/NA12877_wgs_20.bam'#input$getFullPaths("BAM")
  #localBamFile = .getBamLocally(bamFile)
  genomeSeq = '/srv/GT/analysis/lopitz/GATK-tutorial_data/ref/human_g1k_b37_20.fasta' #param$ezRef["refFastaFile"]
  sampleName = 'NA12877_wgs_20' #names(bamFile)
  cmd = paste0(javaCall, " -jar ", PICARD_JAR, " AddOrReplaceReadGroups",
               " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", bamFile,
               " O=withRg.bam SORT_ORDER=coordinate",
               " RGID=RGID_", sampleName, " RGPL=illumina RGSM=", sampleName, " RGLB=RGLB_", sampleName, " RGPU=RGPU_", sampleName,
               " VERBOSITY=WARNING")
  ezSystem(cmd)
  
  if(param$markDuplicates){
    cmd = paste0(javaCall, " -jar ", PICARD_JAR, " MarkDuplicates ",
                 " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", "withRg.bam",
                 " O=", "dedup.bam",
                 " REMOVE_DUPLICATES=false",
                 " ASSUME_SORTED=true",
                 " VALIDATION_STRINGENCY=SILENT",
                 " METRICS_FILE=" ,"dupmetrics.txt",
                 " VERBOSITY=WARNING")
    ezSystem(cmd)
    ezSystem('mv dedup.bam withRg.bam')
  }
  
  ezSystem(paste(SAMTOOLS, "index", "withRg.bam"))
  #BaseRecalibration
  baseRecalibration1 = paste(javaCall,"-jar", GATK_JAR, " -T BaseRecalibrator")
  knownSitesCMD = ''
  for (j in 1:length(knownSites)){
    knownSitesCMD = paste(knownSitesCMD,paste("--knownSites", knownSites[j], collapse=','))
  }
  
  cmd = paste(baseRecalibration1, "-R", genomeSeq,
              "-I withRg.bam",
              knownSitesCMD,
              "--out recal.table", 
              "-nct", param$cores)
  
  if(!is.null(param$targetFile)){
    cmd = paste(cmd,
                "-L", param$targetFile)
  }
  ezSystem(cmd)
  
  baseRecalibration2 = paste(javaCall,"-jar", GATK_JAR, " -T PrintReads")
  cmd = paste(baseRecalibration2, "-R", genomeSeq,
              "-I withRg.bam",
              "-BQSR recal.table",
              "-o recal.bam",
              "-nct", param$cores)
  
  if(!is.null(param$targetFile)){
    cmd = paste(cmd,
                "-L", param$targetFile)
  }
  ezSystem(cmd)

  ########### haplotyping
  haplotyperCall = paste(javaCall,"-jar", GATK_JAR, " -T HaplotypeCaller")
  cmd = paste(haplotyperCall, "-R", genomeSeq,
              "-I recal.bam",
              "-ERC GVCF",
              "-variant_index_type LINEAR -variant_index_parameter 128000",
              "-o", paste0(sampleName, "-HC_calls.vcf"))
  
  if(!is.null(param$targetFile)){
    cmd = paste(cmd,
                "-L", param$targetFile)
  }
  if(param$getRealignedBam){
    cmd = paste(cmd,
              "-bamout", paste0(sampleName, "-realigned.bam"),
              "-forceActive",
              "-disableOptimizations", 
              "-ip 100")
  } else {
    cmd = paste(cmd,
                "-nct", param$cores)
  }
  ezSystem(cmd)
  
  #TODO: recalibrate Variants: 1.Snps, 2. InDels
  #java –jar GenomeAnalysisTK.jar –T VariantRecalibrator \
  #–R human.fasta \
  #–input raw.SNPs.vcf \
  #–resource: {see next slide} \
  #–an DP –an QD –an FS –an MQRankSum {...} \
  #–mode SNP \
  #–recalFile raw.SNPs.recal \
  #–tranchesFile raw.SNPs.tranches \
  #–rscriptFile recal.plots.R
  
  #java –jar GenomeAnalysisTK.jar –T ApplyRecalibraHon \
  #–R human.fasta \
  #–input raw.vcf \
  #–mode SNP \
  #–recalFile raw.SNPs.recal \
  #–tranchesFile raw.SNPs.tranches \
  #–o recal.SNPs.vcf \
  #–ts_filter_level 99.0
  
  
  ## filter the vcf file
  gc()
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodGatkDnaHaplotyper(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppGatkDnaHaplotyper <-
  setRefClass("EzAppGatkDnaHaplotyper",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodGatkDnaHaplotyper
                  name <<- "EzAppGatkDnaHaplotyper"
                  appDefaults <<- rbind(getRealignedBam = ezFrame(Type="logical",  DefaultValue=FALSE, Description="for IGV check"),
                                        targetFile = ezFrame(Type="character",  DefaultValue="", Description="restrict to targeted genomic regions"),
                                        markDuplicates = ezFrame(Type="logical",  DefaultValue=TRUE, Description="not recommended for gene panels, exomes"), 
                                        vcfFilt.minAltCount = ezFrame(Type="integer",  DefaultValue=8,  Description="minimum coverage for the alternative variant"),
                                        vcfFilt.minReadDepth = ezFrame(Type="integer",  DefaultValue=20,  Description="minimum read depth"))
                }
              )
  )
