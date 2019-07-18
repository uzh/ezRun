###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodGatkDnaHaplotyper = function(input=NA, output=NA, param=NA){
  bamFile = input$getFullPaths("BAM")
  ezSystem(paste("rsync -va", bamFile, "local.bam"))
  ezSystem(paste("rsync -va", paste0(bamFile, ".bai"), "local.bam.bai"))
  knownSites = list.files(param$ezRef["refVariantsDir"],pattern='vcf.gz$',full.names = T)
  dbsnpFile = knownSites[grep('dbsnp.*vcf.gz$', knownSites)]
  javaCall = paste0("java", " -Djava.io.tmpdir=. -Xmx", param$ram, "g")
  
  genomeSeq = param$ezRef["refFastaFile"]
  sampleName = names(bamFile)
  if(param$addReadGroup){
    cmd = paste0(javaCall, " -jar ", Sys.getenv("Picard_jar"), " AddOrReplaceReadGroups",
               " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=local.bam",
               " O=withRg.bam SORT_ORDER=coordinate",
               " RGID=RGID_", sampleName, " RGPL=illumina RGSM=", sampleName, " RGLB=RGLB_", sampleName, " RGPU=RGPU_", sampleName,
               " VERBOSITY=WARNING")
    ezSystem(cmd) } else {
    ezSystem('mv local.bam withRg.bam')
  }
  
  if(param$markDuplicates){
    cmd = paste0(javaCall, " -jar ", Sys.getenv("Picard_jar"), " MarkDuplicates ",
                 " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=withRg.bam",
                 " O=dedup.bam",
                 " REMOVE_DUPLICATES=false",
                 " ASSUME_SORTED=true",
                 " VALIDATION_STRINGENCY=SILENT",
                 " METRICS_FILE=" ,"dupmetrics.txt",
                 " VERBOSITY=WARNING")
    ezSystem(cmd)
    ezSystem('mv dedup.bam withRg.bam')
  }
  
  ezSystem(paste("samtools", "index", "withRg.bam"))
  
  if(param$splitNtrim){
    gatk = paste(javaCall, "-jar", Sys.getenv("GATK_jar"))
    cmd = paste(gatk, "-T SplitNCigarReads", "-R", genomeSeq,
                "-I", "withRg.bam",
                "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS",
                "-o splitNtrim.bam") 
    ezSystem(cmd)
    ezSystem('mv splitNtrim.bam withRg.bam')
    ezSystem(paste("samtools", "index", "withRg.bam"))
  }
  
  if(param$targetFile != ''){
    param$targetFile = file.path(TARGET_ENRICHMENT_DESIGN_DIR, param$targetFile)
  }
  
  #BaseRecalibration is done only if known sites are available
  if(param$knownSitesAvailable){
  baseRecalibration1 = paste(javaCall,"-jar", Sys.getenv("GATK_jar"), " -T BaseRecalibrator")
  #knownSitesCMD = ''
  #for (j in 1:length(knownSites)){
  #  knownSitesCMD = paste(knownSitesCMD,paste("--knownSites", knownSites[j], collapse=','))
  #}
  
  cmd = paste(baseRecalibration1, "-R", genomeSeq,
              "-I withRg.bam",
              "--knownSites", dbsnpFile,
              "--out recal.table", 
              "-nct", param$cores)
  
  if(param$targetFile != ''){
    cmd = paste(cmd,
                "-L", param$targetFile)
  }
  ezSystem(cmd)
  
  
  baseRecalibration2 = paste(javaCall,"-jar", Sys.getenv("GATK_jar"), " -T PrintReads")
  cmd = paste(baseRecalibration2, "-R", genomeSeq,
              "-I withRg.bam",
              "-BQSR recal.table",
              "-o recal.bam",
              "-nct", param$cores)
  
  if(param$targetFile != ''){
    cmd = paste(cmd,
                "-L", param$targetFile)
  }
  ezSystem(cmd) } else {
    ezSystem('mv withRg.bam recal.bam')
  }
  ezSystem(paste("samtools", "index", "recal.bam"))
  ########### haplotyping
  haplotyperCall = paste(javaCall,"-jar", Sys.getenv("GATK_jar"), " -T HaplotypeCaller")
  outputFile = paste0(sampleName, "-HC_calls.g.vcf")
  cmd = paste(haplotyperCall, "-R", genomeSeq,
              "-I recal.bam",
              "--emitRefConfidence GVCF",
              "--max_alternate_alleles 2",
              "-o", outputFile)
  
  if(param$knownSitesAvailable){
    cmd = paste(cmd,
                "--dbsnp", dbsnpFile)
  }
  
  if(param$targetFile != ''){
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
  ezSystem(paste("bgzip","-c",outputFile, ">",paste0(outputFile,".gz")))
  ezSystem(paste("tabix","-p vcf",paste0(outputFile,".gz")))
  
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
                  appDefaults <<- rbind(addReadGroup = ezFrame(Type="logical",  DefaultValue=FALSE, Description="add ReadGroup to BAM"),
                                        getRealignedBam = ezFrame(Type="logical",  DefaultValue=FALSE, Description="for IGV check"),
                                        splitNtrim = ezFrame(Type="logical",  DefaultValue=FALSE, Description="for RNA-Seq data"),
                                        targetFile = ezFrame(Type="character",  DefaultValue="", Description="restrict to targeted genomic regions"),
                                        markDuplicates = ezFrame(Type="logical",  DefaultValue=TRUE, Description="not recommended for gene panels, exomes"))
                }
              )
  )
