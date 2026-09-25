###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodGatkDnaHaplotyper = function(input = NA, output = NA, param = NA) {
  bamFile = input$getFullPaths("BAM")
  ezSystem(paste("rsync -va", bamFile, "local.bam"))
  ezSystem(paste("rsync -va", paste0(bamFile, ".bai"), "local.bam.bai"))
  knownSites = list.files(
    param$ezRef["refVariantsDir"],
    pattern = 'vcf.gz$',
    full.names = T
  )
  dbsnpFile = knownSites[grep('dbsnp.*vcf.gz$', knownSites)]
  javaCall = paste0("java", " -Djava.io.tmpdir=. -Xmx", param$ram, "g")
  gatk = file.path(Sys.getenv("GATK"), 'gatk')

  genomeSeq = param$ezRef["refFastaFile"]
  sampleName = input$getNames()
  if (param$addReadGroup) {
    cmd = paste0(
      javaCall,
      " -jar ",
      Sys.getenv("Picard_jar"),
      " AddOrReplaceReadGroups",
      " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000",
      " I=local.bam",
      " O=withRg.bam SORT_ORDER=coordinate",
      " RGID=",
      sampleName,
      " RGPL=illumina RGSM=",
      sampleName,
      " RGLB=RGLB_",
      sampleName,
      " RGPU=RGPU_",
      sampleName,
      " VERBOSITY=WARNING"
    )
    ezSystem(cmd)
  } else {
    ezSystem('mv local.bam withRg.bam')
  }

  if (param$markDuplicates) {
    cmd = paste0(
      javaCall,
      " -jar ",
      Sys.getenv("Picard_jar"),
      " MarkDuplicates ",
      " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000",
      " I=withRg.bam",
      " O=dedup.bam",
      " REMOVE_DUPLICATES=false",
      " ASSUME_SORTED=true",
      " VALIDATION_STRINGENCY=SILENT",
      " METRICS_FILE=",
      "dupmetrics.txt",
      " VERBOSITY=WARNING"
    )
    ezSystem(cmd)
    ezSystem('mv dedup.bam withRg.bam')
  }

  ezSystem(paste("samtools", "index", "withRg.bam"))

  # if(param$splitNtrim){
  #    gatk = paste(javaCall, "-jar", Sys.getenv("GATK_jar"))
  #    cmd = paste(gatk, "-T SplitNCigarReads", "-R", genomeSeq,
  #                "-I", "withRg.bam",
  #                "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS",
  #                "-o splitNtrim.bam")
  #    ezSystem(cmd)
  #    ezSystem('mv splitNtrim.bam withRg.bam')
  #    ezSystem(paste("samtools", "index", "withRg.bam"))
  #  }

  if (param$targetFile != '') {
    param$targetFile = file.path(TARGET_ENRICHMENT_DESIGN_DIR, param$targetFile)
  }

  #BaseRecalibration is done only if known sites are available
  if (param$knownSitesAvailable) {
    baseRecalibration = paste(gatk, 'BaseRecalibrator')
    #knownSitesCMD = ''
    #for (j in 1:length(knownSites)){
    #  knownSitesCMD = paste(knownSitesCMD,paste("--knownSites", knownSites[j], collapse=','))
    #}

    cmd = paste(
      baseRecalibration,
      "-R",
      genomeSeq,
      "-I withRg.bam",
      "--known-sites",
      dbsnpFile,
      "--output recal.table"
    )

    if (param$targetFile != '') {
      cmd = paste(cmd, "-L", param$targetFile)
    }
    ezSystem(cmd)

    ApplyBQSR = paste(gatk, 'ApplyBQSR')
    cmd = paste(
      ApplyBQSR,
      "-R",
      genomeSeq,
      "-I withRg.bam",
      "--bqsr-recal-file recal.table",
      "--output recal.bam"
    )

    if (param$targetFile != '') {
      cmd = paste(cmd, "-L", param$targetFile)
    }
    ezSystem(cmd)
  } else {
    ezSystem('mv withRg.bam recal.bam')
  }
  ezSystem(paste("samtools", "index", "recal.bam"))
  ########### haplotyping
  haplotyperCall = paste(gatk, 'HaplotypeCaller')
  outputFile = paste0(sampleName, "-HC_calls.g.vcf.gz")
  cmd = paste(
    haplotyperCall,
    "-R",
    genomeSeq,
    "-I recal.bam",
    "-ERC GVCF",
    "--max-alternate-alleles 2",
    "-output",
    outputFile
  )

  if (param$knownSitesAvailable) {
    cmd = paste(cmd, "--dbsnp", dbsnpFile)
  }

  if (param$targetFile != '') {
    cmd = paste(cmd, "-L", param$targetFile)
  }
  if (param$getRealignedBam) {
    cmd = paste(cmd, "-bamout", paste0(sampleName, "-realigned.bam"))
  } else {
    cmd = paste(cmd, "--native-pair-hmm-threads", param$cores)
  }
  ezSystem(cmd)

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodGatkDnaHaplotyper(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppGatkDnaHaplotyper <-
  setRefClass(
    "EzAppGatkDnaHaplotyper",
    contains = "EzApp",
    methods = list(
      ## HaplotypeCaller defaults checked against GATK 4.6.1.0 (the Variants/GATK module default) --help.
      ## GATK HaplotypeCaller (GVCF) + samtools unconditional; BQSR (no separate reference) gated on
      ## knownSitesAvailable; Picard only when addReadGroup is true (AddOrReplaceReadGroups, R
      ## default FALSE) or markDuplicates is not false (MarkDuplicates, R default TRUE).
      citation = function(param = list()) {
        c(
          "McKenna, A. et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research 20, 1297-1303 (2010). https://doi.org/10.1101/gr.107524.110",
          "Poplin, R. et al. Scaling accurate genetic variant discovery to tens of thousands of samples. bioRxiv (2017) [preprint, not peer-reviewed]. https://doi.org/10.1101/201178",
          "DePristo, M.A. et al. A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics 43, 491-498 (2011). https://doi.org/10.1038/ng.806",
          if (isTRUE(as.logical(param$addReadGroup)) || !isFALSE(as.logical(param$markDuplicates))) "Picard Toolkit. Broad Institute. https://broadinstitute.github.io/picard/",
          "Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25(16), 2078-2079 (2009). https://doi.org/10.1093/bioinformatics/btp352"
        )
      },
      methods_facts = function(param = list()) {
        c(
          ## ezMethodGatkDnaHaplotyper (app-gatkDnaHaplotyper.R:23-45)
          if (isTRUE(as.logical(param$addReadGroup))) "Read groups were set with Picard AddOrReplaceReadGroups (ID and SM the sample name, LB RGLB_<sample>, PU RGPU_<sample>, PL illumina) and the BAM was coordinate-sorted.",
          ## app-gatkDnaHaplotyper.R:47-65
          if (isTRUE(as.logical(param$markDuplicates))) "Duplicate reads were flagged, not removed, with Picard MarkDuplicates (REMOVE_DUPLICATES=false); HaplotypeCaller then skips flagged duplicates through its default read filters.",
          ## app-gatkDnaHaplotyper.R:85-124, 139-141
          if (isTRUE(as.logical(param$knownSitesAvailable))) "Base quality scores were recalibrated with GATK BaseRecalibrator and ApplyBQSR using the dbSNP VCF of the reference build as the only known-sites set, and HaplotypeCaller took variant IDs from the same dbSNP file (--dbsnp).",
          if (isFALSE(as.logical(param$knownSitesAvailable))) "No base quality score recalibration was applied.",
          ## app-gatkDnaHaplotyper.R:80-82, 102-104, 117-119, 143-145
          if (ezIsSpecified(param$targetFile)) "Variant calling (and base recalibration, when run) was restricted to the intervals of targetFile (-L).",
          ## app-gatkDnaHaplotyper.R:126-151
          "Variants were called per sample with GATK HaplotypeCaller in GVCF mode (-ERC GVCF) with at most 2 alternate alleles per site (--max-alternate-alleles 2, GATK default 6); in GVCF mode HaplotypeCaller sets the calling confidence to 0 (the default 30 applies at joint genotyping), and minimum base quality 10 and minimum mapping quality 20 were GATK defaults.",
          "The output is an unfiltered per-sample gVCF intended for joint genotyping; no variant filtering was applied in this step."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodGatkDnaHaplotyper
        name <<- "EzAppGatkDnaHaplotyper"
        appDefaults <<- rbind(
          addReadGroup = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "add ReadGroup to BAM"
          ),
          getRealignedBam = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "for IGV check"
          ),
          splitNtrim = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "for RNA-Seq data"
          ),
          targetFile = ezFrame(
            Type = "character",
            DefaultValue = "",
            Description = "restrict to targeted genomic regions"
          ),
          markDuplicates = ezFrame(
            Type = "logical",
            DefaultValue = TRUE,
            Description = "not recommended for gene panels, exomes"
          )
        )
      }
    )
  )
