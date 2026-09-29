###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodGatkRnaHaplotyper = function(input = NA, output = NA, param = NA) {
  library(rtracklayer)
  standardCallConfidence <- 20
  bamFile <- getBamLocally(input$getFullPaths("BAM"))
  #knownSites = list.files(param$ezRef["refVariantsDir"],pattern='vcf.gz$',full.names = T)
  #dbsnpFile = knownSites[grep('dbsnp.*vcf.gz$', knownSites)]
  dbsnpFile = param$dbsnpFile
  stopifnot(param$markDuplicates == FALSE)
  #javaCall = paste0("java", " -Djava.io.tmpdir=. -Xmx", param$ram, "g")
  gatkCall = paste0(
    'gatk --java-options "-Djava.io.tmpdir=. -Xmx',
    param$ram,
    'g -Xms4g" '
  )

  ## define the exome intervals
  exomeInterVals <- "exome.bed"
  gtf <- import(param$ezRef["refFeatureFile"])
  ## TODO: we only work on the main chromosomes and use only lncRNA and protein coding
  if (any(grepl("chr", seqnames(gtf)))) {
    gtf <- gtf[grepl("chr", seqnames(gtf)), ]
  }
  gtf <- gtf[
    gtf$type == "exon" & gtf$gene_biotype %in% c("lncRNA", "protein_coding")
  ]
  write.table(
    data.frame(seqnames(gtf), start(gtf) - 1, end(gtf)),
    file = exomeInterVals,
    quote = F,
    sep = "\t",
    col.names = F,
    row.names = F
  )

  genomeSeq = param$ezRef["refFastaFile"]
  sampleName = names(bamFile)
  if (param$addReadGroup) {
    cmd = paste(
      gatkCall,
      " AddOrReplaceReadGroups",
      "-I",
      bamFile,
      "-O",
      "withRg.bam",
      "-LB",
      sampleName,
      "-PL",
      "ILLUMINA",
      "-PU",
      sampleName,
      "-SM",
      sampleName
    )
    ezSystem(cmd)
  } else {
    ezSystem('mv local.bam withRg.bam')
  }

  ezSystem("samtools index withRg.bam")
  cmd = paste(
    gatkCall,
    "SplitNCigarReads",
    "-R",
    genomeSeq,
    "-I",
    "withRg.bam",
    "-L",
    exomeInterVals,
    ## this is the default! "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS",
    "-O splitNtrim.bam"
  )
  ezSystem(cmd)
  ezSystem("samtools index splitNtrim.bam")

  #BaseRecalibration is done only if known sites are available
  if (ezIsSpecified(dbsnpFile) && file.exists(dbsnpFile)) {
    cmd = paste(
      gatkCall,
      "BaseRecalibrator",
      "-R",
      genomeSeq,
      "-I splitNtrim.bam",
      "-L",
      exomeInterVals,
      "--known-sites",
      dbsnpFile,
      "--use-original-qualities", ## as specified in the best practices
      "--output recal.table"
    )
    ezSystem(cmd)

    cmd = paste(
      gatkCall,
      'ApplyBQSR',
      "-R",
      genomeSeq,
      "-I splitNtrim.bam",
      "-L",
      exomeInterVals,
      "--add-output-sam-program-record",
      "--bqsr-recal-file recal.table",
      "--output recal.bam"
    )
    ezSystem(cmd)
  } else {
    ezSystem('mv splitNtrim.bam recal.bam')
  }
  ezSystem("samtools index recal.bam")

  ########### haplotyping
  outputFile = basename(output$getColumn("GVCF")) #  paste0(sampleName, ".g.vcf.gz")
  cmd = paste(
    gatkCall,
    'HaplotypeCaller',
    "-R",
    genomeSeq,
    "-I recal.bam",
    "-L",
    exomeInterVals,
    "-ERC GVCF",
    "-dont-use-soft-clipped-bases",
    # TODO: conflicts with -ERC GVCF "--standard-min-confidence-threshold-for-calling", standardCallConfidence,
    ifelse(
      ezIsSpecified(dbsnpFile) && file.exists(dbsnpFile),
      paste("--dbsnp", dbsnpFile),
      ""
    ),
    "-O",
    outputFile
  )
  ezSystem(cmd)

  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodGatkDnaHaplotyper(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppGatkRnaHaplotyper <-
  setRefClass(
    "EzAppGatkRnaHaplotyper",
    contains = "EzApp",
    methods = list(
      methods_facts = function(param = list()) {
        known <- length(param) > 0
        dbsnp <- ezIsSpecified(param$dbsnpFile)
        c(
          ## ezMethodGatkRnaHaplotyper: exome.bed built from the GTF, passed with -L to every GATK step after read groups
          paste0("The interval file exome.bed was built inside the job from the exon features of lncRNA and protein_coding genes (gene_biotype) of the reference GTF, keeping only sequences whose name contains 'chr' when any does; no exome capture or target file was used. It was passed with -L to SplitNCigarReads and HaplotypeCaller", if (!(known && !dbsnp)) " and, when it ran, to base recalibration", "."),
          ## ezMethodGatkRnaHaplotyper: AddOrReplaceReadGroups -LB/-PL/-PU/-SM, no -ID
          if (isTRUE(as.logical(param$addReadGroup))) "Read groups were replaced with GATK AddOrReplaceReadGroups, setting library (-LB), platform unit (-PU) and sample (-SM) to the sample name and platform (-PL) to ILLUMINA; no read-group ID (-ID) was passed.",
          ## ezMethodGatkRnaHaplotyper: stopifnot(param$markDuplicates == FALSE)
          "Duplicate reads were not marked: the app stops when markDuplicates is set.",
          ## ezMethodGatkRnaHaplotyper: SplitNCigarReads with -R, -I, -L only
          "Reads spanning introns were split at N CIGAR operations with GATK SplitNCigarReads, given only the reference, the BAM and the intervals.",
          ## ezMethodGatkRnaHaplotyper: BQSR only if (ezIsSpecified(dbsnpFile) && file.exists(dbsnpFile))
          if (!known) "Base quality score recalibration runs only when dbsnpFile names an existing file; otherwise the SplitNCigarReads BAM is only renamed recal.bam.",
          if (known && dbsnp) paste0("When the file ", param$dbsnpFile, " existed, base quality scores were recalibrated with GATK BaseRecalibrator (--known-sites ", param$dbsnpFile, ", --use-original-qualities) and ApplyBQSR, and HaplotypeCaller took variant IDs from it (--dbsnp); without the file no recalibration was done."),
          if (known && !dbsnp) "No base quality score recalibration was done because dbsnpFile was empty: the SplitNCigarReads BAM was only renamed recal.bam, and HaplotypeCaller got no --dbsnp file.",
          ## ezMethodGatkRnaHaplotyper: HaplotypeCaller -ERC GVCF -dont-use-soft-clipped-bases; the confidence option is commented out
          "Variants were called per sample with GATK HaplotypeCaller in GVCF mode (-ERC GVCF), ignoring soft-clipped bases (-dont-use-soft-clipped-bases); no calling-confidence threshold was passed.",
          "The output is one unfiltered gVCF per sample for joint genotyping; this app applies no variant filtering."
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodGatkRnaHaplotyper
        name <<- "EzAppGatkRnaHaplotyper"
        appDefaults <<- rbind(
          addReadGroup = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "add ReadGroup to BAM"
          ),
          #getRealignedBam = ezFrame(Type="logical",  DefaultValue=FALSE, Description="for IGV check"),
          #splitNtrim = ezFrame(Type="logical",  DefaultValue=FALSE, Description="for RNA-Seq data"),
          #targetFile = ezFrame(Type="character",  DefaultValue="", Description="restrict to targeted genomic regions"),
          markDuplicates = ezFrame(
            Type = "logical",
            DefaultValue = FALSE,
            Description = "not recommended for gene panels, exomes"
          )
        )
      }
    )
  )
