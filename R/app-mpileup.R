###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMpileup = function(input=NA, output=NA, param=NA){

  require("VariantAnnotation")
  
  reportDir = basename(output$getColumn("Report"))
  htmlFile = basename(output$getColumn("Html"))
  vcfOutputFile = output$getColumn("VCF")
  
  bamFiles = input$getFullPaths("BAM")
  bamDataset = input$meta
  genomeSeq = param$ezRef["refFastaFile"]
  nBamsInParallel = min(4, param$cores)
  bamFilesClean = ezMclapply(names(bamFiles), function(sampleName){
    javaCall = paste0("java", " -Djava.io.tmpdir=. -Xmx", floor(param$ram/nBamsInParallel), "g")
    setwdNew(paste(sampleName, "proc", sep="-"))
    bf = bamFiles[sampleName]
    obf = file.path(getwd(), basename(bf))
    if (!is.null(param$region) && param$region != ""){
      bamParam = ScanBamParam(what=scanBamWhat())
      reg = splitRegion(param$region)
      if (is.na(reg$end)){
        seqLengths = ezBamSeqLengths(bf)
        reg$start = 1
        reg$end = seqLengths[reg$seq]
      }
      bamWhich(bamParam) = GRanges(seqnames=reg$seq, ranges=IRanges(start=reg$start, end=reg$end))
      filterBam(bf, "local.bam", param=bamParam)
      ezSystem(paste("samtools", "index", "local.bam"))
    } else {
      ezSystem(paste("cp", bf, "local.bam"))
      ezSystem(paste("cp", paste0(bf, ".bai"), "local.bam.bai"))
    }
    cmd = paste0(javaCall, " -jar ", Sys.getenv("Picard_jar"), " AddOrReplaceReadGroups",
                " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", "local.bam",
                " O=withRg.bam SORT_ORDER=coordinate",
                " RGID=RGID_", sampleName, " RGPL=illumina RGSM=", sampleName, " RGLB=RGLB_", sampleName, " RGPU=RGPU_", sampleName,
                " VERBOSITY=WARNING",
                " > addreplace.stdout 2> addreplace.stderr")
    ezSystem(cmd)
    file.remove("local.bam")
    
    cmd = paste0(javaCall, " -jar ", Sys.getenv("Picard_jar"), " MarkDuplicates ",
                " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", "withRg.bam",
                " O=", "dedup.bam",
                " REMOVE_DUPLICATES=false", ## do not remove, do only mark
                " ASSUME_SORTED=true",
                " VALIDATION_STRINGENCY=SILENT",
                " METRICS_FILE=" ,"dupmetrics.txt",
                " VERBOSITY=WARNING",
                " >markdup.stdout 2> markdup.stderr")
    ezSystem(cmd)
    file.remove("withRg.bam")
    
    #     ezSystem(paste("samtools", "index", "dedup.bam"))
    #     gatk = paste(javaCall, "-jar", "$GATK_jar")
    #     cmd = paste(gatk, "-T SplitNCigarReads", "-R", genomeSeq,
    #                 "-I", "dedup.bam",
    #                 "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS",
    #                 "-o", obf,
    #                 "> splitncigars.stdout 2> splitncigars.stderr") 
    #     ezSystem(cmd)
    #     file.remove("dedup.bam")
    ezSystem(paste("mv", "dedup.bam", obf))
    ezSystem(paste("samtools", "index", obf))
    setwd("..")
    return(obf)
  }, mc.cores=nBamsInParallel, mc.preschedule=FALSE)
  
  mpileupCmd = paste("bcftools", "mpileup","-Ou",
                     "-f", param$ezRef["refFastaFile"],
                     param$mpileupOptions,
                     ifelse(param$region == "", "", paste("--region", param$region)),
                     paste(bamFilesClean, collapse=" "))
  callCmd = paste("bcftools", "call",
                  "-Ou",
                  param$callOptions,
                  "-") ## read from stdin
  filterCmd = paste("bcftools", "filter",
                    "--output-type z",
                    "--output", basename(vcfOutputFile),
                    param$filterOptions,
                    "-") ## read from stdin
  ezSystem(paste(mpileupCmd, "|", callCmd, "|", filterCmd))
  indexTabix(basename(vcfOutputFile),format = "vcf")
  
  ## write an igv link
  # writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
  #                 vcfUrls = paste(PROJECT_BASE_URL, vcfOutputFile, sep="/") )
  # writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", vcfOutputFile),
  #              sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
  chromSizes = ezChromSizesFromVcf(basename(vcfOutputFile))
  genotype = geno(readVcf( basename(vcfOutputFile), genome="genomeDummy"))
  gt = genotype$GT
  gt[genotype$DP < param$minReadDepth] = "lowCov" ## those calls will become NA in subsequent analyses
  
  setwdNew(reportDir)
  makeRmdReport(
    input=input,
    output = output, param = param, chromSizes=chromSizes, gt=gt,
    rmdFile = "Mpileup.Rmd", reportTitle = param$comparison
  )
  setwd("..")
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodMpileup(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppMpileup <-
  setRefClass("EzAppMpileup",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMpileup
                  name <<- "EzAppMpileup"
                  appDefaults <<- rbind(region=ezFrame(Type="character",  DefaultValue="",  Description="should analysis be done on a region of the genom only chr:start-end"),
                                    mpileupOptions=ezFrame(Type="character",  DefaultValue="",	Description="options to mpileup command"),
                                    callOptions=ezFrame(Type="character", DefaultValue="", Description="options to call command"),
                                    filterOptions=ezFrame(Type="character", DefaultValue="", Description="options to filter command"),
                                    minReadDepth=ezFrame(Type="integer", DefaultValue="10", Description="use for clustering only SNV with coverage higher than"))
                }
              )
  )
