###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


.runGatk = function(input=NA, output=NA, param=NA){
  
  param = fillWithDefaults(param) ## TODOMF: function doesn't exist
  options(cores=param$cores)
  ref = param$ezRef["refFastaFile"]
  gatk = paste("java -Xmx8g -Djava.io.tmpdir=.", "-jar", Sys.getenv("GATK_jar"))
  inputBam = ezFullPaths(param$dataRoot, input$BAM) ## TODO: Refactor when updating function
  localBam = basename(input$BAM)
  realigned = sub(".bam$", ".realigned.bam", localBam)
  stopifnot(realigned != localBam)
  intervals = basename(output$"Intervals [File]")
  vcf = basename(output$"VCF [File]")
  vcfReport = basename(output$"VCF Report [File]")
  
  
  cmd = paste0("java -Xmx8g -Djava.io.tmpdir=. -jar ", Sys.getenv("Picard_jar"), 
               " AddOrReplaceReadGroups ",
              " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", inputBam,
              " O=", "tmp.bam", ' SORT_ORDER=coordinate',
              " RGID=RGID_", input$Name, " RGPL=illumina RGSM=RGSM_", input$Name, " RGLB=RGLB_", input$Name, " RGPU=RGPU_", input$Name,
              " VERBOSITY=WARNING",
              " > addreplace.out")
  ezSystem(cmd)
  ezSystem(paste("samtools", "index", "tmp.bam"))
  cmd = paste0("java -Xmx8g -Djava.io.tmpdir=. -jar ", Sys.getenv("Picard_jar"), 
               " ReorderSam ",
              " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", "tmp.bam",
              " O=", "ordered.bam",
              " REFERENCE=" , ref,
              " VERBOSITY=WARNING",
              " > reorder.out")
  ezSystem(cmd)
  cmd = paste0("java -Xmx8g -Djava.io.tmpdir=. -jar ", Sys.getenv("Picard_jar"), 
               " MarkDuplicates",
              " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", "ordered.bam",
              " O=", localBam,
              " REMOVE_DUPLICATES=true",
              " ASSUME_SORTED=true",
              " METRICS_FILE=" ,"dupmetrics.txt",
              " VERBOSITY=WARNING",
              " >markdup.out")
  ezSystem(cmd)
  ezSystem(paste("samtools", "index", localBam))
  
  
  ########FINDING POSSIBLE INDELS ####
  cmd = paste(gatk, "--num_threads", ezThreads(), "-T  RealignerTargetCreator", 
              "-R", ref, "-o", intervals,
              "-I", localBam,
              ">", "realigntarget.out")
  ezSystem(cmd)
  
  ### REALINGING AROUND POSSIBLE INDELS ###
  cmd = paste(gatk, "-T  IndelRealigner", "-R", ref, "-targetIntervals", intervals,
              "-I", localBam, "-o", realigned, ">", "realign.out")
  ezSystem(cmd)
  
  ### DETECTING VARIANTS ###
  cmd = paste(gatk, "--num_threads", ezThreads(), "-T  UnifiedGenotyper", "-R", ref, "-glm BOTH",
              "-I", realigned, "-o", vcf, ">", "genotype.out")
  ezSystem(cmd)
  
  ### PRELIMINARY EVALUATION ###
  cmd = paste(gatk, "--num_threads", ezThreads(), "-T  VariantEval", "-R", ref,
              "-eval", vcf, "-o", vcfReport, ">", "vcfreport.out")
  ezSystem(cmd)
  return("Success")
}

## calling a haplotype

.gatkHaplotypeCaller = function(bamFile, genomeSeq, knownVariants, vcfOutput){
  gatk = paste("java -Xmx8g -Djava.io.tmpdir=. ", "-jar", Sys.getenv("GATK_jar"))
  cmd = paste(gatk, "-T HaplotypeCaller", "-R", genomeSeq,
              "-I", bamFile,
              "--dbsnp", knownVariants,
              "--recoverDanglingHeads --dontUseSoftClippedBases",
              "_-stand_call_conf 20 --stand_emit_conf 20",
              "-o", vcfOutput) 
}

