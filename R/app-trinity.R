###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodTrinity = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){
  
  if (ezIsSpecified(param$samples)){
    input = input$subset(param$samples)
  }
  
  trimmedInput = ezMethodTrim(input = input, param = param)
  param$dataRoot = ""
        
  if (param$paired){
    read1 = paste(trimmedInput$getColumn("Read1"), collapse=",")
    read2 = paste(trimmedInput$getColumn("Read2"), collapse=",")
    readOpt = paste("--left", read1, "--right", read2)
    libOpt = switch(param$strandMode, sense="--SS_lib_type FR", antisense="--SS_lib_type RF", both="")
  } else {
    read1 = paste(trimmedInput$getColumn("Read1"), collapse=",")
    readOpt = paste("--single", read1)
    libOpt = switch(param$strandMode, sense="--SS_lib_type F", antisense="--SS_lib_type R", both="")
  }
  
  cmd = paste("Trinity", "--seqType fq", readOpt,
              "--max_memory", paste0(param$ram, "G"), "--bflyCalculateCPU", ## "--bflyHeapSpaceMax", paste0(round(as.numeric(param$ram)/4), "G"),
              "--CPU", ezThreads(),
              libOpt,
              param$trinityOpt,
              "--output", "trinity", ">", "trinity.stdout")
  ezSystem(cmd)
  pathTranscripts = "trinity/Trinity.fasta"
  
  # Stats and QC
  cmd = paste("TrinityStats.pl", pathTranscripts, ">", basename(output$getColumn("Stats")))
  ezSystem(cmd)
  
  dir.create("abundance")
  # Note that if you change the method, you would also need to change the names of the abundance
  # files, e.g. salmon produces quant.sf while kallisto creates abundance.tsv
  abundanceMethod = "salmon"
  abundancePrefix = "transcript"
  for (nm in trimmedInput$getNames()) {
    sampleDs = trimmedInput$subset(nm)
    if (param$paired) {
      readOpt = paste("--left", sampleDs$getColumn("Read1"), "--right", sampleDs$getColumn("Read2"))
    } else {
      readOpt = paste("--single", sampleDs$getColumn("Read1"))
    }
    cmd = paste("align_and_estimate_abundance.pl", 
                "--transcripts", pathTranscripts,
                "--seqType", "fq",
                "--est_method", abundanceMethod,
                "--prep_reference",
                "--trinity_mode",
                readOpt,
                "--output_dir", file.path("abundance", nm))
    ezSystem(cmd)
  }
  abundanceFiles <- file.path("abundance", trimmedInput$getNames(), "quant.sf")
  cmd = paste("abundance_estimates_to_matrix.pl",
              "--est_method", abundanceMethod,
              "--out_prefix", abundancePrefix,
              "--name_sample_by_basedir",
              paste(abundanceFiles, collapse = " "))
  ezSystem(cmd)
  cmd = paste("$Trinity/util/misc/contig_ExN50_statistic.pl", 
              paste0(abundancePrefix, ".TMM.EXPR.matrix"),
              pathTranscripts,
              ">", basename(output$getColumn("ExN50")))
  ezSystem(cmd)
  
  # Rename output files
  ezSystem(paste("mv", pathTranscripts, basename(output$getColumn("Fasta")) ))
  ezSystem(paste("mv", paste0(abundancePrefix, ".counts.matrix"), basename(output$getColumn("Abundance Counts"))))
  ezSystem(paste("mv", paste0(abundancePrefix, ".TPM.not_cross_norm"), basename(output$getColumn("Abundance TPM"))))
  ezSystem(paste("mv", paste0(abundancePrefix, ".TMM.EXPR.matrix"), basename(output$getColumn("Abundance TMM"))))
  
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodTrinity(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run 
##' @seealso \code{\link{ezMethodTrim}}
EzAppTrinity <-
  setRefClass("EzAppTrinity",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodTrinity
                  name <<- "EzAppTrinity"
                  appDefaults <<- rbind(trinityOpt = ezFrame(Type="character",  DefaultValue="--min_kmer_cov 2",  Description="trinity commandline options"))
                }
              )
  )
