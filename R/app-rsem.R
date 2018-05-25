###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodRSEM = function(input=NA, output=NA, param=NA){
  require(readr)
  require(dplyr)
  
  sampleName = input$getNames()
  ref = getRSEMReference(param)
  
  opt = param$cmdOptions
  if (ezIsSpecified(param$"bowtie-e")){
    opt = paste(param$cmdOptions, "--bowtie-e", param$"bowtie-e")
  }
  samtoolsSortMem = paste0(round( (as.numeric(param$ram) / ezThreads() /2)*1000), "M")
  ciMemory = round((as.numeric(param$ram) - 4) * 1000)
  opt = paste(opt, "-p", ezThreads(), "--sort-bam-memory-per-thread", samtoolsSortMem, 
              "--ci-memory", ciMemory)
  
  trimmedInput = ezMethodTrim(input = input, param = param)
  
  strandOpt = switch(param$strandMode,
                     "sense"="--forward-prob 1.0",
                     "antisense"="--forward-prob 0.0",
                     "both"="")
  if (!is.null(param$keepBam) && param$keepBam){
    opt = sub("--no-bam-output", "", opt)
    opt = paste("--output-genome-bam", opt)
  }
  if (param$paired){
    opt = sub('--sort-bam-by-read-name', '', opt)
    cmd = paste("rsem-calculate-expression", opt, strandOpt,
                "--paired-end", trimmedInput$getColumn("Read1"), trimmedInput$getColumn("Read2"),
                ref, sampleName, "2> rsem.stderr", "> rsem.stdout")
  } else {
    cmd = paste("rsem-calculate-expression", opt, strandOpt,
                trimmedInput$getColumn("Read1"),
                ref, sampleName, "2> rsem.stderr", "> rsem.stdout")
  }
  ezSystem(cmd)
  
  file.remove(trimmedInput$getColumn("Read1"))
  if(param$paired)
    file.remove(trimmedInput$getColumn("Read2"))
  
  file.remove(c("rsem.stderr", "rsem.stdout"))
  if (!is.null(param$keepBam) && param$keepBam){
    localBam = paste0(sampleName, ".genome.bam")
    bamFile = basename(output$getColumn("BAM"))
    ezSortIndexBam(localBam, bamFile, ram = param$ram, min(param$cores, 8))
  }
  
  #transcriptResult = ezRead.table(paste0(sampleName, ".isoforms.results"), header=TRUE)
  transcriptResult <- read_tsv(paste0(sampleName, ".isoforms.results"))
  if (ezIsSpecified(param$transcriptFasta)){
    #result = transcriptResult[ , c("transcript_id", "gene_id", "length")]
    result <- select(transcriptResult, "transcript_id", "gene_id", "length")
    #colnames(result) = c("gene_id", "width")
  } else {
    #result = ezRead.table(param$ezRef["refAnnotationFile"], colClasses="character")
    result <- read_tsv(param$ezRef["refAnnotationFile"])
    stopifnot(all(transcriptResult$transcript_id %in% result$transcript_id))
    #stopifnot(rownames(transcriptResult) %in% rownames(result))
    result <- result[match(transcriptResult$transcript_id, result$transcript_id),
                     intersect(colnames(result),
                               c("transcript_id", "gene_id", "width", "gc"))]
  }
  result$transcriptCount = transcriptResult$expected_count
  #   result$transcriptExpression = transcriptResult$TPM /1e6 * sum(transcriptResult$expected_count)
  result$TPM = transcriptResult$TPM
  result$FPKM = transcriptResult$FPKM
  if (!is.null(transcriptResult$posterior_mean_count)){
    result$transcriptCountPosteriorEstimate = transcriptResult$posterior_mean_count
    #     result$transcriptExpressionLowerBound = transcriptResult$TPM_ci_lower_bound /1e6 * sum(transcriptResult$pme_expected_count)
    #     result$transcriptExpressionUpperBound = transcriptResult$TPM_ci_upper_bound /1e6 * sum(transcriptResult$pme_expected_count)
    result$TPM_PosteriorEstimate = transcriptResult$pme_TPM
    result$FPKM_PosteriorEstimate = transcriptResult$pme_FPKM
  }
  write_tsv(result, path=basename(output$getColumn("Count")))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodRSEM(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
##' @seealso \code{\link{getRSEMReference}}
##' @seealso \code{\link{ezMethodTrim}}
EzAppRSEM <-
  setRefClass("EzAppRSEM",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodRSEM
                  name <<- "EzAppRSEM"
                  appDefaults <<- rbind("bowtie-e"=ezFrame(Type="integer",  DefaultValue="200",  Description="maximum sum of mismatch base qualities"),
                                    keepBam=ezFrame(Type="logical",  DefaultValue=FALSE,  Description="if the alignment files should be kept"),
                                    transcriptFasta=ezFrame(Type="character", DefaultValue='', Description="full path to a transcript fasta file (e.g. trinity output); must be in a writeable directory;"))
                  
                }
              )
  )

getRSEMReference = function(param){
  if (ezIsSpecified(param$transcriptFasta)){
    refBase = file.path(getwd(), "RSEMIndex/transcripts") #paste0(file_path_sans_ext(param$trinityFasta), "_RSEMIndex/transcripts")
  } else {
    if(ezIsSpecified(param$transcriptTypes)){
      rsemBase <- paste(sort(param$transcriptTypes), collapse="-")
      ## This is a combination of transcript types to use.
    }else{
      rsemBase <- ""
    }
    refBase = ifelse(param$ezRef["refIndex"] == "", 
                     sub(".gtf$", paste0("_", rsemBase, "_RSEMIndex/transcripts"),
                         param$ezRef["refFeatureFile"]),
                     param$ezRef["refIndex"])
  }
  lockFile = file.path(dirname(refBase), "lock")
  i = 0
  while(file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT){
    ### somebody else builds and we wait
    Sys.sleep( 60)
    i = i + 1
  }
  if (file.exists(lockFile)){
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles = list.files(dirname(refBase), basename(refBase))
  if (length(refFiles) > 0 ){
    ## we assume the index is built and complete
    return(refBase)
  }
  ## we have to build the reference
  wd = getwd()
  dir.create(dirname(refBase))
  setwd(dirname(refBase))
  ezWrite(Sys.info(), con=lockFile)
  on.exit(file.remove(lockFile))
  on.exit(setwd(wd), add=TRUE)
  
  prepareReference = "rsem-prepare-reference"
  job = ezJobStart("rsem build")
  if (ezIsSpecified(param$transcriptFasta)){
    ## check if the file comes from trinity
    trxNames = sub(" .*", "", names(readDNAStringSet(param$transcriptFasta, nrec=100)))
    if (all(grepl("^comp.+_c.+_s.+", trxNames))){
      mapFile = paste0(refBase, "_geneMap.txt")
      cmd = paste("extract-transcript-to-gene-map-from-trinity", param$transcriptFasta, mapFile)
      ezSystem(cmd)
      cmd = paste(prepareReference, "--bowtie", "-q", "--transcript-to-gene-map", mapFile, param$transcriptFasta, "transcripts")
      ezSystem(cmd)    
    } else {
      cmd = paste(prepareReference, "--bowtie", "-q", param$transcriptFasta, "transcripts")
      ezSystem(cmd)    
    }
  } else{
    if(ezIsSpecified(param$transcriptTypes)){
      seqAnno = ezFeatureAnnotation(param$ezRef@refAnnotationFile,
                                    dataFeatureType="transcript")
      transcriptsUse = rownames(seqAnno)[seqAnno$type %in% param$transcriptTypes]
      
      gtf <- ezReadGff(param$ezRef@refFeatureFile)
      transcripts <- ezGffAttributeField(gtf$attributes,
                                         field="transcript_id", 
                                         attrsep="; *", valuesep=" ")
      gtf = gtf[transcripts %in% transcriptsUse, ]
      gtfFile <- sub(".gtf$", paste0("_", rsemBase, ".gtf"),
                     param$ezRef["refFeatureFile"])
      write.table(gtf, gtfFile, quote=FALSE, sep="\t", 
                  row.names=FALSE, col.names=FALSE)
    }else{
      gtfFile <- param$ezRef["refFeatureFile"]
    }
    cmd = paste(prepareReference, "--bowtie", "--gtf", gtfFile, #param$ezRef["refFeatureFile"],
                param$ezRef["refFastaFile"],
                "transcripts")
    ezSystem(cmd)
    ezSystem(paste("ln -s", gtfFile, "."))
  }
  ezWriteElapsed(job, "done")
  return(refBase)
}
