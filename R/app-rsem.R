###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName RSEM
##' @seealso \code{\link{EzAppRSEM}}
##' @seealso \code{\link{getRSEMReference}}
##' @seealso \code{\link{ezMethodTrim}}
ezMethodRSEM = function(input=NA, output=NA, param=NA){
  
  sampleName = input$getNames()
  Sys.setenv(PATH=paste(BOWTIE_DIR, dirname(SAMTOOLS), Sys.getenv("PATH"), sep=":"))  
  ref = getRSEMReference(param)
  
  opt = param$cmdOptions
  if (ezIsSpecified(param$"bowtie-e")){
    opt = paste(param$cmdOptions, "--bowtie-e", param$"bowtie-e")
  }
  samtoolsSortMem = paste0(round( (as.numeric(param$ram) / ezThreads() /2)*1000), "M")
  ciMemory = round((as.numeric(param$ram) - 4) * 1000)
  opt = paste(opt, "-p", ezThreads(), "--bowtie-path", BOWTIE_DIR, "--samtools-sort-mem", samtoolsSortMem, "--ci-memory", ciMemory)
  
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
    cmd = paste(file.path(RSEM_DIR, "rsem-calculate-expression"), opt, strandOpt,
                "--paired-end", trimmedInput$getColumn("Read2"), trimmedInput$getColumn("Read1"),
                ref, sampleName, "2> rsem.stderr", "> rsem.stdout")
  } else {
    cmd = paste(file.path(RSEM_DIR, "rsem-calculate-expression"), opt, strandOpt,
                trimmedInput$getColumn("Read1"),
                ref, sampleName, "2> rsem.stderr", "> rsem.stdout")
  }
  ezSystem(cmd)
  if (!is.null(param$keepBam) && param$keepBam){
    localBam = paste0(sampleName, ".genome.sorted.bam")
    localBai = paste0(sampleName, ".genome.sorted.bam.bai")
    bamFile = basename(output$getColumn("BAM"))
    baiFile = basename(output$getColumn("BAI"))
    ezSystem(paste("mv", localBam, bamFile))
    ezSystem(paste("mv", localBai, baiFile))
  }
  
  transcriptResult = ezRead.table(paste0(sampleName, ".isoforms.results"), header=TRUE)
  if (ezIsSpecified(param$transcriptFasta)){
    result = transcriptResult[ , c("gene_id", "length")]
    colnames(result) = c("gene_id", "width")
  } else {
    result = ezRead.table(param$ezRef["refAnnotationFile"], colClasses="character")
    result = result[ , intersect(colnames(result), "gene_id", "width", "gc"), drop=FALSE]
  }
  transcriptResult = transcriptResult[rownames(result), ]
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
  ezWrite.table(result, file=basename(output$getColumn("Count")))
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodRSEM()
##' @seealso \code{\link{ezMethodRSEM}}
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

##' @template getref-template
##' @templateVar methodName RSEM
##' @param param a list of parameters:
##' \itemize{
##'   \item{transcriptFasta}{ an optional character specifying the path to a fasta file. If specified, the reference will be prepared using it.}
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refFeatureFile}{ a character specifying the path to the annotation feature file (.gtf).}
##'   \item{ezRef@@refFastaFile}{ a character specifying the path to the fasta file.}
##' }
getRSEMReference = function(param){
  
  if (ezIsSpecified(param$transcriptFasta)){
    refBase = file.path(getwd(), "RSEMIndex/transcripts") #paste0(file_path_sans_ext(param$trinityFasta), "_RSEMIndex/transcripts")
  } else {
    refBase = ifelse(param$ezRef["refIndex"] == "", 
                     sub(".gtf$", "_RSEMIndex/transcripts", param$ezRef["refFeatureFile"]),
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
  
  prepareReference = file.path(RSEM_DIR, "rsem-prepare-reference")
  job = ezJobStart("rsem build")
  if (ezIsSpecified(param$transcriptFasta)){
    ## check if the file comes from trinity
    trxNames = sub(" .*", "", names(readDNAStringSet(param$transcriptFasta, nrec=100)))
    if (all(grepl("^comp.+_c.+_s.+", trxNames))){
      mapFile = paste0(refBase, "_geneMap.txt")
      cmd = paste(file.path(RSEM_DIR, "extract-transcript-to-gene-map-from-trinity"), param$transcriptFasta, mapFile)
      ezSystem(cmd)
      cmd = paste(prepareReference, "--bowtie", "-q", "--bowtie-path", BOWTIE_DIR, "--transcript-to-gene-map", mapFile, param$transcriptFasta, "transcripts")
      ezSystem(cmd)    
    } else {
      cmd = paste(prepareReference, "--bowtie", "-q", "--bowtie-path", BOWTIE_DIR, param$transcriptFasta, "transcripts")
      ezSystem(cmd)    
    }
  } else{
    cmd = paste(prepareReference, "--bowtie", "--gtf", param$ezRef["refFeatureFile"], 
                "--bowtie-path", BOWTIE_DIR, param$ezRef["refFastaFile"],
                "transcripts")
    ezSystem(cmd)
    ezSystem(paste("ln -s", param$ezRef["refFeatureFile"], "."))
  }
  ezWriteElapsed(job, "done")
  return(refBase)
}
