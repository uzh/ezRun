ezMethodTophat = function(input=NA, output=NA, param=NA){
  
  ref = getBowtie2Reference(param)
  gtf = param$ezRef["refFeatureFile"]
  bamFile = output$getColumn("BAM")
  trimmedInput = ezMethodFastpTrim(input = input, param = param)
  refBase = sub(".gtf$", "_BOWTIE2Index/transcripts", gtf)
  lockFile = file.path(dirname(refBase), "lock")
  i = 0
  while(file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT){
    ### somebody else builds and we wait at most 180min
    Sys.sleep( 60)
    i = i + 1
  }
  if (file.exists(lockFile)){
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min:", lockFile))
  }
  ## there is no lock file
  refFiles = list.files(dirname(refBase), basename(refBase))
  if(length(refFiles) == 0){
    ## tophat does not hav a command for bulding the index
    ## it can only save the index during a mapping
    ## so we build the index during a mapping of the first 10090 reads
    dir.create(path=dirname(lockFile), recursive=TRUE)
    ezWrite(Sys.info(), con=lockFile)
    gtfOpt = paste("--GTF", gtf)
    head1 = "head1_tmp.fastq"
    ezSystem(paste("head -n 1000", trimmedInput$getColumn("Read1"), ">", head1))
    ## use default mapping with no further options
    cmd = paste("tophat", "-o delme", "--num-threads", param$cores, 
                gtfOpt, "--transcriptome-index", refBase, ref, head1, "2> tophat.log")
    ezSystem(cmd)
    file.remove(lockFile)
    file.remove(head1)
  }
  strandOpt = paste("--library-type", getTuxedoLibraryType(param$strandMode))
  cmd = paste("tophat", "-o .", param$cmdOptions, "-z pigz", "--num-threads", param$cores, strandOpt,
              "--transcriptome-index", refBase, ref, trimmedInput$getColumn("Read1"),
              if(param$paired) trimmedInput$getColumn("Read2"), "2> tophat.log")
  ezSystem(cmd)
  ezSortIndexBam("accepted_hits.bam", basename(bamFile), ram=param$ram,
                 removeBam=TRUE, cores=param$cores)
  
  ## check the strandedness
  bedFile = getReferenceFeaturesBed(param)
  ezSystem(paste("infer_experiment.py", "-r", bedFile, "-i", basename(bamFile), "-s 1000000"), stopOnFailure=FALSE)
  
  ## write an igv link
  if (param$writeIgvLink){ 
    if ("IGV" %in% output$colNames){
      writeIgvHtml(param, output)
    }
  }
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodTophat(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
##' @seealso \code{\link{getBowtie2Reference}}
##' @seealso \code{\link{ezMethodFastpTrim}}
EzAppTophat <-
  setRefClass("EzAppTophat",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodTophat
                  name <<- "EzAppTophat"
                  appDefaults <<- rbind(writeIgvSessionLink=ezFrame(Type="logical", DefaultValue="TRUE", Description="should an IGV link be generated"))
                }
              )
  )
