###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


### mapping to transcriptome coordinates


ezMethodBowtie2Transcriptome = function(input=NA, output=NA, param=NA){
  
  ref = getBowtie2TranscriptomeReference(param)
  bamFile = output$getColumn("trBAM")
  trimmedInput = ezMethodTrim(input = input, param = param)
  defOpt = paste("-p", ezThreads())
  strandOpt = switch(param$strandMode,
                     sense="--norc",
                     antisense="--nofw",
                     both="")
  
  cmd = paste("bowtie2", param$cmdOptions, defOpt, strandOpt,
              "-x", ref, if(param$paired) "-1", trimmedInput$getColumn("Read1"), 
              if(param$paired) paste("-2", trimmedInput$getColumn("Read2")),
              "2> bowtie.log", "|", "samtools", "view -S -b -", " > bowtie.bam")
  ezSystem(cmd)
  ezSortIndexBam("bowtie.bam", basename(bamFile), ram=param$ram, removeBam=TRUE, cores=ezThreads())
  return("Success")
}

##' @template getref-template
##' @templateVar methodName Bowtie2
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refBuildDir}{ a character specifying the directory of the reference build.}
##'   \item{ezRef@@refFeatureFile}{ the gtf file with the transcript coordinates}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getBowtie2TranscriptomeReference = function(param){
  
  refBase = ifelse(param$ezRef["refIndex"] == "", 
                   sub(".gtf$", "_BOWTIE2Index/transcripts", param$ezRef["refFeatureFile"]),
                   param$ezRef["refIndex"])
  ## check the ref
  lockFile = file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))){
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con=lockFile)
    wd = getwd()
    setwd(dirname(refBase))

    require(GenomicFeatures)

    txdb = makeTxDbFromGFF(param$ezRef@refFeatureFile,
                           dataSource="FGCZ", taxonomyId = 2759)# organism=organism, chrominfo=NULL) ## taxonomy id is for eukariots
    genomeFa = Rsamtools::FaFile(param$ezRef@refFastaFile)
    
    trSeqFile = file.path(getwd(), "transcriptSeq.fa")
    exonRgList = exonsBy(txdb, by="tx", use.names=TRUE)
    trSeqs = extractTranscriptSeqs(genomeFa, exonRgList)
    Biostrings::writeXStringSet(trSeqs, trSeqFile)
    cmd = paste("bowtie2-build", "--threads", as.numeric(param$cores), "-f", basename(trSeqFile), basename(refBase))
    ezSystem(cmd)
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refBase)))
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
  if (length(refFiles) < 3 ){
    ## we assume the index is built and complete
    stop(paste("index not available: ", refBase))
  }
  return(refBase)
}

##' @template app-template
##' @templateVar method ezMethodBowtie2(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
##' @seealso \code{\link{getBowtie2Reference}}
##' @seealso \code{\link{ezMethodTrim}}
EzAppBowtie2Transcriptome <-
  setRefClass("EzAppBowtie2Transcriptome",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBowtie2Transcriptome
                  name <<- "EzAppBowtie2Transcriptome"
                  ##appDefaults <<- rbind(writeIgvSessionLink=ezFrame(Type="logical", DefaultValue="FALSE", Description="should an IGV link be generated"))
                }
              )
  )
