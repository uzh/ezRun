###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @template method-template
##' @templateVar methodName Tophat
##' @seealso \code{\link{EzAppTophat}}
##' @seealso \code{\link{getBowtie2Reference}}
##' @seealso \code{\link{ezMethodTrim}}
ezMethodTophat = function(input=NA, output=NA, param=NA){
  
  Sys.setenv(PATH=paste(BOWTIE2_DIR, BOWTIE_DIR, dirname(SAMTOOLS), Sys.getenv("PATH"), sep=":"))
  ref = getBowtie2Reference(param)
  gtf = param$ezRef["refFeatureFile"]
  bamFile = output$getColumn("BAM")
  trimmedInput = ezMethodTrim(input = input, param = param)
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
    cmd = paste(file.path(TOPHAT_DIR, "tophat"), "-o delme", "--num-threads", ezThreads(), 
                gtfOpt, "--transcriptome-index", refBase, ref, head1, "2> tophat.log")
    ezSystem(cmd)
    file.remove(lockFile)
    file.remove(head1)
  }
  strandOpt = paste("--library-type", getTuxedoLibraryType(param$strandMode))
  cmd = paste(file.path(TOPHAT_DIR, "tophat"), "-o .", param$cmdOptions, "-z pigz", "--num-threads", ezThreads(), strandOpt,
              "--transcriptome-index", refBase, ref, trimmedInput$getColumn("Read1"),
              ifelse(param$paired, trimmedInput$getColumn("Read2"), ""), "2> tophat.log")
  ezSystem(cmd)
  ezSortIndexBam("accepted_hits.bam", basename(bamFile), ram=param$ram, removeBam=TRUE, cores=ezThreads())
  
  ## write an igv link
  if (param$writeIgvSessionLink){
    writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
                    bamUrls = paste(PROJECT_BASE_URL, bamFile, sep="/") )
    writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", bamFile),
                 sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
  }
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodTophat()
##' @seealso \code{\link{ezMethodTophat}}
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

##' @template method-template
##' @templateVar methodName Bowtie2
##' @seealso \code{\link{EzAppBowtie2}}
##' @seealso \code{\link{getBowtie2Reference}}
##' @seealso \code{\link{ezMethodTrim}}
ezMethodBowtie2 = function(input=NA, output=NA, param=NA){
  
  ref = getBowtie2Reference(param)
  #Sys.setenv(BOWTIE2_INDEXES=dirname(ref))
  #message("bowtie Dir:", Sys.getenv("BOWTIE2_INDEXES"))
  bamFile = output$getColumn("BAM")
  trimmedInput = ezMethodTrim(input = input, param = param)
  defOpt = paste("-p", ezThreads())
  
  cmd = paste(file.path(BOWTIE2_DIR, "bowtie2"), param$cmdOptions, defOpt, 
              "-x", ref, trimmedInput$getColumn("Read1"), ifelse(param$paired, trimmedInput$getColumn("Read2"), ""),
              "2> bowtie.log", "|", SAMTOOLS, "view -S -b -", " > bowtie.bam")
  ezSystem(cmd)
  ezSortIndexBam("bowtie.bam", basename(bamFile), ram=param$ram, removeBam=TRUE, cores=ezThreads())
  
  ## write an igv link
  if (param$writeIgvSessionLink){
    writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
                    bamUrls = paste(PROJECT_BASE_URL, bamFile, sep="/") )
    writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", bamFile),
                 sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName Bowtie2
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refBuildDir}{ a character specifying the directory of the reference build.}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getBowtie2Reference = function(param){
  
  refBase = ifelse(param$ezRef["refIndex"] == "", 
                   file.path(param$ezRef["refBuildDir"], "Sequence/BOWTIE2Index/genome"),
                   param$ezRef["refIndex"])
  ## check the ref
  lockFile = file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))){
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con=lockFile)
    wd = getwd()
    setwd(dirname(refBase))
    
    fastaFile = param$ezRef["refFastaFile"]
    ezSystem(paste("ln -s", fastaFile, "."))
    cmd = paste(file.path(BOWTIE2_DIR, "bowtie2-build"), "-f", basename(fastaFile), basename(refBase))
    ezSystem(cmd)
    #ezWriteElapsed(job, "done")
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
##' @templateVar method ezMethodBowtie2()
##' @seealso \code{\link{ezMethodBowtie2}}
EzAppBowtie2 <-
  setRefClass("EzAppBowtie2",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBowtie2
                  name <<- "EzAppBowtie2"
                  appDefaults <<- rbind(writeIgvSessionLink=ezFrame(Type="logical", DefaultValue="TRUE", Description="should an IGV link be generated"))
                }
              )
  )

##' @template method-template
##' @templateVar methodName Bowtie
##' @seealso \code{\link{EzAppBowtie}}
##' @seealso \code{\link{getBowtieReference}}
##' @seealso \code{\link{ezMethodTrim}}
ezMethodBowtie = function(input=NA, output=NA, param=NA){
    
  ref = getBowtieReference(param)
  bamFile = output$getColumn("BAM")  
  trimmedInput = ezMethodTrim(input = input, param = param)
  defOpt = paste("--chunkmbs 256", "--sam", "-p", ezThreads())
    cmd = paste(file.path(BOWTIE_DIR, "bowtie"), param$cmdOptions, defOpt, 
              ref, trimmedInput$getColumn("Read1"), ifelse(param$paired, trimmedInput$getColumn("Read2"), ""),
              "2> bowtie.log", "|", SAMTOOLS, "view -S -b -", " > bowtie.bam")
  ezSystem(cmd)
  ezSortIndexBam("bowtie.bam", basename(bamFile), ram=param$ram, removeBam=TRUE, cores=ezThreads())
  
  ## write an igv link
  if (param$writeIgvSessionLink){ 
    writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
                    bamUrls = paste(PROJECT_BASE_URL, bamFile, sep="/") )
    writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", bamFile),
                 sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName Bowtie
##' @inheritParams getBowtie2Reference
getBowtieReference = function(param){
  
  refBase = ifelse(param$ezRef["refIndex"] == "", 
                   file.path(param$ezRef["refBuildDir"], "Sequence/BOWTIEIndex/genome"),
                   param$ezRef["refIndex"])
  ## check the ref
  lockFile = file.path(dirname(refBase), "lock")
  if (!file.exists(dirname(refBase))){
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refBase))
    ezWrite(Sys.info(), con=lockFile)
    wd = getwd()
    setwd(dirname(refBase))
    
    fastaFile = param$ezRef["refFastaFile"]
    #job = ezJobStart("bowtie index")
    ezSystem(paste("ln -s", fastaFile, "."))
    cmd = paste(file.path(BOWTIE_DIR, "bowtie-build"), "-f", basename(fastaFile), basename(refBase))
    ezSystem(cmd)
    #ezWriteElapsed(job, "done")
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
##' @templateVar method ezMethodBowtie()
##' @seealso \code{\link{ezMethodBowtie}}
EzAppBowtie <-
  setRefClass("EzAppBowtie",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBowtie2
                  name <<- "EzAppBowtie"
                  appDefaults <<- rbind(writeIgvSessionLink=ezFrame(Type="logical", DefaultValue="TRUE", Description="should an IGV link be generated"))
                }
              )
  )

##' @template method-template
##' @templateVar methodName STAR
##' @seealso \code{\link{EzAppSTAR}}
##' @seealso \code{\link{getSTARReference}}
##' @seealso \code{\link{ezMethodTrim}}
ezMethodSTAR = function(input=NA, output=NA, param=NA){

  refDir = getSTARReference(param)
  bamFile = output$getColumn("BAM")
  trimmedInput = ezMethodTrim(input = input, param = param)
  
  cmd = paste(STAR, "--genomeLoad NoSharedMemory --genomeDir", refDir,  "--sjdbOverhang 150", "--readFilesIn",
              trimmedInput$getColumn("Read1"), ifelse(param$paired, trimmedInput$getColumn("Read2"), ""),
              "--runThreadN", ezThreads(), param$cmdOptions, "--outStd BAM_Unsorted --outSAMtype BAM Unsorted",
              ">  Aligned.out.bam")## writes the output file Aligned.out.bam
  ##"|", SAMTOOLS, "view -S -b -", " >", "Aligned.out.bam")
  ezSystem(cmd)
  nSortThreads = min(ezThreads(), 8)
  ezSortIndexBam("Aligned.out.bam", basename(bamFile), ram=param$ram, removeBam=TRUE, cores=nSortThreads)
  if (param$getChimericJunctions){
    ezSystem(paste("mv Chimeric.out.junction", basename(output$getColumn("Chimerics"))))
  }  
  ## write an igv link
  if (param$writeIgvSessionLink){ 
    writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
                    bamUrls = paste(PROJECT_BASE_URL, bamFile, sep="/") )
    writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", bamFile),
                 sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName STAR
##' @param param a list of parameters:
##' \itemize{
##'   \item{ezRef@@refIndex}{ a character specifying the location of the index that is used in the alignment.}
##'   \item{ezRef@@refFeatureFile}{ a character specifying the file path to the annotation feature file (.gtf).}
##'   \item{ezRef@@refChromDir}{ a character specifying the file path to the directory of the chromosome information.}
##'   \item{ram}{ an integer specifying how many gigabytes of RAM to use.}
##'   \item{ezRef@@refFastaFile}{ a character specifying the file path to the fasta file.}
##' }
getSTARReference = function(param){
  
  if (ezIsSpecified(param$ezRef["refIndex"])){
    refDir = param$ezRef["refIndex"]
  } else {
    if (!ezIsSpecified(param$ezRef["refFeatureFile"])){
      stop("not refFeatureFile defined")
    }
    refDir = sub(".gtf$", "_STARIndex", param$ezRef["refFeatureFile"])
  }
  
  lockFile = file.path(refDir, "lock")
  i = 0
  while(file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT){
    ### somebody else builds and we wait
    Sys.sleep(90)
    i = i + 1
  }
  if(file.exists(lockFile)){
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  ## there is no lock file
  refFiles = list.files(refDir)
  if(length(refFiles) > 0){
    ## we assume the index is built and complete
    return(refDir)
  }
  
  ## no lock file and no refFiles, so we build the reference
  dir.create(refDir)
  ezWrite(Sys.info(), con=lockFile)
  
  if (length(list.files(param$ezRef["refChromDir"], ".fa$")) > 50){
    binOpt = "--genomeChrBinNbits 16"
  } else {
    binOpt = ""
  }
  
  job = ezJobStart("STAR genome build")
  cmd = paste(STAR, "--runMode genomeGenerate --genomeDir", refDir, binOpt,
              "--limitGenomeGenerateRAM", format(param$ram * 1e9, scientific=FALSE),
              "--genomeFastaFiles", param$ezRef["refFastaFile"], 
              "--sjdbGTFfile", param$ezRef["refFeatureFile"], "--sjdbOverhang 150", "--runThreadN", ezThreads())
  ezSystem(cmd)
  file.remove(lockFile)
  ezWriteElapsed(job, "done")
  file.remove("Log.out")
  return(refDir)
}

##' @template app-template
##' @templateVar method ezMethodSTAR()
##' @seealso \code{\link{ezMethodSTAR}}
EzAppSTAR <- 
  setRefClass("EzAppSTAR",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSTAR
                  name <<- "EzAppSTAR"
                  appDefaults <<- rbind(getChimericJunctions=ezFrame(Type="logical",  DefaultValue="FALSE",	Description="should chimeric reads be returned"),
                                        writeIgvSessionLink=ezFrame(Type="logical", DefaultValue="TRUE", Description="should an IGV link be generated"))
                }
              )
  )

##' @template method-template
##' @templateVar methodName BWA
##' @seealso \code{\link{EzAppBWA}}
##' @seealso \code{\link{getBWAReference}}
##' @seealso \code{\link{ezMethodTrim}}
ezMethodBWA = function(input=NA, output=NA, param=NA){
  
  refIdx = getBWAReference(param)
  bamFile = output$getColumn("BAM")
  trimmedInput = ezMethodTrim(input = input, param = param)
  cmd = paste(BWA, param$algorithm, param$cmdOptions, "-t", ezThreads(),
              refIdx, trimmedInput$getColumn("Read1"), ">", "read1.sai", "2> bwa.log")
  ezSystem(cmd)
  if (param$algorithm == "aln"){
    if (param$paired){
      cmd = paste(BWA, param$algorithm, param$cmdOptions, "-t", ezThreads(),
                  refIdx, trimmedInput$getColumn("Read2"), ">", "read2.sai", "2> bwa.log")
      ezSystem(cmd)
      cmd = paste(BWA, "sampe", refIdx, "read1.sai", "read2.sai", trimmedInput$getColumn("Read1"), trimmedInput$getColumn("Read2"), "2> bwa.log", "|",
                  SAMTOOLS, "view -S -b -", " > aligned.bam", "2> bwa.log")
      ezSystem(cmd)
    } else {
      cmd = paste(BWA, "samse", refIdx, "read1.sai", trimmedInput$getColumn("Read1"), "|",
                  SAMTOOLS, "view -S -b -", " > aligned.bam", "2> bwa.log")
      ezSystem(cmd)
    }
  } else {
    if(param$algorithm == "bwasw" && param$paired){
      stop("paired is not supported for algorithm bwasw")
    }
    cmd = paste(BWA, param$algorithm, param$cmdOptions, "-t", ezThreads(),
                refIdx, trimmedInput$getColumn("Read1"), ifelse(param$paired, trimmedInput$getColumn("Read2"), ""),
                "|", SAMTOOLS, "view -S -b -", " > aligned.bam", "2> bwa.log")
    ezSystem(cmd)
  }
  ezSortIndexBam("aligned.bam", basename(bamFile), ram=param$ram, removeBam=TRUE, cores=ezThreads())
  
  ## write an igv link
  if (param$writeIgvSessionLink){ 
    writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
                    bamUrls = paste(PROJECT_BASE_URL, bamFile, sep="/") )
    writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", bamFile),
                 sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
  }
  return("Success")
}

##' @template getref-template
##' @templateVar methodName BWA
##' @inheritParams getBowtie2Reference
getBWAReference = function(param){
  
  refPath = ifelse(param$ezRef["refIndex"] == "", 
                   file.path(param$ezRef["refBuildDir"], "Sequence/BWAIndex/genome.fa"),
                   param$ezRef["refIndex"])
  ## check the ref
  lockFile = file.path(dirname(refPath), "lock")
  if (!file.exists(dirname(refPath))){
    ## no lock file and no refFiles, so we build the reference
    dir.create(dirname(refPath))
    ezWrite(Sys.info(), con=lockFile)
    wd = getwd()
    setwd(dirname(refPath))
    
    fastaFile = param$ezRef["refFastaFile"]
    #job = ezJobStart("bwa index")
    ezSystem(paste("ln -s", fastaFile, "."))
    cmd = paste(BWA, "index", "-a", "bwtsw", basename(fastaFile))
    ezSystem(cmd)
    #ezWriteElapsed(job, "done")
    setwd(wd)
    file.remove(lockFile)
  }
  stopifnot(file.exists(dirname(refPath)))
  i = 0
  while(file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT){
    ### somebody else builds and we wait
    Sys.sleep( 60)
    i = i + 1
  }
  if (file.exists(lockFile)){
    stop(paste("reference building still in progress after", INDEX_BUILD_TIMEOUT, "min"))
  }
  if (!file.exists(paste0(refPath, ".sa"))){
    stop(paste("sa index not found for:", refPath))
  }  
  return(refPath)
}

##' @template app-template
##' @templateVar method ezMethodBWA()
##' @seealso \code{\link{ezMethodBWA}}
EzAppBWA <- 
  setRefClass("EzAppBWA",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBWA
                  name <<- "EzAppBWA"
                  appDefaults <<- rbind(algorithm=ezFrame(Type="character",  DefaultValue="mem",  Description="bwa's alignment algorithm. One of aln, bwasw, mem."),
                                        writeIgvSessionLink=ezFrame(Type="logical", DefaultValue="TRUE", Description="should an IGV link be generated"))
                }
              )
  )

##' @template method-template
##' @templateVar methodName Bismark
##' @seealso \code{\link{EzAppBismark}}

ezMethodBismark = function(input=NA, output=NA, param=NA){
  ##TODO: create reference if not existing
  ref = file.path('/srv/GT/reference',dirname(dirname(param[['refBuild']])),'Sequence/WholeGenomeFasta')
  bamFile = output$getColumn("BAM")
  trimmedInput = ezMethodTrim(input = input, param = param)
  defOpt = paste("-p", max(2,ezThreads()/2))
  cmd = paste(file.path(BISMARK_DIR, "bismark"), param$cmdOptions ,defOpt, ref, '-1',
              trimmedInput$getColumn("Read1"), ifelse(param$paired, paste('-2',trimmedInput$getColumn("Read2")), ""),  
              "2> bismark.log")
  
  ezSystem(cmd)
  bamFileNameBismark = list.files('.',pattern='bam$')
  reportFileNameBismark = list.files('.',pattern='report.txt$')
  ezSystem(paste('mv ', reportFileNameBismark, paste0(names(bamFile),'.report.txt')))
  cmd = paste(SAMTOOLS, "view -S -b ",bamFileNameBismark, " > bismark.bam")
  ezSystem(cmd)
  ezSortIndexBam("bismark.bam", basename(bamFile), ram=param$ram, removeBam=TRUE, cores=ezThreads())
              
              ## write an igv link
              #if (param$writeIgvSessionLink){
              #  writeIgvSession(genome = getIgvGenome(param), refBuild=param$ezRef["refBuild"], file=basename(output$getColumn("IGV Session")),
              #                  bamUrls = paste(PROJECT_BASE_URL, bamFile, sep="/") )
              #  writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), projectId = sub("\\/.*", "", bamFile),
              #               sessionUrl = paste(PROJECT_BASE_URL, output$getColumn("IGV Session"), sep="/"))
              #}
              return("Success")
}

##' @template app-template
##' @templateVar method ezMethodBismark()
##' @seealso \code{\link{ezMethodBismark}}
EzAppBismark <-
  setRefClass("EzAppBismark",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodBismark
                  name <<- "EzAppBismark"
                  #appDefaults <<- ''
                }
              )
  )