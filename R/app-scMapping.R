###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

### EzAppSingleCellSTAR
EzAppSingleCellSTAR <- 
  setRefClass("EzAppSingleCellSTAR",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSingleCellSTAR
                  name <<- "EzAppSingleCellSTAR"
                  appDefaults <<- rbind(getJunctions=ezFrame(Type="logical",
                                                             DefaultValue="FALSE",
                                                             Description="should junctions be returned"),
                                        writeIgvSessionLink=ezFrame(Type="logical",
                                                                    DefaultValue="TRUE",
                                                                    Description="should an IGV link be generated"),
                                        markDuplicates=ezFrame(Type="logical", 
                                                               DefaultValue="FALSE", 
                                                               Description="should duplicates be marked with picard"),
                                        checkStrandness=ezFrame(Type="logical", 
                                                                DefaultValue="TRUE", 
                                                                Description="should strandness be checked"),
                                        randomSleep=ezFrame(Type="logical",  
                                                            DefaultValue="FALSE",
                                                            Description="should there be a random sleep to avoid to much network traffic when loading the STAR index"),
                                        twopassMode=ezFrame(Type="logical",
                                                            DefaultValue="FALSE",
                                                            Description="1-pass mapping or basic 2-pass mapping"),
                                        controlSeqs=ezFrame(Type="charVector",
                                                            DefaultValue="",
                                                            Description="control sequences to add")
                  )
                }
              )
  )
### STAR for single cell data: reads in a unmapped bam
###
ezMethodSingleCellSTAR = function(input=NA, output=NA, param=NA){
  
  refDir = getSTARReference(param)
  bamFile = output$getColumn("BAM")
  
  isSingleBam <- !is.na(input$readType()) && input$readType() == "bam"
  if(isSingleBam){
    ## Read 1 is uBam
    fastqInput <- ezMethodBam2Fastq(input = input, param = param,
                                    OUTPUT_PER_RG=TRUE)
    
  }else{
    ## The read data is in CellDataset and input is fastq files
    fastqInput <- EzDataset(file=input$getFullPaths("CellDataset"),
                            dataRoot=DEFAULT_DATA_ROOT)
  }
  
  trimmedInput <- ezMethodTrim(input = fastqInput, param = param)
  
  ## Merge and clean prepross logs
  preprocessLogFns <- paste0(trimmedInput$getNames(), "_preprocessing.log")
  preprocessLogs <- lapply(preprocessLogFns, readLines)
  file.remove(preprocessLogFns)
  writeLines(unlist(preprocessLogs),
             con=paste0(input$getNames(), "_preprocessing.log"))
  
    
  # Clean converted fastqs
  if(isSingleBam){
    file.remove(fastqInput$getFullPaths("Read1"))
    if (param$paired){
      file.remove(fastqInput$getFullPaths("Read2"))
    }
  }
    
  ## fastq to bam
  trimmedBamFn <- tempfile(pattern = "trimmedBam", tmpdir=getwd(),
                           fileext = ".bam")
  if (param$paired){
    fastqs2bam(fastqFns=trimmedInput$getFullPaths("Read1"),
               fastq2Fns=trimmedInput$getFullPaths("Read2"),
               readGroupNames=trimmedInput$getNames(),
               bamFn=trimmedBamFn, mc.cores=param$cores)
    file.remove(c(trimmedInput$getFullPaths("Read1"),
                  trimmedInput$getFullPaths("Read2"))
                )
  }else{
    fastqs2bam(fastqFns=trimmedInput$getFullPaths("Read1"),
               readGroupNames=trimmedInput$getNames(),
               bamFn=trimmedBamFn, mc.cores=param$cores)
    file.remove(trimmedInput$getFullPaths("Read1"))
  }
  ## We can concatenate the fastqs for the aligner
  ## But it will takes more space. So we delete and convert here.
    
  inputTrimmed <- input$copy()
  inputTrimmed$setColumn("Read1", trimmedBamFn)
  inputTrimmed$dataRoot <- NULL
  alignerInput <- ezMethodBam2Fastq(input=inputTrimmed, param = param,
                                    OUTPUT_PER_RG=FALSE)
    
  if(param$cmdOptions == "")
    param$cmdOptions <- "--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --outSAMattributes All"
  
  if (!grepl("outSAMattributes", param$cmdOptions)){
    param$cmdOptions = paste(param$cmdOptions, "--outSAMattributes All")
  }
  
  genomeFn <- param$ezRef@refFastaFile
  
  if(ezIsSpecified(param$controlSeqs)){
    ## control sequences
    controlSeqsLocalFn <- tempfile(pattern="controlSeqs", tmpdir=getwd(),
                                   fileext = ".fa")
    writeXStringSet(getControlSeqs(param$controlSeqs), filepath=controlSeqsLocalFn)
    on.exit(file.remove(controlSeqsLocalFn), add = TRUE)
    
    genomeLocalFn <- tempfile(pattern="genome", tmpdir=getwd(),
                              fileext = ".fa")
    file.copy(from=genomeFn, to=genomeLocalFn)
    writeXStringSet(getControlSeqs(param$controlSeqs), filepath=genomeLocalFn,
                    append=TRUE)
    dictFile = sub(".fa$", ".dict", genomeLocalFn)
    cmd = paste("java -Djava.io.tmpdir=. ", " -jar", 
                Sys.getenv("Picard_jar"), "CreateSequenceDictionary",
                paste0("R=", genomeLocalFn), paste0("O=", dictFile))
    ezSystem(cmd)
    
    genomeFn <- genomeLocalFn
    on.exit(file.remove(c(genomeLocalFn, dictFile)), add=TRUE)
  }
  
  cmd = paste("STAR", " --genomeDir", refDir,  "--sjdbOverhang 150", 
              "--readFilesIn", alignerInput$getColumn("Read1"), 
              if(param$paired) alignerInput$getColumn("Read2"),
              "--twopassMode", ifelse(param$twopassMode, "Basic", "None"),
              if(ezIsSpecified(param$controlSeqs))
                paste("--genomeFastaFiles", controlSeqsLocalFn),
              "--runThreadN", param$cores, param$cmdOptions,
              "--outStd BAM_Unsorted --outSAMtype BAM Unsorted",
              ">  Aligned.out.bam")## writes the output file Aligned.out.bam
  ##"|", "samtools", "view -S -b -", " >", "Aligned.out.bam")
  ezSystem(cmd)
  file.remove(alignerInput$getColumn("Read1"))
  if(param$paired)
    file.remove(alignerInput$getColumn("Read2"))
  
  on.exit(file.remove(c("Log.progress.out", "Log.out", 
                        "Log.std.out")), add=TRUE) ## clean star log files
  on.exit(unlink(c("_STARgenome", "_STARpass1"),
                 recursive = TRUE, force = TRUE), add=TRUE)
  
  ## Merge unmapped and mapped bam to recover the tags
  
  mergeBamAlignments(alignedBamFn="Aligned.out.bam",
                     unmappedBamFn=inputTrimmed$getFullPaths("Read1"),
                     outputBamFn="Aligned.out.merged.bam",
                     fastaFn=genomeFn)
  file.remove("Aligned.out.bam")
  file.rename(from="Aligned.out.merged.bam", to="Aligned.out.bam")
  file.remove(trimmedBamFn)
  
  nSortThreads = min(param$cores, 8)
  ## if the index is loaded in shared memory we have to use only 10% of the scheduled RAM
  if (grepl("--genomeLoad LoadAndKeep", param$cmdOptions)){
    sortRam = param$ram / 10
  } else {
    sortRam = param$ram
  }
  
  file.rename('Log.final.out', to = basename(output$getColumn("STARLog")))
  
  if (!is.null(param$markDuplicates) && param$markDuplicates){
    ezSortIndexBam("Aligned.out.bam", "sorted.bam", ram=sortRam, removeBam=TRUE, 
                   cores=nSortThreads)
    dupBam(inBam="sorted.bam", outBam=basename(bamFile),
           operation="mark", cores=param$cores)
    file.remove("sorted.bam")
  } else {
    ezSortIndexBam("Aligned.out.bam", basename(bamFile), ram=sortRam, 
                   removeBam=TRUE, cores=nSortThreads)
  }
  
  if (param$getJunctions){
    ezSystem(paste("mv SJ.out.tab", basename(output$getColumn("Junctions"))))
    ezSystem(paste("mv Chimeric.out.junction", 
                   basename(output$getColumn("Chimerics"))))
  }else{
    on.exit(file.remove(c("SJ.out.tab", "Chimeric.out.junction",
                          "Chimeric.out.sam")), add=TRUE)
  }
  
  ## check the strandedness
  if (!is.null(param$checkStrandness) && param$checkStrandness){
    cat(Sys.getenv("PATH"), "\n")
    bedFile = getReferenceFeaturesBed(param)
    ezSystem(paste("infer_experiment.py", "-r", bedFile,
                   "-i", basename(bamFile), "-s 1000000"))
  }
  
  ## write an igv link
  if (param$writeIgvSessionLink){ 
    writeIgvSession(genome = getIgvGenome(param), 
                    refBuild=param$ezRef["refBuild"], 
                    file=basename(output$getColumn("IGV Session")),
                    bamUrls = paste(PROJECT_BASE_URL, bamFile, sep="/") )
    writeIgvJnlp(jnlpFile=basename(output$getColumn("IGV Starter")), 
                 projectId = sub("\\/.*", "", bamFile),
                 sessionUrl = paste(PROJECT_BASE_URL, 
                                    output$getColumn("IGV Session"), sep="/"))
  }

  return("Success")
}
