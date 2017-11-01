###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodSingleCellCounts = function(input=NA, output=NA, param=NA, htmlFile="00index.html"){

  input = EzDataset(file=input$getFullPaths("ReadDataset"), dataRoot=param$dataRoot)
  if (ezIsSpecified(param$samples)){
    input = input$subset(param$samples)
  }
  countDir = basename(output$getColumn('CountFolder'))
  dir.create(countDir)
  bamMeta = input$meta[ , !input$columnHasTag("File")]
  bamMeta[["BAM [File]"]] = paste0(getwd(), "/", input$getNames(), "/", input$getNames(), ".bam")
  bamMeta[["BAI [File]"]] = paste0(getwd(), "/", input$getNames(), "/", input$getNames(), ".bam.bai")
  bamMeta[["Read Count"]] = ceiling(bamMeta[["Read Count"]] / param$subsampleReads)
  bamMeta[["STARLog [File]"]] = paste0(getwd(), "/", input$getNames(), "/", input$getNames(), "_STAR.log")
  bamOutput = EzDataset(meta=bamMeta, param$dataRoot)
  bamParam = param
  bamParam$mail = ""
  countMeta = input$meta[ , !input$columnHasTag("File")]
  countMeta[['Count [File]']] = paste0(rownames(countMeta), '.txt')
  countMeta[['Stats [File]']] = paste0(rownames(countMeta), '-stats.txt')
  ezWrite.table(countMeta, file=basename(output$getColumn('CountDataset')), head='Name')
  countDs = EzDataset(meta=countMeta)
  
  switch(param$mapMethod,
         STAR={
           refDir = getSTARReference(param)
           mappingApp = EzAppSTAR$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --outSAMattributes All")
           if (!grepl("genomeLoad", bamParam$cmdOptions)){
             bamParam$cmdOptions = paste("--genomeLoad LoadAndKeep" , bamParam$cmdOptions)
           }
           # LoadAndKeep is not compatible with twopassMode
           bamParam$twopassMode = F
           if (!grepl("outSAMattributes", bamParam$cmdOptions)){
             bamParam$cmdOptions = paste(bamParam$cmdOptions, "--outSAMattributes All")
           }

           on.exit({
             if (ezIsSpecified(param$ezRef["refIndex"])){
               refDir = param$ezRef["refIndex"]
             } else {
               if (!ezIsSpecified(param$ezRef["refFeatureFile"])){
                 stop("refFeatureFile not defined")
               }
               refDir = sub(".gtf$", "_STARIndex", param$ezRef["refFeatureFile"])
             }
             ezSystem(paste("STAR", '--genomeDir', refDir, '--genomeLoad Remove'))
           })
         },
         bowtie={
           mappingApp = EzAppBowtie$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "")
         },
         bowtie2={
           mappingApp = EzAppBowtie2$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "--no-unal")
         },
         tophat={
           mappingApp = EzAppTophat$new()
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "--mate-inner-dist 100 --mate-std-dev 150")
         },
         "bwa-mem"={
           mappingApp = EzAppBWA$new()
           bamParam$algorithm = "mem"
           bamParam$cmdOptions = ifelse(bamParam$mapOptions != "", bamParam$mapOptions,
                                        "")
         },
         stop("unsupported mapMethod: ", param$mapMethod)
  )
  
  
  for (nm in input$getNames()){
    setwdNew(nm)
    readI = input$subset(nm)
    bamI = bamOutput$subset(nm)
    countI = countDs$subset(nm)
    mappingApp$run(input=readI, output=bamI, param=bamParam)
    EzAppFeatureCounts$new()$run(input=bamI, output=countI, param=bamParam)
    file.rename(basename(countI$getColumn("Count")), paste0("../", countDir, "/", basename(countI$getColumn("Count"))))
    file.rename(basename(countI$getColumn("Stats")), paste0("../", countDir, "/", basename(countI$getColumn("Stats"))))
    ezSystem(paste("mv", basename(bamI$getColumn("BAM")), file.path("..", countDir) ))
    ezSystem(paste("mv", basename(bamI$getColumn("BAI")), file.path("..", countDir) ))
    ezSystem(paste("mv", basename(bamI$getColumn("STARLog")), file.path("..", countDir) ))
    setwd('..')
    unlink(nm, recursive=TRUE, force=TRUE)
    
  }
  
  countFiles = paste0(countDir, "/", basename(countDs$getColumn("Count")))
  templateFile = file.path(countDir, countDs$subset( input$getNames()[1] )$getColumn("Count"))
  x = ezRead.table(templateFile)
  counts = ezMatrix(0, rows=rownames(x), cols=countDs$getNames())
  for (nm in input$getNames()) {
    countFile = file.path(countDir, countDs$subset(nm)$getColumn("Count"))
    x = ezRead.table(countFile)
    stopifnot(setequal(rownames(x), rownames(counts)))
    counts[rownames(x), nm] = x[ , "matchCounts"]
  }
  ezWrite.table(counts, head=paste0(param$featureLevel, "_id"), 
                file=basename(output$getColumn('CountMatrix')))
  
  # Determine cell cycle phases. The training data is only available for Hsap and Mmus Ensembl
  trainData = NULL
  if (startsWith(param$refBuild, "Homo_sapiens/Ensembl")) {
    trainData = readRDS(system.file("exdata", "human_cycle_markers.rds", package = "scran", mustWork=TRUE))
  } else if (startsWith(param$refBuild, "Mus_musculus/Ensembl")) {
    trainData = readRDS(system.file("exdata", "mouse_cycle_markers.rds", package = "scran", mustWork=TRUE))
  }
  if (!is.null(trainData)) {
    cellCycleData = scran::cyclone(counts, trainData)
    cellPhase = data.frame(Name = colnames(counts), Phase = cellCycleData$phases)
    write.table(cellPhase, file = basename(output$getColumn('CellCyclePhase')), quote = F, sep = "\t", row.names = F)
  }
  return("SUCCESS")
}

##' @template app-template
##' @templateVar method ezMethodSingleCellCounts(input=NA, output=NA, param=NA)
##' @description Use this reference class to run 
EzAppSingleCellCounts <-
  setRefClass("EzAppSingleCellCounts",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSingleCellCounts
                  name <<- "EzAppSingleCellCounts"
                  appDefaults <<- rbind(mapMethod=ezFrame(Type="character",	DefaultValue="STAR",	Description="the mapper to use"),
                                        mapOptions=ezFrame(Type="character", DefaultValue="", Description="options passed to the mapper"),
                                        writeIgvSessionLink=ezFrame(Type="logical", DefaultValue=FALSE, Description="whether to write IGV session links.")
                                        )
                }
              )
  )

### EzAppSingleCellSTAR
EzAppSingleCellSTAR <- 
  setRefClass("EzAppSingleCellSTAR",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodSingleCellSTAR
                  name <<- "ezMethodSingleCellSTAR"
                  appDefaults <<- rbind(getJunctions=ezFrame(Type="logical",  DefaultValue="FALSE",	Description="should junctions be returned"),
                                        writeIgvSessionLink=ezFrame(Type="logical", DefaultValue="TRUE", Description="should an IGV link be generated"),
                                        markDuplicates=ezFrame(Type="logical", DefaultValue="FALSE", Description="should duplicates be marked with picard"),
                                        checkStrandness=ezFrame(Type="logical", DefaultValue="TRUE", Description="should strandness be checked"),
                                        randomSleep=ezFrame(Type="logical",  DefaultValue="FALSE",  Description="should there be a random sleep to avoid to much network traffic when loading the STAR index"),
                                        twopassMode=ezFrame(Type="logical", DefaultValue="FALSE", Description="1-pass mapping or basic 2-pass mapping")
                  )
                }
              )
  )
### STAR for single cell data: reads in a unmapped bam
###
ezMethodSingleCellSTAR = function(input=NA, output=NA, param=NA){
  
  refDir = getSTARReference(param)
  bamFile = output$getColumn("BAM")
  if(!is.null(param$randomSleep)){
    if(param$randomSleep){
      randomNumber = runif(1, min = 0, max = 1)
      if(randomNumber <= 1/3) {
        cat('Wait 15m \n')
        Sys.sleep( 900) 
      } else if(randomNumber > 1/3 & randomNumber <= 2/3) {
        cat('Wait 30m \n')
        Sys.sleep( 1800)
      }
    }
  }
  
  if(input$readType() == "bam"){
    fastqInput <- ezMethodBam2Fastq(input = input, param = param)
    trimmedInput <- ezMethodTrim(input = fastqInput, param = param)
    file.remove(fastqInput$getColumn("Read1"))
    if(param$paired)
      file.remove(fastqInput$getColumn("Read2"))
  }else{
    trimmedInput <- ezMethodTrim(input = input, param = param)
  }
  if(param$cmdOptions == "")
    param$cmdOptions <- "--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --outSAMattributes All"
  
  if (!grepl("outSAMattributes", param$cmdOptions)){
    param$cmdOptions = paste(param$cmdOptions, "--outSAMattributes All")
  }
  cmd = paste("STAR", " --genomeDir", refDir,  "--sjdbOverhang 150", 
              "--readFilesIn", trimmedInput$getColumn("Read1"), 
              if(param$paired) trimmedInput$getColumn("Read2"),
              "--twopassMode", ifelse(param$twopassMode, "Basic", "None"),
              "--runThreadN", ezThreads(), param$cmdOptions, 
              "--outStd BAM_Unsorted --outSAMtype BAM Unsorted",
              ">  Aligned.out.bam")## writes the output file Aligned.out.bam
  ##"|", "samtools", "view -S -b -", " >", "Aligned.out.bam")
  ezSystem(cmd)
  file.remove(trimmedInput$getColumn("Read1"))
  if(param$paired)
    file.remove(trimmedInput$getColumn("Read2"))
  
  on.exit(file.remove(c("Log.progress.out", "Log.out", 
                        "Log.std.out")), add=TRUE) ## clean star log files
  
  
  ## Merge unmapped and mapped bam to recover the tags
  if(input$readType() == "bam"){
    mergeBamAlignments(alignedBamFn="Aligned.out.bam",
                       unmappedBamFn=input$getFullPaths("Read1"),
                       outputBamFn="Aligned.out.merged.bam",
                       fastaFn=param$ezRef@refFastaFile)
    file.remove("Aligned.out.bam")
    file.rename(from="Aligned.out.merged.bam", to="Aligned.out.bam")
  }
  
  nSortThreads = min(ezThreads(), 8)
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
    javaCall = paste0("java", " -Djava.io.tmpdir=. -Xmx", 
                      min(floor(param$ram), 10), "g")
    cmd = paste0(javaCall, " -jar ", "$Picard_jar", " MarkDuplicates ",
                 " TMP_DIR=. MAX_RECORDS_IN_RAM=2000000", " I=", "sorted.bam",
                 " O=", basename(bamFile),
                 " REMOVE_DUPLICATES=false", ## do not remove, do only mark
                 " ASSUME_SORTED=true",
                 " VALIDATION_STRINGENCY=SILENT",
                 " METRICS_FILE=" ,"dupmetrics.txt",
                 " VERBOSITY=WARNING",
                 " >markdup.stdout 2> markdup.stderr")
    on.exit(file.remove(c("markdup.stdout", "markdup.stderr")), add=TRUE)
    
    ezSystem(cmd)
    ezSystem(paste("samtools", "index", basename(bamFile)))
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
