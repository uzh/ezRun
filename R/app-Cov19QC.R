###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCov19QC <- function(input = NA, output = NA, param = NA, htmlFile = "00index.html"){
    if (ezIsSpecified(param$samples)){
        input <- input$subset(param$samples)
    }  
    setwdNew(param[['name']])
    ###1. Get Adapter Dimer Fraction (per Sample Mapping against AdapterSeq)
    adapterResult <- getAdapterStats(param, input)
    
    ###2. Preprocess with fastp
    param[['cmdOptionsFastp']] <- paste0('--reads_to_process=', param$readsUsed)
    inputProc <- ezMethodFastpTrim(input = input, param = param)
    
    ###3. Run Bowtie2 Mapping to get basic bamStats
    mappingResult <- mapToCovidGenome(param, input = inputProc, workDir="mappingResult")
    
    ###4. Report results (% Primer Dimer, % Mapped Reads, AvgCoverage, StdCoverage)
    dataset <- data.frame(Name = rownames(input$meta), input$meta, stringsAsFactors = FALSE, check.names = FALSE)
    result <- merge(adapterResult, mappingResult, by.x = 'Name', by.y = 'Name', all.x = TRUE, all.y = TRUE)
    result <- merge(result, dataset, by.x = 'Name', by.y = 'Name', all.x = TRUE, all.y = TRUE)
    ezWrite.table(result, 'result.txt', row.names = FALSE)
    
    file <- file.path(system.file("templates", package="ezRun"), 
                      "Cov19QC.Rmd")
    file.copy(from=file, to=basename(file), overwrite=TRUE)
    rmarkdown::render(input=basename(file), envir = new.env(),
                      output_dir=".", output_file=htmlFile,
                      quiet=TRUE)
    ezSystem('rm *.fastq.gz adapters.fa mappingResult/*.bai')
    return('success')
}

getAdapterStats <- function(param, input){
    inputFiles <- input$getFullPaths("Read1")
    if(param$Adapter1!='') {
            pattern <- param$Adapter1
        } else {
            pattern <- 'CTGTCTCTTATACACATCT'
        }
    result <- data.frame(Name = names(inputFiles), adapterDimerFraction = 0)
    i <- 0
    for (nm in names(inputFiles)) {
        i <- i + 1
        readCount <- min(c(20000, input$meta[['Read Count']][i]))
        file <- inputFiles[nm]
        myCmd <- paste('zcat', file, '|head -n 80000|grep', paste0('"^.*',pattern,'"'),'-o |while read LINE; do echo $LINE | wc -m; done')
        adapterPositions <- as.numeric(system(myCmd, intern = TRUE))-nchar(pattern)
        result[['adapterDimerFraction']][i] <- sum(adapterPositions <= 10, na.rm = TRUE)/readCount
    }
    return(result)
}

mapToCovidGenome <- function(param, input, workDir){
    R1_files = input$getFullPaths("Read1")
    R2_files = input$getFullPaths("Read2")
    
    result <- data.frame(Name = names(R1_files), mappingRate = 0, avgCov = 0, sdCov = 0)
    ref <- getBowtie2Reference(param)
    defOpt <- paste("-p", param$cores)
    setwdNew(workDir)
    i <- 0
    for (nm in names(R1_files)) {
        i <- i + 1
        nReads <- min(c(as.numeric(param$readsUsed), input$meta[['Read Count']][i]))
        bamFile <- paste0(nm, ".bam")
        cmd <- paste(
            "bowtie2", param$cmdOptions, defOpt,
            "-x", ref, if (param$paired) "-1", R1_files[nm],
            if (param$paired) paste("-2", R2_files[nm]),
            "2>", paste0(nm, "_bowtie2.log"), "|",
            "samtools", "view -S -b - > bowtie.bam")
        ezSystem(cmd)
        ezSortIndexBam("bowtie.bam", bamFile,
                       ram = param$ram, removeBam = TRUE,
                       cores = param$cores)
        
        mappingStats <- getBamMultiMatching(param, bamFile, nReads)
        result[['mappingRate']][i] <- 100* (sum(mappingStats[2:length(mappingStats)])/nReads)
        reads <- ezReadBamFileAsGRanges(bamFile, chromosomes = NULL, pairedEndReads = param$paired,
                                        max.fragment.width = 5000, min.mapq = 10, remove.duplicate.reads = FALSE)
        cov <- coverage(reads)
        avgCov <- mean(unlist(cov), na.rm = TRUE)
        sdCov <- sd(unlist(cov), na.rm = TRUE)
        result[['avgCov']][i] <- avgCov
        result[['sdCov']][i] <- sdCov
        
        file.remove(bamFile)
    }
    setwd('..')
    return(result)
}

##' @template app-template
##' @templateVar method ezMethodFastqScreen(input=NA, output=NA, param=NA, htmlFile="00index.html")
##' @description Use this reference class to run
##' @section Functions:

EzAppCov19QC <-
    setRefClass("EzAppCov19QC",
                contains = "EzApp",
                methods = list(
                    initialize = function() {
                        "Initializes the application using its specific defaults."
                        runMethod <<- ezMethodCov19QC
                        name <<- "EzAppCov19QC"
                        appDefaults <<- rbind(Adapter1 = ezFrame(Type = "character", DefaultValue = "CTGTCTCTTATACACATCT",
                                                                      Description = "Adapter Sequence for dimer calc"))
                    }
                )
    )