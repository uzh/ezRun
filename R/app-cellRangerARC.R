###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodCellRangerARC <- function(input = NA, output = NA, param = NA) {
  sampleName <- input$getNames()
  
  # Setup directories
  RNADataDir <- input$getColumn("RNADataDir") %>% 
    strsplit(",") %>% unlist() %>% sort() %>%  ## support multiple komma-separated entreies
    file.path(input$dataRoot, .)
  ATACDataDir <- input$getColumn("ATACDataDir") %>% 
    strsplit(",") %>% unlist() %>% sort() %>%  ## support multiple komma-separated entreies
    file.path(input$dataRoot, .)
  
  #1. extract tar files if they are in tar format
  ## RNA
  if (all(grepl("\\.tar$", RNADataDir))) {
    runRNADataDir <- tar2Fastq(RNADataDir, prefix="GEX")
    ## read folders have structure <prefix>--<run name>/<orig sample name>/<orig sample name>_<sampe lane pattern>*.fastq.gz
    ## if the "orig sample name> is different from the sampleName, then we rename
    sampleFqDirs <- list.dirs(runRNADataDir, full.names=TRUE, recursive=FALSE)
    idx <- which(basename(sampleFqDirs) != sampleName)
    cwd <- getwd()
    for (i in idx){
      myDir <- sampleFqDirs[i]
      setwd(myDir)
      cmd <- paste('rename',
                   paste0('s/', basename(myDir),'/',sampleName, '/g'),
                   paste0(basename(myDir),'*.gz'))
      ezSystem(cmd)
      setwd(cwd)
    }
  } else {
    stop("Require rna inputs to be provided in .tar files.")
  }
  
  ## ATAC
  if (all(grepl("\\.tar$", ATACDataDir))) {
    runATACDataDir <- tar2Fastq(ATACDataDir, prefix="ATAC")
    ## NOTE: renaming of ATAC samples is not supported
    if (ezIsSpecified(param$atacLengthR1)){
      fqFiles <- list.files(runATACDataDir, "_R1.*fastq.gz", full.names = TRUE)
      sapply(fqFiles, trimFastq, length=param$atacLengthR1)
    }
    if (ezIsSpecified(param$atacLengthR2)){
      fqFiles <- list.files(runATACDataDir, "_R3.*fastq.gz", full.names = TRUE)
      if (length(fqFiles) == 0){
        fqFiles <- list.files(runATACDataDir, "_R2.*fastq.gz", full.names = TRUE)
      }
      sapply(fqFiles, trimFastq, length=param$atacLengthR2)
    }
  } else {
    stop("Require atac inputs to be provided in .tar files.")
  }

  cellRangerARCFolder <- str_sub(sampleName, 1, 45)# %>% str_c("-cellRanger-arc")
  #3.Generate the cellranger command with the required arguments
  #3.1. Obtain the ARC reference
  refDir <- getCellRangerARCReference(param)
  
  #3.5. Create library file that contains the sample and atac dirs location
  libraryFn <- createARCLibraryFile(runRNADataDir, runATACDataDir)
  
  #3.6. Command
  cmd <- paste(
    "cellranger-arc count", 
    paste0("--id=", cellRangerARCFolder),
    paste0("--reference=", refDir),
    paste0("--libraries=", libraryFn),
    paste0("--localmem=", param$ram),
    paste0("--localcores=", param$cores),
    if (ezIsSpecified(param$expectedCells)) {paste0("--expect-cells=", param$expectedCells)},
    ifelse(ezIsSpecified(param$excludeIntrons) && param$excludeIntrons, "--gex-exclude-introns", "")
  )
  
  #4. Add additional cellranger-arc options if specified
  if (ezIsSpecified(param$cmdOptions)) {
    cmd <- paste(cmd, param$cmdOptions)
  }
  
  #5. Execute the command
  ezSystem(cmd)
  

  #7. Delete temp files and rename the final cellranger-arc output folder
  #unlink(dirname(RNADataDir), recursive = TRUE)
  if (exists("ATACDataDir")){
    #unlink(basename(ATACDataDir))
  }
  file.rename(file.path(cellRangerARCFolder, "outs"), basename(output$getColumn("ResultDir")))
  file.copy(libraryFn, basename(output$getColumn("ResultDir")))
  #unlink(cellRangerARCFolder, recursive = TRUE)
  if (ezIsSpecified(param$controlSeqs)) 
    unlink(refDir, recursive = TRUE)
  

   if(!param$keepBam){
     filesToRemove <- file.path(basename(output$getColumn("ResultDir")), 
                             c("gex_possorted_bam.bam", "gex_possorted_bam.bam.bai", 
                               "atac_possorted_bam.bam", "atac_possorted_bam.bam.bai"))
     unlink(filesToRemove)
   }
  
  return("Success")
}


createARCLibraryFile <- function(RNADataDir, ATACDataDir) {#}, sampleName, peakName) {
  libraryFn <- "ARC-library.csv"
  sampleName <- list.dirs(RNADataDir, recursive = FALSE) %>% basename() %>% unique()
  peakName <- list.dirs(ATACDataDir, recursive = FALSE) %>% basename() %>% unique()
  ## all tar files must have the same samplename
  stopifnot(length(sampleName) == 1)
  stopifnot(length(peakName) == 1)
  rnaLib <- tibble(fastqs=normalizePath(RNADataDir), sample=sampleName, library_type="Gene Expression")
  atacLib <- tibble(fastqs=normalizePath(ATACDataDir), sample=peakName, library_type="Chromatin Accessibility")
  write_csv(rbind(rnaLib, atacLib), libraryFn)
  return(libraryFn)
}




getCellRangerARCReference <- function(param) {
  require(rtracklayer)
  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)
  
  if (ezIsSpecified(param$controlSeqs) | ezIsSpecified(param$secondRef)) {
    refDir <- file.path(getwd(), "10X_customised_Ref")
  } else {
    if (ezIsSpecified(param$transcriptTypes)) {
      cellRangerBase <- paste(sort(param$transcriptTypes), collapse = "-")
      ## This is a combination of transcript types to use.
    } else {
      cellRangerBase <- ""
    }
    refDir <- sub(
      "\\.gtf$", paste0("_10XARC_SC_", cellRangerBase, "_Index"),
      param$ezRef["refFeatureFile"]
    )
  }
  
  lockFile <- paste0(refDir, ".lock")
  i <- 0
  while (file.exists(lockFile) && i < INDEX_BUILD_TIMEOUT) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste(
      "reference building still in progress after",
      INDEX_BUILD_TIMEOUT, "min"
    ))
  }
  ## there is no lock file
  if (file.exists(refDir)) {
    ## we assume the index is built and complete
    return(refDir)
  }
  
  ## we have to build the reference
  setwd(dirname(refDir))
  ezWrite(Sys.info(), con = lockFile)
  on.exit(file.remove(lockFile), add = TRUE)
  
  job <- ezJobStart("10X cellranger-arc build")
  
  if (ezIsSpecified(param$controlSeqs)) {
    ## make reference genome
    genomeLocalFn <- tempfile(
      pattern = "genome", tmpdir = getwd(),
      fileext = ".fa"
    )
    file.copy(from = param$ezRef@refFastaFile, to = genomeLocalFn)
    writeXStringSet(getControlSeqs(param$controlSeqs),
                    filepath = genomeLocalFn,
                    append = TRUE
    )
    on.exit(file.remove(genomeLocalFn), add = TRUE)
  } else if(ezIsSpecified(param$secondRef)){
    ## make reference genome
    genomeLocalFn <- tempfile(
      pattern = "genome", tmpdir = getwd(),
      fileext = ".fa"
    )
    file.copy(from = param$ezRef@refFastaFile, to = genomeLocalFn)
    secondaryRef  <- readDNAStringSet(param$secondRef)
    writeXStringSet(secondaryRef,
                    filepath = genomeLocalFn,
                    append = TRUE
    )
    on.exit(file.remove(genomeLocalFn), add = TRUE) 
  } else {
    genomeLocalFn <- param$ezRef@refFastaFile
  }
  
  ## make gtf
  gtfFile <- tempfile(
    pattern = "genes", tmpdir = getwd(),
    fileext = ".gtf"
  )
  if (ezIsSpecified(param$transcriptTypes)) {
    export.gff2(gtfByTxTypes(param, param$transcriptTypes),
                con = gtfFile
    )
  } else {
    file.copy(from = param$ezRef@refFeatureFile, to = gtfFile)
  }
  if (ezIsSpecified(param$controlSeqs)|ezIsSpecified(param$secondRef)) {
    extraGR <- makeExtraControlSeqGR(param)
    gtfExtraFn <- tempfile(
      pattern = "extraSeqs", tmpdir = getwd(),
      fileext = ".gtf"
    )
    on.exit(file.remove(gtfExtraFn), add = TRUE)
    export.gff2(extraGR, con = gtfExtraFn)
    ezSystem(paste("cat", gtfExtraFn, ">>", gtfFile))
  }
  
  cmd <- paste(
    "cellranger-arc mkref",
    paste0("--genome=", basename(refDir)),
    paste0("--fasta=", genomeLocalFn),
    paste0("--genes=", gtfFile),
    paste0("--nthreads=", param$cores)
  )
  ezSystem(cmd)
  
  file.remove(gtfFile)
  
  return(refDir)
}


##' @author Paul, Gueguen
##' @template app-template
##' @templateVar method ezMethodCellRanger(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppCellRangerARC <-
  setRefClass("EzAppCellRangerARC",
              contains = "EzApp",
              methods = list(
                initialize = function() {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodCellRangerARC
                  name <<- "EzAppCellRangerARC"
                  appDefaults <<- rbind(
                    chemistry = ezFrame(
                      Type = "character",
                      DefaultValue = "auto",
                      Description = "Assay configuration."
                    ),
                    expectedCells = ezFrame(
                      Type = "numeric",
                      DefaultValue = 10000,
                      Description = "Expected number of cells."
                    ),
                    excludeIntrons = ezFrame(
                      Type = "logical",
                      DefaultValue = TRUE,
                      Description = "Count reads on introns."
                    ),
                    controlSeqs = ezFrame(
                      Type = "charVector",
                      DefaultValue = "",
                      Description = "control sequences to add"
                    ),
                    bamStats = ezFrame(
                      Type = "logical",
                      DefaultValue = TRUE,
                      Description = "compute per cell alignment stats"
                    ),
                    runVeloCyto = ezFrame(
                      Type = "logical",
                      DefaultValue = FALSE,
                      Description = "run velocyto and generate loom file"
                    ),
                    keepBam = ezFrame(
                      Type = "logical",
                      DefaultValue = TRUE,
                      Description = "keep bam file produced by CellRangerARC"
                    )
                  )
                }
              )
  )
