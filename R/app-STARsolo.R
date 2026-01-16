###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSTARsolo = 
    setRefClass("EzAppSTARsolo",
            contains = "EzApp",
            methods = list(
                initialize = function()
                {
                    "Initializes the application using its specific defaults."
                    runMethod <<- ezMethodSTARsolo
                    name <<- "EzAppSTARsolo"
                    appDefaults <<- rbind(controlSeqs=ezFrame(Type="charVector",
                                                              DefaultValue="",
                                                              Description="control sequences to add"),
                                          ## STARsolo parameters
                                          soloType=ezFrame(Type="character",
                                                           DefaultValue="CB_UMI_Simple",
                                                           Description="CB_UMI_Simple (a.k.a. Droplet), CB_UMI_Complex (e.g. Droplet)."),
                                          soloCBwhitelist=ezFrame(Type="character",
                                                           DefaultValue="SC3Pv3",
                                                           Description="Barcode whitelist. Choose: SC3Pv1, SC3Pv2, SC3Pv3."),
                                          soloCBstart=ezFrame(Type="character",
                                                           DefaultValue="1",
                                                           Description="Cell barcode start base."),
                                          soloCBlen=ezFrame(Type="character",
                                                           DefaultValue="16",
                                                           Description="Cell barcode length."),
                                          soloUMIstart=ezFrame(Type="character",
                                                           DefaultValue="17",
                                                           Description="UMI start base."),
                                          soloUMIlen=ezFrame(Type="character",
                                                             DefaultValue="auto",
                                                             Description="UMI length. 10 if SC3Pv3 not selected, 12 otherwise."),
                                          soloUMIfiltering=ezFrame(Type="character",
                                                                   DefaultValue="MultiGeneUMI",
                                                                   Description="type of UMI filtering.
                                                                            -               ... basic filtering: remove UMIs with N and homopolymers (similar to CellRanger 2.2.0)
                                                                            MultiGeneUMI    ... remove lower-count UMIs that map to more than one gene (introduced in CellRanger 3.x.x)"),
                                          soloCBmatchWLtype=ezFrame(Type="character",
                                                                    DefaultValue="1MM_multi_pseudocounts",
                                                                    Description="matching the Cell Barcodes to the WhiteList.
                                                                                Exact                   only exact matches allowed
                                                                                1MM                     only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact match.
                                                                                1MM_multi               multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches.
                                                                                                        Allowed CBs have to have at least one read with exact match. Similar to CellRanger 2.2.0
                                                                                1MM_multi_pseudocounts  same as 1MM_Multi, but pseudocounts of 1 are added to all whitelist barcodes.
                                                                                                        Similar to CellRanger 3.x.x"),
                                          keepAlignment=ezFrame(Type="logical",
                                                                 DefaultValue="TRUE",
                                                                 Description=""),
                                          soloCellFilter=ezFrame(Type="character",
                                                                 DefaultValue="EmptyDrops_CR",
                                                                 Description="cell filtering type and parameters
                                                                        EmptyDrops_CR   ... EmptyDrops filtering in CellRanger flavor. Please cite the original EmptyDrops paper: A.T.L Lun et al, Genome Biology, 20, 63 (2019): https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y
                                                                        CellRanger2.2   ... simple filtering of CellRanger 2.2, followed by thre numbers: number of expected cells, robust maximum percentile for UMI count, maximum to minimum ratio for UMI count
                                                                        TopCells        ... only report top cells by UMI count, followed by the excat number of cells
                                                                        None            ... do not output filtered cells")
                                         )
                }
            )
)

ezMethodSTARsolo = function(input=NA, output=NA, param=NA){
  ## check if a new index is needed, then create it with STAR.
  refDir <- getSTARSoloReference(param)
  
  sampleName = input$getNames()
  sampleDirs = strsplit(input$getColumn("RawDataDir"), ",")[[sampleName]]
  sampleDirs <- file.path(input$dataRoot, sampleDirs)
  if (all(grepl("\\.tar$", sampleDirs))) {
    runDirs <- c()
    for (i in 1:length(sampleDirs)) {
      runDirs[i] <- paste0("run_", i)
      res <- untar(sampleDirs[i], exdir=runDirs[i], tar=system("which tar", intern=TRUE))
      if (res > 0) {
        stop(sprintf("There was an error unpacking '%s' into '%s'. See warnings.",
                     sampleDirs[i], runDirs[i]))
      }
    }
  }
  # parse solo features
  soloFeatures <- unlist(strsplit(param$soloFeatures, ","))
  if (length(soloFeatures) == 0) {
    stop("No Solo Features specified!")
  }
  
  # create STARsolo command
  cmd = makeSTARsoloCmd(param, refDir, sampleName, runDirs, soloFeatures)
  
  # run STARsolo shell command
  ezSystem(cmd)
  
  # filter raw counts with EmptyDrop
  require(Matrix)
  outputDirs = file.path(sampleName,'Solo.out',soloFeatures)
  names(outputDirs) <- soloFeatures
  ## Remove alignments file if not setted differently
  samFn = paste0(sampleName,'/Aligned.out.sam')
  if(param[['keepAlignment']]){
    ezSystem(paste0("mv ", sampleName, "/Aligned.sortedByCoord.out.bam", " ",
                    sampleName, "/possorted_genome_bam.bam"))
    ezSystem(paste0("samtools index ", sampleName, "/possorted_genome_bam.bam"))
  }else{
    ezSystem(paste0("rm ", sampleName, "/Aligned.sortedByCoord.out.bam"))
  }
  ### gzip all files
  for (outputDir in outputDirs) {
    for (subDir in c("raw", "filtered")) {
      ezSystem(paste("pigz -p 4",file.path(outputDir, subDir, "*")))
    }
  }

  ### save in new directory
  rawDir <- paste0(sampleName,'/raw_feature_bc_matrix')
  filteredDir = paste0(sampleName,'/filtered_feature_bc_matrix')
  # We have to make subdirectories for every solo feature
  dir.create(rawDir)
  dir.create(filteredDir)
  
  # Process rest of features
  for (soloFeatureParam in soloFeatures) {
    outputDir <- outputDirs[soloFeatureParam]
    subRawDir <- file.path(rawDir, soloFeatureParam)
    subFilteredDir <- file.path(filteredDir, soloFeatureParam)
    if (soloFeatureParam == "Velocyto") {
      writeVelocyto(outputDir, subRawDir, subFilteredDir)
    } else {
      writeGenericFeature(outputDir, subRawDir, subFilteredDir)
    }
  }
  unlink(file.path(sampleName,'Solo.out'), recursive = TRUE, force=TRUE, expand=FALSE)
  return("Success")
}

# create STARsolo shell command
makeSTARsoloCmd = function(param, refDir, sampleName, sampleDirs, soloFeatures){
  require('tools')

  ## decide which chemistry whitelist to take
  barcodeInclusionListPath <- getBarcodeInclusionListPath(param)
  stopifnot("Could not find barcode whitelist" = file.exists(barcodeInclusionListPath))
  
  ## decide soloUMIlen
  if(param[['soloUMIlen']]=='auto') {
    if (param[['soloCBwhitelist']] %in% c('SC3Pv3', 'SC3Pv4')) {
      soloUMIlen = '12'
    } else{
      soloUMIlen = '10'
    }
  } else{
    soloUMIlen = param[['soloUMIlen']]
  }
  
  ## format soloFeatures list
  soloFeatures <- paste(soloFeatures, collapse=" ")
  
  ## create readFilesIn
  ### list files: the names are always the same for standard runs. Only the presence of indexes is optional.
  fastqfiles = sort(list.files(sampleDirs,full.names = TRUE,pattern = '.fastq', recursive = TRUE))
  cDNAfastqs = paste(grep('_R2',fastqfiles, value = TRUE), collapse = ',')
  barcodesfastqs = paste(grep('_R1',fastqfiles, value = TRUE), collapse = ',')
  readFilesIn = paste(cDNAfastqs, barcodesfastqs, collapse = ' ')
  
  ## create readFilesCommand
  if (unique(file_ext(fastqfiles)) == 'gz') {
    readFilesCommand = 'zcat'
  } else{
    readFilesCommand = 'cat'
  }
  
  # create output folder for the sample (in CellRanger this was automated)
  if (!dir.exists(sampleName)) {
    dir.create(sampleName)
  }
  
  # create full STARsolo command
  cmd = paste(
    "STAR",
    ## STAR general parameters
    paste0('--runThreadN ', param[['cores']]),
    paste0('--outFileNamePrefix ', sampleName, '/'),
    
    paste0('--readFilesIn ', readFilesIn),
    paste0('--readFilesCommand ', readFilesCommand),
    paste0('--genomeDir ', refDir),
    
    ## STARsolo parameters
    paste0('--soloType ', param[['soloType']]),
    paste0('--soloCBwhitelist ', barcodeInclusionListPath),
    paste0('--soloCBstart ', param[['soloCBstart']]),
    paste0('--soloCBlen ', param[['soloCBlen']]),
    paste0('--soloUMIstart ', param[['soloUMIstart']]),
    paste0('--soloUMIlen ', soloUMIlen),
    
    paste0('--soloUMIfiltering ', param[['soloUMIfiltering']]),
    paste0('--soloCBmatchWLtype ', param[['soloCBmatchWLtype']]),
    paste0('--soloCellFilter ', param[['soloCellFilter']]),
    paste0('--soloFeatures ', soloFeatures),
    "--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM",
    "--outSAMtype BAM SortedByCoordinate",
    "--outBAMcompression 6",
    "--limitBAMsortRAM",
    format(param$ram * 0.4 * 1e9, scientific = FALSE) ## use only 80% of the available RAM
  )
  
  if (ezIsSpecified(param$cmdOptions)) {
    cmd = paste(cmd, param$cmdOptions)
  }
  
  return(cmd)
}

checkGzippedAndUnzip <- function(filePath) {
  filePathGz <- paste0(filePath, ".gz")
  if (file.exists(filePathGz)) {
    # copy the file locally and unzip it
    tempFilePath <- basename(filePath)
    ezSystem(paste("gunzip -c", shQuote(filePathGz), ">", shQuote(tempFilePath)))
    return(tempFilePath)
  }
  return("")
}

getBarcodeInclusionListPath <- function(param) {
  soloCBwhitelist = list(
    SC3Pv1 = paste0(Sys.getenv("CellRanger"), '/lib/python/cellranger/barcodes/737K-april-2014_rc.txt'),
    SC3Pv2 = paste0(Sys.getenv("CellRanger"), '/lib/python/cellranger/barcodes/737K-august-2016.txt'),
    SC3Pv3 = paste0(Sys.getenv("CellRanger"), '/lib/python/cellranger/barcodes/3M-february-2018.txt'),
    SC3Pv4 = paste0(Sys.getenv("CellRanger"), '/lib/python/cellranger/barcodes/3M-3pgex-may-2023_TRU.txt')
  )
  
  soloWhiteListPath <- soloCBwhitelist[[param[['soloCBwhitelist']]]]
  ezWrite(sprintf("Trying %s", soloWhiteListPath))
  if (file.exists(soloWhiteListPath)) {
    return(soloWhiteListPath)
  }
  tempPath <- checkGzippedAndUnzip(soloWhiteListPath)
  ezWrite(sprintf("Trying %s", tempPath))
  if (tempPath != "") {
    return(tempPath)
  }
  # in some versions of CellRanger, the whitelist file has a different suffix
  soloWhiteListPathVariation <- sprintf("%s_TRU.txt", file_path_sans_ext(soloWhiteListPath))
  ezWrite(sprintf("Trying %s", soloWhiteListPathVariation))
  if (file.exists(soloWhiteListPathVariation)) {
    return(soloWhiteListPathVariation)
  }
  tempPath <- checkGzippedAndUnzip(soloWhiteListPathVariation)
  ezWrite(sprintf("Trying %s", tempPath))
  if (tempPath != "") {
    return(tempPath)
  }
  return("")
}

# Write Velocyto outputs
writeVelocyto <- function(outputDir, subRawDir, subFilteredDir) {
  subFeatures <- c("spliced", "unspliced", "ambiguous")
  # Create parent directory
  dir.create(subRawDir)
  dir.create(subFilteredDir)
  for (subFeature in subFeatures) {
    # Get root directory to raw and filtered counts
    rawOutputDir <- file.path(outputDir, "raw")
    filtOutputDir <- file.path(outputDir, "filtered")
    # Copy the appropriate files over
    writeVelocytoFeature(subFeature, rawOutputDir, subRawDir)
    writeVelocytoFeature(subFeature, filtOutputDir, subFilteredDir)
  }
}

# Write velocyto feature outputs
writeVelocytoFeature <- function(subFeature, fromDir, destDir) {
  # Make the destination sub-directory
  destSubDir <- file.path(destDir, subFeature)
  dir.create(destSubDir)
  
  # Copy the raw files over
  file.copy(from=Sys.glob(file.path(fromDir, "*.tsv.gz")), 
            to=destSubDir)
  file.copy(from=file.path(fromDir, paste0(subFeature, ".mtx.gz")), 
            to=file.path(destSubDir, "matrix.mtx.gz"))
}

# Write generic feature outputs
writeGenericFeature <- function(outputDir, subRawDir, subFilteredDir) {
  # Create parent directory
  dir.create(subRawDir)
  dir.create(subFilteredDir)
  # Get root directory to raw and filtered counts
  rawOutputDir <- file.path(outputDir, "raw")
  filtOutputDir <- file.path(outputDir, "filtered")
  # Copy files over
  file.copy(from=Sys.glob(file.path(rawOutputDir, "*.gz")), to=subRawDir)
  file.copy(from=Sys.glob(file.path(filtOutputDir, "*.gz")), to=subFilteredDir)
}

getSTARSoloReference <- function(param) {
  require(rtracklayer)
  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)
  
  if (ezIsSpecified(param$controlSeqs)) {
    refDir <- file.path(getwd(), "STARSolo_customised_Ref")
  } else {
    if (ezIsSpecified(param$transcriptTypes)) {
      starSoloBase <- paste(sort(param$transcriptTypes), collapse = "-")
      ## This is a combination of transcript types to use.
    } else {
      starSoloBase <- ""
    }
    refDir <- sub(
      "\\.gtf$", paste0("_STARSolo_", starSoloBase, "_Index"),
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
  if (file.exists(file.path(refDir, "SAindex"))) {
    ## we assume the index is built and complete
    return(refDir)
  }
  
  ## we have to build the reference
  setwd(dirname(refDir))
  ezWrite(Sys.info(), con = lockFile)
  on.exit(file.remove(lockFile), add = TRUE)
  
  job <- ezJobStart("STAR genome build")
  
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
  if (ezIsSpecified(param$controlSeqs)) {
    extraGR <- makeExtraControlSeqGR(param)
    gtfExtraFn <- tempfile(
      pattern = "extraSeqs", tmpdir = getwd(),
      fileext = ".gtf"
    )
    on.exit(file.remove(gtfExtraFn), add = TRUE)
    export.gff2(extraGR, con = gtfExtraFn)
    ezSystem(paste("cat", gtfExtraFn, ">>", gtfFile))
  }
  
  # STARSolo-specific things
  fai <- ezRead.table(paste0(param$ezRef["refFastaFile"], ".fai"), header = FALSE)
  colnames(fai) <- c("LENGTH", "OFFSET", "LINEBASES", "LINEWDITH")
  if (nrow(fai) > 50) {
    binOpt <- "--genomeChrBinNbits 16"
  } else {
    binOpt <- ""
  }
  
  genomeLength <- sum(fai$LENGTH)
  readLength <- 150 ## assumption
  indexNBasesOpt <- paste("--genomeSAindexNbases", min(13, floor(log2(genomeLength) / 2 - 1)))
  if (binOpt == "") {
    genomeChrBinNbits <- paste("--genomeChrBinNbits", floor(min(
      18,
      log2(max(genomeLength / nrow(fai), readLength))
    )))
  } else {
    genomeChrBinNbits <- ""
  }
  
  cmd <- paste(
    "STAR", "--runMode genomeGenerate --genomeDir", refDir, binOpt, indexNBasesOpt, genomeChrBinNbits,
    "--limitGenomeGenerateRAM", format(param$ram * 1e9, scientific = FALSE),
    "--genomeFastaFiles", genomeLocalFn,
    "--outTmpDir", tempfile(pattern="", fileext="_STARtmp", tmpdir=getwd()),
    "--sjdbGTFfile", gtfFile, "--sjdbOverhang 150", "--runThreadN", param$cores, "--genomeSAsparseD 2"
  )
  ezSystem(cmd)
  ezWriteElapsed(job, "done")
  file.remove(gtfFile)
  
  return(refDir)
}
