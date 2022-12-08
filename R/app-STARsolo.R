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
                                                                 DefaultValue="None",
                                                                 Description="cell filtering type and parameters
                                                                        CellRanger2.2   ... simple filtering of CellRanger 2.2, followed by thre numbers: number of expected cells, robust maximum percentile for UMI count, maximum to minimum ratio for UMI count
                                                                        TopCells        ... only report top cells by UMI count, followed by the excat number of cells
                                                                        None            ... do not output filtered cells")
                                         )
                }
            )
)

ezMethodSTARsolo = function(input=NA, output=NA, param=NA){
  # analysis vars
  ## check if a new index is needed, then create it with CellRanger and the last version of STAR.
  refDir <- getSTARReference(param)
  
  sampleName = input$getNames()
  sampleDirs = strsplit(input$getColumn("RawDataDir"), ",")[[sampleName]]
  sampleDirs <- file.path(input$dataRoot, sampleDirs)
  if (all(grepl("\\.tar$", sampleDirs))) {
    runDirs <- c()
    for (i in 1:length(sampleDirs)) {
      runDirs[i] <- paste0("run_", i)
      untar(sampleDirs[i], exdir=runDirs[i], tar=system("which tar", intern=TRUE))
    }
  }

  
  # create STARsolo command
  cmd = makeSTARsoloCmd(param, refDir, sampleName, runDirs)
  
  # run STARsolo shell command
  ezSystem(cmd)
  
  # filter raw counts with EmptyDrop
  require(Matrix)
  require(DropletUtils)
  if (grepl('--soloFeatures', param$cmdOptions)) {
    soloParams <- limma::strsplit2(param$cmdOptions, '--')
    soloFeatureParam <- sub('soloFeatures ', '', soloParams[grep('soloFeatures',soloParams)])
    soloFeatureParams <- stringr::str_split(soloFeatureParam, " ")[[1]]
    # Make sure "Genes" is at the beginning so we detect drops based on it
    if ("Gene" %in% soloFeatureParams) {
      soloFeatureParams <- c("Gene", soloFeatureParams[soloFeatureParams!="Gene"])
    }
  } else {
    soloFeatureParams <- "Gene"
  }
  outputDirs = file.path(sampleName,'Solo.out',soloFeatureParams,'raw')
  names(outputDirs) <- soloFeatureParams
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
    ezSystem(paste("pigz -p 4",paste0(outputDir,"/*")))
  }

  outputDir <- outputDirs[1]
  ## Filter raw STARsolo counts with EmptyDrops
  sce = read10xCounts(samples = outputDir, col.names = TRUE, version = '3')
  rawCounts = sce@assays@data@listData$counts
  ### cell calling
  called = defaultDrops(rawCounts)
  ### save in new directory
  rawDir <- paste0(sampleName,'/raw_feature_bc_matrix')
  filteredDir = paste0(sampleName,'/filtered_feature_bc_matrix')
  if (length(outputDirs) > 1) {
    # We have to make subdirectories for every solo feature
    dir.create(rawDir)
    dir.create(filteredDir)
    
    # First write the main results
    subRawDir <- file.path(rawDir, soloFeatureParams[1])
    subFilteredDir <- file.path(filteredDir, soloFeatureParams[1])
    
    # Write the raw and filtered counts of first feature
    writeRawAndFiltered10XCounts(sce, subRawDir, subFilteredDir, called)
    
    # Process rest of features
    for (soloFeatureParam in soloFeatureParams[-1]) {
      outputDir <- outputDirs[soloFeatureParam]
      subRawDir <- file.path(rawDir, soloFeatureParam)
      subFilteredDir <- file.path(filteredDir, soloFeatureParam)
      if (soloFeatureParam == "Velocyto") {
        writeVelocyto(outputDir, subRawDir, subFilteredDir, called)
      } else {
        writeGenericMode(outputDir, subRawDir, subFilteredDir, called)
      }
    }
  } else {
    writeRawAndFiltered10XCounts(sce, rawDir, filteredDir, called)
  }
  return("Success")
}
# create STARsolo shell command
makeSTARsoloCmd = function(param, refDir, sampleName, sampleDirs){
  require('tools')

  ## decide which chemistry whitelist to take
  soloCBwhitelist = list(
    SC3Pv1 = paste0(Sys.getenv("CellRanger"), '/lib/python/cellranger/barcodes/737K-april-2014_rc.txt'),
    SC3Pv2 = paste0(Sys.getenv("CellRanger"), '/lib/python/cellranger/barcodes/737K-august-2016.txt'),
    SC3Pv3 = paste0(Sys.getenv("CellRanger"), '/lib/python/cellranger/barcodes/3M-february-2018.txt')
  )
  
  ## decide soloUMIlen
  if(param[['soloUMIlen']]=='auto') {
    if (param[['soloCBwhitelist']] == 'SC3Pv3') {
      soloUMIlen = '12'
    } else{
      soloUMIlen = '10'
    }
  } else{
    soloUMIlen = param[['soloUMIlen']]
  }
  
  ## create readFilesIn
  ### list files: the names are always the same for standard runs. Only the presence of indexes is optional.
  fastqfiles = sort(list.files(sampleDirs,full.names = TRUE,pattern = '.fastq', recursive = TRUE))
  cDNAfastqs = paste(grep('R2',fastqfiles, value = TRUE), collapse = ',')
  barcodesfastqs = paste(grep('R1',fastqfiles, value = TRUE), collapse = ',')
  readFilesIn = trimws(paste(cDNAfastqs, barcodesfastqs, collapse = ' ')) # remove last space, incase indexfastqs is empty.
  
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
    paste0('--soloCBwhitelist ', soloCBwhitelist[[param[['soloCBwhitelist']]]]),
    paste0('--soloCBstart ', param[['soloCBstart']]),
    paste0('--soloCBlen ', param[['soloCBlen']]),
    paste0('--soloUMIstart ', param[['soloUMIstart']]),
    paste0('--soloUMIlen ', soloUMIlen),
    
    paste0('--soloUMIfiltering ', param[['soloUMIfiltering']]),
    paste0('--soloCBmatchWLtype ', param[['soloCBmatchWLtype']]),
    paste0('--soloCellFilter ', param[['soloCellFilter']]),
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

# Write Velocyto outputs
writeVelocyto <- function(outputDir, subRawDir, subFilteredDir, calledCells) {
  subFeatures <- c("spliced", "unspliced", "ambiguous")
  # Create parent directory
  dir.create(subRawDir)
  dir.create(subFilteredDir)
  for (subFeature in subFeatures) {
    outputDirSplit <- file.path(outputDir, "..", subFeature)
    # We create the directory manually
    dir.create(outputDirSplit)
    # Copy the files over
    file.copy(from=Sys.glob(file.path(outputDir, "*.tsv.gz")), to=outputDirSplit)
    file.copy(from=file.path(outputDir, paste0(subFeature, ".mtx.gz")), 
              to=file.path(outputDirSplit, "matrix.mtx.gz"))
    # Call the generic mode
    subRawFeatureDir <- file.path(subRawDir, subFeature)
    subFilteredFeatureDir <- file.path(subFilteredDir, subFeature)
    writeGenericMode(outputDirSplit, subRawFeatureDir, subFilteredFeatureDir, calledCells)
  }
}

# Write outputs from other solo features
writeGenericMode <- function(outputDir, rawDir, filteredDir, calledCells) {
  sce = read10xCounts(samples = outputDir, col.names = TRUE, version = '3')
  writeRawAndFiltered10XCounts(sce, rawDir, filteredDir, calledCells)
}

# Takes in a raw-count SCE and outputs the 
writeRawAndFiltered10XCounts <- function(rawSce, rawDir, filteredDir, calledCells) {
  ## Filter raw STARsolo counts with previous called
  rawCounts = rawSce@assays@data@listData$counts
  filteredCounts = rawCounts[,calledCells]
  geneID = rawSce@rowRanges@elementMetadata@listData[["ID"]]
  geneSymbol = rawSce@rowRanges@elementMetadata@listData[["Symbol"]]
  
  write10xCounts(path = rawDir, x = rawCounts, 
                 gene.id = geneID, gene.symbol = geneSymbol, version = '3')
  write10xCounts(path = filteredDir, x = filteredCounts, 
                 gene.id = geneID, gene.symbol = geneSymbol, version = '3')
}
