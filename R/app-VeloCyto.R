###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @template app-template
##' @templateVar method ezMethodVeloCyto(input=NA, output=NA, param=NA)
##' @description Use this reference class to run velocyto on CellRanger outputs
##' @author Lennart Opitz
EzAppVeloCyto <-
  setRefClass(
    "EzAppVeloCyto",
    contains = "EzApp",
    methods = list(
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodVeloCyto
        name <<- "EzAppVeloCyto"
        appDefaults <<- rbind(
          outputDir = ezFrame(
            Type = "character",
            DefaultValue = ".",
            Description = "Output directory"
          )
        )
      }
    )
  )


ezMethodVeloCyto <- function(input = NA, output = NA, param = NA) {
  if (
    ("SCDataOrigin" %in% input$colNames) &&
      input$getColumn("SCDataOrigin") == 'BDRhapsody'
  ) {
    result <- runVelocytoBD(input, output, param)
    return(result)
  }

  ###########
  ## assuming 10X
  gtfFile <- param$ezRef["refFeatureFile"]
  sampleName <- input$getNames()

  ##Copy data to scratch
  cellRangerPath <- file.path(input$dataRoot, input$getColumn("ResultDir"))
  cmd <- paste('cp -R', cellRangerPath, '.')
  ezSystem(cmd)

  sampleDir <- basename(cellRangerPath)
  ## restore the original outs directory
  cmd <- paste(
    'rsync -av --remove-source-files',
    paste0(sampleDir, '/*'),
    paste0(sampleDir, '/outs')
  )
  ezSystem(cmd)

  cwd <- getwd()
  sampleBam <- list.files(
    '.',
    pattern = 'sample_alignments.bam$',
    recursive = TRUE
  )
  sampleCram <- list.files(
    '.',
    pattern = 'sample_alignments.cram$',
    recursive = TRUE
  )
  sampleAlignPath <- c(sampleBam, sampleCram)

  if (length(sampleAlignPath) == 1L) {
    #CellRanger Multi Output
    sampleAlignFn <- basename(sampleAlignPath)
    fileExt <- tools::file_ext(sampleAlignFn)
    setwd(dirname(sampleAlignPath))
    system(sprintf('mv %s possorted_genome_bam.%s', sampleAlignFn, fileExt))
    system(sprintf('samtools index possorted_genome_bam.%s', fileExt))
    system('mv sample_filtered_feature_bc_matrix filtered_feature_bc_matrix')
    system(sprintf('mv * %s', file.path(cwd, sampleDir, "outs")))
    setwd(cwd)
  }

  convertCramToBam(
    file.path(sampleDir, 'outs', 'possorted_genome_bam.cram'),
    file.path(sampleDir, 'outs', 'possorted_genome_bam.bam'),
    cores = param$cores
  )

  # Run velocyto
  cmd <- paste('velocyto', 'run10x', sampleDir, gtfFile, '-@', param$cores)
  ezSystem(cmd)
  file.copy(file.path(sampleName, 'velocyto', paste0(sampleName, '.loom')), '.')
  ezSystem(paste('rm -Rf ', sampleName))
  return('Success')
}


runVelocytoBD <- function(input, output, param) {
  gtfFile <- param$ezRef["refFeatureFile"]
  sampleName <- input$getNames()

  ## convert the bam file
  gstoreBamFile <- input$getFullPaths("AlignmentFile")

  ## convert tags
  # samtools view -h /home/ubuntu/data/RNAVelo/Combined_Cartridge-1_Bioproduct_filtered_fixed3.bam \
  # | sed 's/MA:Z:/UB:Z:/' \
  # | samtools view -Sb -@6 -o /home/ubuntu/data/RNAVelo/Combined_Cartridge-1_Bioproduct_final.bam
  #
  localBamFile <- "bd.bam"
  cmd <- paste(
    "samtools view -h",
    gstoreBamFile,
    "| grep -v XF:Z:__intergenic", ## ignore all reads that don't align to genes
    "| grep -v  XF:Z:SampleTag", ## ignore sample tag reads
    "| sed s/MA:Z:/UB:Z:/ ",
    "| samtools view --tag UB -b -o",
    localBamFile
  )
  ezSystem(cmd)

  # Activate conda environment and run velocyto
  # velocyto run \
  # -b /home/ubuntu/data/RNAVelo/barcodes_C1.tsv \
  # -o /home/ubuntu/data/RNAVelo/ \
  # -m /home/ubuntu/data/RNAVelo/hg38_rmsk.gtf \
  # --samtools-threads 8 \
  # --samtools-memory 12000 \
  # /home/ubuntu/data/RNAVelo/Combined_Cartridge-1_Bioproduct_final.bam \
  # /home/ubuntu/data/RNAVelo/annotation.gtf
  #
  barcodesFile <- "barcodes.tsv"
  ezSystem(paste(
    "zcat",
    file.path(input$getFullPaths("CountMatrix"), "barcodes.tsv.gz"),
    ">",
    barcodesFile
  ))
  cmd <- paste(
    ". /usr/local/ngseq/miniforge3/etc/profile.d/conda.sh",
    "&& conda activate gi_velocyto",
    "&& velocyto run",
    "-b",
    barcodesFile,
    "-e",
    input$getNames(),
    "-o",
    ".",
    # optional the msk file "-m"
    "--samtools-memory",
    floor(param$ram * 0.7 / param$cores * 1000),

    '--samtools-threads',
    param$cores,
    localBamFile,
    gtfFile
  )
  ezSystem(cmd)
  return('Success')
}


##' @title Convert CRAM to BAM if needed
##' @description Checks if the input file is a CRAM file and converts it to BAM format using samtools.
##' Requires samtools to be available in the PATH.
##' @param inputFile path to the input file (can be BAM or CRAM)
##' @param outputBam desired output BAM file path
##' @param cores number of CPU cores to use (must be a positive integer)
##' @return Returns the path to the BAM file (either original or converted)
##' @details This function will stop with an error if:
##' \itemize{
##'   \item The input file does not exist
##'   \item The cores parameter is not a valid positive integer
##'   \item The BAM conversion fails
##'   \item The BAM indexing fails
##' }
convertCramToBam <- function(inputFile, outputBam, cores = 1) {
  # Check if the input file is a CRAM based on extension
  if (grepl("\\.cram$", inputFile, ignore.case = TRUE)) {
    ezLog("Detected CRAM file, converting to BAM format...")

    # Convert CRAM to BAM using samtools view
    # Note: -b flag outputs BAM format, -o specifies output file
    cmd <- paste(
      "samtools view -b -@",
      as.integer(cores),
      "-o",
      shQuote(outputBam),
      shQuote(inputFile)
    )
    ezSystem(cmd)

    # Verify the BAM file was created successfully
    if (!file.exists(outputBam)) {
      stop(paste0("Failed to create BAM file: ", outputBam))
    }

    # Index the resulting BAM file
    ezLog("Indexing BAM file...")
    cmd <- paste("samtools index", shQuote(outputBam))
    ezSystem(cmd)

    # Verify the index was created successfully
    if (!file.exists(paste0(outputBam, ".bai"))) {
      stop(paste0("Failed to create BAM index for: ", outputBam))
    }

    ezLog("CRAM to BAM conversion completed successfully")
    return(outputBam)
  } else {
    # Input is already BAM, return as-is
    return(inputFile)
  }
}
