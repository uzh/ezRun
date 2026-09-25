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
      ## velocyto (run10x for 10x, run for BD Rhapsody) and samtools unconditional.
      citation = function() {
        c(
          "La Manno, G. et al. RNA velocity of single cells. Nature 560, 494-498 (2018). https://doi.org/10.1038/s41586-018-0414-6",
          "Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25(16), 2078-2079 (2009). https://doi.org/10.1093/bioinformatics/btp352"
        )
      },
      methods_facts = function(param = list()) {
        c(
          ## ezMethodVeloCyto (app-VeloCyto.R:92); velocyto run10x --help and source checked in the gi_velocyto env (velocyto 0.17.17)
          "For 10x data, spliced and unspliced counts were obtained with velocyto run10x on the CellRanger output and the genes.gtf of refBuild, restricted to CellRanger's filtered cell barcodes, with samtools sorting on cores threads and velocyto's default 2048 MB per thread.",
          ## app-VeloCyto.R:92 and :141-159 (no -l, -M or -m given); velocyto logic.py Default = Permissive10X
          "Reads were assigned to spliced, unspliced and ambiguous molecules with velocyto's default logic (Permissive10X), counting only uniquely mapped reads (--multimap off) and without a repeat-mask annotation.",
          ## app-VeloCyto.R:62-84
          "For CellRanger Multi output, the per-sample alignment file and the per-sample filtered matrix were renamed to the cellranger count layout before velocyto was run, so the per-sample cell calls were used.",
          ## convertCramToBam (app-VeloCyto.R:179-216)
          "CRAM alignments were converted to BAM with samtools before counting.",
          ## app-VeloCyto.R:92 (no -t) vs :141-159; velocyto run10x/run --help
          "The loom layers were stored as uint16 for 10x data (run10x default) and uint32 for BD Rhapsody data (velocyto run default).",
          ## app-VeloCyto.R:92 and :141-159 build the command without param$cmdOptions; param$featureLevel is not read
          "The cmdOptions and featureLevel parameters were not passed to velocyto and had no effect.",
          ## runVelocytoBD (app-VeloCyto.R:113-122)
          "For BD Rhapsody data (SCDataOrigin = BDRhapsody), reads tagged XF:Z:__intergenic or XF:Z:SampleTag were removed from the BD alignment file, the BD molecule tag MA was renamed to UB, and only reads carrying a UB tag were kept.",
          ## runVelocytoBD (app-VeloCyto.R:134-159)
          "BD Rhapsody data were counted with velocyto run restricted to the cell barcodes of the BD filtered count matrix (CB tag), with samtools sorting on cores threads and 70% of the job memory divided across threads."
        )
      },
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

  cramPath <- file.path(sampleDir, 'outs', 'possorted_genome_bam.cram')
  bamPath <- file.path(sampleDir, 'outs', 'possorted_genome_bam.bam')
  alignFile <- if (file.exists(cramPath)) cramPath else bamPath
  convertCramToBam(alignFile, bamPath, cores = param$cores)

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
