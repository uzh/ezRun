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
  setRefClass("EzAppVeloCyto",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
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

##' @title Convert CRAM to BAM if needed
##' @description Checks if the input file is a CRAM file and converts it to BAM format
##' @param inputFile path to the input file (can be BAM or CRAM)
##' @param outputBam desired output BAM file path
##' @param cores number of CPU cores to use
##' @return Returns the path to the BAM file (either original or converted)
convertCramToBamIfNeeded <- function(inputFile, outputBam, cores = 1) {
  # Validate inputs
  if (!file.exists(inputFile)) {
    stop(paste0("Input file does not exist: ", inputFile))
  }
  if (!is.numeric(cores) || cores < 1) {
    stop(paste0("Invalid cores parameter: ", cores))
  }
  
  # Check if the input file is a CRAM based on extension
  if (grepl("\\.cram$", inputFile, ignore.case = TRUE)) {
    ezWrite("Detected CRAM file, converting to BAM format...")
    
    # Convert CRAM to BAM using samtools view
    # Note: -b flag outputs BAM format, -o specifies output file
    cmd <- paste("samtools view -b -@", cores, "-o", shQuote(outputBam), shQuote(inputFile))
    ezSystem(cmd)
    
    # Verify the BAM file was created successfully
    if (!file.exists(outputBam)) {
      stop(paste0("Failed to create BAM file: ", outputBam))
    }
    
    # Index the resulting BAM file
    ezWrite("Indexing BAM file...")
    cmd <- paste("samtools index", shQuote(outputBam))
    ezSystem(cmd)
    
    # Verify the index was created successfully
    if (!file.exists(paste0(outputBam, ".bai"))) {
      stop(paste0("Failed to create BAM index for: ", outputBam))
    }
    
    ezWrite("CRAM to BAM conversion completed successfully")
    return(outputBam)
  } else {
    # Input is already BAM, return as-is
    return(inputFile)
  }
}


ezMethodVeloCyto <- function(input=NA, output=NA, param=NA){
  
  if (("SCDataOrigin" %in% input$colNames) && 
      input$getColumn("SCDataOrigin") == 'BDRhapsody'){
    result <- runVelocytoBD(input, output, param) 
    return(result)
  }
  
  ###########
  ## assuming 10X
  gtfFile <- param$ezRef["refFeatureFile"]
  sampleName <- input$getNames()
  
  ##Copy data to scratch
  cellRangerPath <- file.path(input$dataRoot,input$getColumn("ResultDir"))
  cmd <- paste('cp -R',  cellRangerPath, '.')
  ezSystem(cmd)
  
  sampleDir <- basename(cellRangerPath)
  ## restore the original outs directory
  cmd <- paste('rsync -av --remove-source-files', paste0(sampleDir,'/*'), paste0(sampleDir,'/outs'))
  ezSystem(cmd)
  
  cwd <- getwd()
  
  # Check for CRAM files and convert to BAM if needed
  sampleCram <- list.files('.', pattern = 'sample_alignments.cram$', recursive=TRUE)
  if(length(sampleCram) == 1L) { #CellRanger Multi Output with CRAM
    ezWrite("Found sample_alignments.cram, converting to BAM...")
    cramDir <- dirname(sampleCram)
    setwd(cramDir)
    convertCramToBamIfNeeded('sample_alignments.cram', 'sample_alignments.bam', param$cores)
    setwd(cwd)
  }
  
  sampleBam <- list.files('.', pattern = 'sample_alignments.bam$', recursive=TRUE)
  if(length(sampleBam) == 1L) { #CellRanger Multi Output
    setwd(dirname(sampleBam))
    system('mv sample_alignments.bam possorted_genome_bam.bam')
    system('samtools index possorted_genome_bam.bam')
    system('mv sample_filtered_feature_bc_matrix filtered_feature_bc_matrix')
    system(sprintf('mv * %s', file.path(cwd, sampleDir, "outs")))
    setwd(cwd)
  }
  
  # Check for possorted_genome_bam.cram and convert if needed
  possortedCram <- list.files('.', pattern = 'possorted_genome_bam.cram$', recursive=TRUE)
  if(length(possortedCram) > 0) {
    ezWrite("Found possorted_genome_bam.cram file(s), converting to BAM...")
    for(cramFile in possortedCram) {
      cramDir <- dirname(cramFile)
      setwd(cramDir)
      outputBam <- sub("\\.cram$", ".bam", basename(cramFile))
      convertCramToBamIfNeeded(basename(cramFile), outputBam, param$cores)
      # Remove the original CRAM file only if conversion was successful
      if(file.exists(outputBam) && file.exists(paste0(outputBam, ".bai"))) {
        tryCatch({
          file.remove(basename(cramFile))
          ezWrite(paste0("Removed original CRAM file: ", basename(cramFile)))
          if(file.exists(paste0(basename(cramFile), ".crai"))) {
            file.remove(paste0(basename(cramFile), ".crai"))
            ezWrite(paste0("Removed CRAM index file: ", basename(cramFile), ".crai"))
          }
        }, error = function(e) {
          ezWrite(paste0("Warning: Failed to remove CRAM file(s): ", e$message))
        })
      }
      setwd(cwd)
    }
  }

  # Activate conda environment and run velocyto
  conda_activate <- paste("source", "/usr/local/ngseq/miniforge3/etc/profile.d/conda.sh", "&&", "conda activate gi_velocyto", "&&")
  cmd <- paste("bash -c \"", conda_activate, 'velocyto run10x', sampleDir, gtfFile, '-@', param$cores, "\"")
  ezSystem(cmd)
  file.copy(file.path(sampleName, 'velocyto', paste0(sampleName,'.loom')), '.')
  ezSystem(paste('rm -Rf ', sampleName))
  return('Success')
}



runVelocytoBD <- function(input, output, param){
  
  gtfFile <- param$ezRef["refFeatureFile"]
  sampleName <- input$getNames()
  
  ## convert the bam file
  gstoreBamFile <- input$getFullPaths("AlignmentFile")
  
  ## Check if input is CRAM and convert to BAM if needed
  if (grepl("\\.cram$", gstoreBamFile, ignore.case = TRUE)) {
    ezWrite("Input file is CRAM, converting to BAM format...")
    convertedBam <- paste0("input_converted_", Sys.getpid(), ".bam")
    gstoreBamFile <- convertCramToBamIfNeeded(gstoreBamFile, convertedBam, param$cores)
  }
  
  ## convert tags
  # samtools view -h /home/ubuntu/data/RNAVelo/Combined_Cartridge-1_Bioproduct_filtered_fixed3.bam \
  # | sed 's/MA:Z:/UB:Z:/' \
  # | samtools view -Sb -@6 -o /home/ubuntu/data/RNAVelo/Combined_Cartridge-1_Bioproduct_final.bam
  # 
  localBamFile <- "bd.bam"
  cmd <- paste("samtools view -h", gstoreBamFile, 
               "| grep -v XF:Z:__intergenic", ## ignore all reads that don't align to genes
               "| grep -v  XF:Z:SampleTag",  ## ignore sample tag reads
               "| sed s/MA:Z:/UB:Z:/ ", 
               "| samtools view --tag UB -b -o", localBamFile)
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
  ezSystem(paste("zcat", file.path(input$getFullPaths("CountMatrix"), "barcodes.tsv.gz"), ">", barcodesFile))
  cmd <- paste(". /usr/local/ngseq/miniforge3/etc/profile.d/conda.sh",
               "&& conda activate gi_velocyto",
               "&& velocyto run",
               "-b", barcodesFile, 
               "-e", input$getNames(),
               "-o", ".",
               # optional the msk file "-m"
               "--samtools-memory", floor(param$ram * 0.7 / param$cores * 1000),
               
               '--samtools-threads', param$cores, localBamFile, gtfFile)
  ezSystem(cmd)
  return('Success')
}
