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
  sampleBam <- list.files('.', pattern = 'sample_alignments.bam$', recursive=TRUE)
  if(length(sampleBam) == 1L) { #CellRanger Multi Output
    setwd(dirname(sampleBam))
    system('mv sample_alignments.bam possorted_genome_bam.bam')
    system('samtools index possorted_genome_bam.bam')
    system('mv sample_filtered_feature_bc_matrix filtered_feature_bc_matrix')
    system(sprintf('mv * %s', file.path(cwd, sampleDir, "outs")))
    setwd(cwd)
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
