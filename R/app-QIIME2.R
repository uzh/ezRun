###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodQIIME2 = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  
  require(rmarkdown)
  dataset = input$meta
  sampleNames = input$getNames() 
  isPaired <- param$paired
  ### read fastq files and prepare inputs
  file1PathInDataset <- input$getFullPaths("Read1")
  if(isPaired){
    file2PathInDataset <- input$getFullPaths("Read2")
  }
  
  manifestCmd <- paste("awk '{print $1'\t'$2}'", dataset, "| tail -n +2 | sed '1i sample-id\tabsolute-filepath' > manifest.tsv")
  if(isPaired){
    manifestCmd <- paste("awk '{print $1'\t'$2'\t'$4}'", dataset, "| tail -n +2 | sed '1i sample-id\tforward-absolute-filepath\treverse-absolute-filepath' > manifest.tsv")
  }
  ezSystem(manifestCmd)
  
  sampleMetadataCmd <- paste("awk '{print $1'\t'$3}'", dataset, "| tail -n +2 | sed '1i sample-id\tGroup' > sample_metadata.tsv")
  ezSystem(sampleMetadataCmd)
  
  updateBatchCmd1 <- paste0("sed -e s/\"TRIM_LEFT\"/", param$trim_left, "/g",
                                 " -e s/\"TRUNC_LEN\"/", param$truncate_len, "/g",
                                 " -e s/\"SAMPLING_DEPTH\"/", param$sampling_depth, "/g ",
                                 file.path(METAGENOMICS_ROOT,UNIFIED_QIIME2_WORKFLOW_SINGLEEND), 
                                 " > ",
                                 UNIFIED_QIIME2_WORKFLOW_SINGLEEND)
  if(isPaired){
    updateBatchCmd1 <- paste0("sed -e s/\"TRIM_LEFT\"/", param$trim_left, "/g",
                                " -e s/\"TRUNC_LEN\"/", param$truncate_len, "/g",
                                " -e s/\"SAMPLING_DEPTH\"/", param$sampling_depth, "/g ",
                                file.path(METAGENOMICS_ROOT,UNIFIED_QIIME2_WORKFLOW_PAIREDEND), 
                                " > ",
                                UNIFIED_QIIME2_WORKFLOW_PAIREDEND)
  }
  ezSystem(updateBatchCmd1)
  
  cmdQIIME2 = paste("sh",UNIFIED_QIIME2_WORKFLOW_SINGLEEND)
  if(isPaired){
    cmdQIIME2 = paste("sh",UNIFIED_QIIME2_WORKFLOW_PAIREDEND)
  }
  ezSystem(cmdQIIME2)
  
  ZipQZVCmd <- paste("for i in $(ls *.qzv); do mv $i ${i}.zip; done")
  ezSystem(ZipQZVCmd)
  unzipQZVCmd <- paste("for i in $(ls *.zip); do unzip $i -d $i.folder; done")
  ezSystem(unzipQZVCmd)
  createDirCmd <- paste("echo -e './dada2_denoising_stats.qzv.zip.folder\n./dada2_rep_set.qzv.zip.folder\n./demux_seqs.qzv.zip.folder\n./table.qzv.zip.folder' > list_of_folders_of_interest")
  ezSystem(createDirCmd)
  moveOUTCmd <- paste("cat list_of_folders_of_interest | while read -r line; do find $line -name 'data' -exec cp -r {} $line/ \;; done")
  ezSystem(moveOUTCmd)
  
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodQIIME2()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppQIIME2 <-
  setRefClass("ezMethodQIIME2",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodQIIME2
                  name <<- "ezMethodQIIME2"
                  appDefaults <<- rbind(trim_left= ezFrame(Type="integer",  DefaultValue="0",Description="Position at which sequences should be trimmed due to low quality"),
                                        truncate_len= ezFrame(Type="integer",  DefaultValue="150",Description="Position at which sequences should be truncated due to decrease in quality"),
                                        sampling_depth= ezFrame(Type="integer",  DefaultValue="1000",Description="Total frequency that each sample should be rarefied to")
                  )
                }
              )
  )


