###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodQIIME2 = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  
  require(rmarkdown)
  require(data.table)
  require(here)
  dataset = input$meta
  sampleNames = input$getNames() 
  isPaired <- param$paired
  ### read fastq files and prepare inputs
  file1PathInDataset <- input$getFullPaths("Read1")
  if(isPaired){
    file2PathInDataset <- input$getFullPaths("Read2")
  }
  
  if (!file.exists("sample_metadata.tsv")) {
    file.create("sample_metadata.tsv")
    sample_metadata <- dataset[, 2, drop = FALSE]
    sample_metadata$sample_id <- rownames(sample_metadata)
    sample_metadata <- sample_metadata[,c(2,1)]
    setnames(sample_metadata, "sample_id", "sample-id")
    setnames(sample_metadata, colnames(sample_metadata)[2], "Group")
    write_tsv(sample_metadata, file = "sample_metadata.tsv")
  } else {
    print("The file exists")
  }
  
  if (!file.exists("manifest.tsv")) {
  
    file.create("manifest.tsv")
    manifest <- dataset[, 1, drop = FALSE]
    manifest$sample_id <- rownames(manifest)
    manifest <- manifest[,c(2,1)]
    setnames(manifest, "sample_id", "sample-id")
    setnames(manifest, colnames(manifest)[2], "absolute-filepath")
    write_tsv(manifest, file = "manifest.tsv")
    if(isPaired){
      manifest <- table_metadata[, c(1,3)]
      manifest$sample_id <- rownames(manifest)
      manifest <- manifest[,c(3,1,2)]
      setnames(manifest, "sample_id", "sample-id")
      setnames(manifest, colnames(manifest)[2], "forward-absolute-filepath")
      setnames(manifest, colnames(manifest)[3], "reverse-absolute-filepath")
      write_tsv(manifest, file = "manifest.tsv")
    }
  } else {
    print("The file exists")
  }
  
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


