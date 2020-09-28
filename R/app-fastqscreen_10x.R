###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodFastqScreen_10x <- function(input=NA, output=NA, param=NA,
                                    htmlFile="00index.html"){
  require(tidyverse)
  dataset <- input$meta
  sampleDirs <- input$getFullPaths("RawDataDir")
  stopifnot(all(grepl("\\.tar$", sampleDirs)))
  
  taredfiles <- lapply(sampleDirs, untar, list=TRUE)
  if(any(str_detect(unlist(taredfiles), "_R3_"))){
    ## ATAC has R3 for real data
    taredfiles_R2 <- sapply(taredfiles,
                            function(x){grep("_R3_", x, value=TRUE) %>% head(1)})
  }else if(any(str_detect(unlist(taredfiles), "_R2_"))){
    ## RNA has R2 for real data
    taredfiles_R2 <- sapply(taredfiles,
                            function(x){grep("_R2_", x, value=TRUE) %>% head(1)})
  }
  for(i in 1:length(sampleDirs)){
    untar(sampleDirs[i], files=taredfiles_R2[i])
  }
  taredfiles_R2 <- normalizePath(taredfiles_R2)
  dataset$`Read1` <- taredfiles_R2
  input <- EzDataset(meta=dataset, dataRoot=NULL)
  
  if(sum(input$meta$`Read Count`) > 1e9){
    input <- ezMethodSubsampleFastq(input=input, param=param)
  }
  inputRaw <- input$copy()
  
  # Preprocessing
  param$trimAdapter = TRUE
  input <- ezMethodFastpTrim(input = input, param = param)
  
  # fastqscreen part
  ## get Adapter contamination from raw data
  confFile = FASTQSCREEN_ADAPTER_CONF
  files_rawData = inputRaw$getFullPaths("Read1")
  resultFiles_rawData = executeFastqscreenCMD(param, confFile = confFile,
                                              files_rawData)
  fastqData_rawData = collectFastqscreenOutput(files_rawData,
                                               resultFiles_rawData)
  tempFns <- list.files(".",
                        pattern=".*(tagged\\.fastq|tagged_filter\\.fastq)$")
  file.remove(tempFns)
  tempFns <- list.files(".", 
                        pattern=".*(screen\\.html|screen\\.png)$")
  file.remove(tempFns)
  file.remove(resultFiles_rawData)
  
  ## PreprocessedData
  confFile = FASTQSCREEN_GENOMICDNA_RIBORNA_CONF
  files_ppData = input$getFullPaths("Read1")
  resultFiles_ppData = executeFastqscreenCMD(param, confFile = confFile, 
                                             files_ppData)
  fastqData_ppData = collectFastqscreenOutput(files_ppData, resultFiles_ppData)
  noHit_files = gsub('.fastq.gz$', '.tagged_filter.fastq.gz', files_ppData)
  readCount = ezFrame(totalReadCount = integer(length(files_ppData)), 
                      unmappedReadCount = integer(length(files_ppData)),
                      row.names=names(files_ppData))
  readCount[ , "totalReadCount"] <- countReadsInFastq(files_ppData)
  readCount[ , "unmappedReadCount"] <- countReadsInFastq(noHit_files)
  file.remove(resultFiles_ppData)
  tempFns <- list.files(".", pattern=".*tagged\\.fastq$")
  file.remove(tempFns)
  tempFns <- list.files(".", 
                        pattern=".*(screen\\.html|screen\\.png)$")
  file.remove(tempFns)
  
  # bowtie2 reference part
  countFiles = executeBowtie2CMD(param, input)
  speciesPercentageTop = collectBowtie2Output(param, countFiles, 
                                              readCount, virusResult=F)
  
  #Always check human data for viruses
  if(grepl('^Human|^Homo',dataset$Species[1])){
    param[['virusCheck']] = T
  }
  
  if(param[['virusCheck']]){
    countFiles = executeBowtie2CMD_Virus(param, noHit_files)
    speciesPercentageTopVirus = collectBowtie2Output(param, countFiles, 
                                                     readCount, virusResult = T)
    dir.create('virusCheck')
    for (i in 1:length(countFiles)){
      ezSystem(paste('mv', file.path(countFiles[i]), 'virusCheck'))
    }
  } else {
    speciesPercentageTopVirus = NULL
  }
  
  file.remove(noHit_files)
  
  rRNA_strandInfo = get_rRNA_Strandness(param, input)
  krakenResult = runKraken(param, input)
  
  file.remove(input$getFullPaths("Read1"))

  # debug
 # save(fastqData_ppData, fastqData_rawData, speciesPercentageTop, krakenResult,
 #      dataset, param, rRNA_strandInfo, file="fastqscreen.rda")
  
  #create report
  setwdNew(basename(output$getColumn("Report")))
  if(param[['virusCheck']]){
      ezSystem(paste('mv', '../virusCheck', '.'))
  }
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "FastqScreen.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  
  ## generate the main reports
  rmarkdown::render(input="FastqScreen.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  
  prepareRmdLib()
  
  unlink(dirname(taredfiles_R1), recursive=TRUE)
  
  return("Success")
}

EzAppFastqScreen_10x <-
  setRefClass("EzAppFastqScreen_10x",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodFastqScreen_10x
                  name <<- "EzAppFastqScreen_10x"
                  appDefaults <<- rbind(nTopSpecies=ezFrame(Type="integer",
                                                            DefaultValue=10,
                                                            Description="number of species to show in the plots"),
                                        confFile=ezFrame(Type="character",
                                                         DefaultValue="",
                                                         Description="the configuration file for fastq screen"),
                                        virusCheck=ezFrame(Type="logical",
                                                           DefaultValue=FALSE,
                                                           Description="check for viruses in unmapped data"),
                                        minAlignmentScore=ezFrame(Type="integer",
                                                                  DefaultValue="-20",
                                                                  Description="the min alignment score for bowtie2"),
                                        trimAdapter=ezFrame(Type="logical",
                                                            DefaultValue=TRUE,
                                                            Description="whether to search for the adapters and trim them")
                                        )
                }
              )
  )
