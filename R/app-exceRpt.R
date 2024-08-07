###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title exceRpt_smallRNA app
##' @description Use this reference class to run exceRpt_smallRNA parameters.
##' @author Miquel Anglada Girotto
EzAppExceRpt =
  setRefClass( "EzAppExceRpt",
               contains = "EzApp",
               methods = list(
                 initialize = function(){
                   "Initializes the application using its specific defaults."
                   runMethod <<- ezMethodExceRpt
                   name <<- "EzAppExceRpt"
                 }
               )
  )

##' @title execute exceRpt_smallRNA.
##' @description R wrapper to execute exceRpt_smallRNA pipeline to count smallRNAs in bulk RNA-seq experiments.
##' @param input ezDataFrame()
##' @param output ezDataFrame()
##' @param param list() exceRpt_smallRNA and environment parameters.
##' @author Miquel Anglada Girotto
ezMethodExceRpt = function(input=NA, output=NA, param=NA){
  ## objects from dataset information
  readFile = input$getFullPaths("Read1") # input file path
  param[['ADAPTER_SEQ']] = input$getColumn('Adapter1')
  stopifnot(length(param[['ADAPTER_SEQ']]) == 1)
  
  outputDir = basename(output$getColumn("excerpt"))
  dir.create(outputDir)
  
  ## objects from sushi form
  param[['SAMPLE_NAME']] <- input$getNames()
  param[['refBuild']] = gsub('.*/', '', param[['refBuild']])
  
  ## create command to run exceRpt_smallRNA
  cmd = pasteExcerptCmd(readFile, outputDir, param)
  if(ezIsSpecified(param[['cmdOptions']])) 
    cmd = paste(cmd, param[['cmdOptions']])
  ezSystem(cmd)
  
  ## Function to delete all files in OUTPUT_DIR that are not in CORE_RESULTS; remove unnecessary
  ## data to create report.
  if(param[['REMOVE_LARGE_INTERMEDIATE_FILES']]=='true'){
    keepOnlyCoreFiles(outputDir)
  }
}


##' @title create appropriate command line call.
##' @description paste character strings to create the adequate shell command to execute exceRpt_smallRNA. 
##' @param param list() exceRpt_smallRNA parameters.
##' @param readFile <string> full path of the *.fastq(.gz) input file.
##' @param outputDir <string> directory in which to save the processed outputs created by exceRpt_smallRNA.
##' @author Miquel Anglada Girotto
pasteExcerptCmd=function(readFile, outputDir, param){
  cmd = paste(
        ## makefile to execute
        paste0('make -f ', file.path(Sys.getenv("exceRpt"), 'exceRpt_smallRNA')),
        paste0('MAX_RAM=', param[['ram']], 'G'),
        paste0('N_THREADS=', param[['cores']]),
        
        ## required options
        paste0('INPUT_FILE_PATH=',readFile),
        paste0('OUTPUT_DIR=',outputDir),
        
        paste0('EXE_DIR=', Sys.getenv("exceRpt")),
        
        ## paths to our modules
        paste0('FASTX_CLIP_EXE=', Sys.which("fastx_clipper")),
        paste0('FASTX_FILTER_EXE=', Sys.which("fastq_quality_filter")),
        paste0('BOWTIE2_EXE=', Sys.which("bowtie2")),
        paste0('SAMTOOLS_EXE=', Sys.which("samtools")),
        paste0('JAVA_EXE=', Sys.which("java")),
        paste0('FASTQC_EXE=', Sys.which("fastqc")),
        paste0('STAR_EXE=', Sys.which("STAR")),
        
        ## main analysis options
        paste0('ADAPTER_SEQ=', param[['ADAPTER_SEQ']]),
        paste0('SAMPLE_NAME=', param[['SAMPLE_NAME']]),
        paste0('MAIN_ORGANISM_GENOME_ID=', param[['refBuild']]),
        paste0('ENDOGENOUS_LIB_PRIORITY=', param[['ENDOGENOUS_LIB_PRIORITY']]),
        
        ## additional analysis options
        paste0('TRIM_N_BASES_5p=',param[['TRIM_N_BASES_5p']]),
        paste0('TRIM_N_BASES_3p=',param[['TRIM_N_BASES_3p']]),
        paste0('MIN_READ_LENGTH=',param[['MIN_READ_LENGTH']]),
        paste0('QFILTER_MIN_QUAL=',param[['QFILTER_MIN_QUAL']]),
        paste0('QFILTER_MIN_READ_FRAC=',param[['QFILTER_MIN_READ_FRAC']]), 
        
        paste0('RANDOM_BARCODE_LENGTH=', param[['RANDOM_BARCODE_LENGTH']]),
        paste0('RANDOM_BARCODE_LOCATION=', "'", param[['RANDOM_BARCODE_LOCATION']], "'"),
        paste0('KEEP_RANDOM_BARCODE_STATS=', param[['KEEP_RANDOM_BARCODE_STATS']]),
        
        ## Alignment options
        paste0('STAR_alignEndsType=', param[['STAR_alignEndsType']]),
        paste0('STAR_outFilterMatchNmin=', param[['STAR_outFilterMatchNmin']]),
        paste0('STAR_outFilterMatchNminOverLread=', param[['STAR_outFilterMatchNminOverLread']]),
        paste0('STAR_outFilterMismatchNmax=', param[['STAR_outFilterMismatchNmax']]),
        
        ## output options
        paste0('REMOVE_LARGE_INTERMEDIATE_FILES=', param[['REMOVE_LARGE_INTERMEDIATE_FILES']])
  )
  return(cmd)
}


##' @title keep only files in *CORE_RESULTS* folder
##' @description to reduce the final output size, this function deletes all files within the "processed_output"
##' @description directory to keep only the core files required to generate the counts and final report.
##' @param data.dir <string> directory in which the function will be executed.
##' @author Miquel Anglada Girotto
keepOnlyCoreFiles = function(data.dir){
  # list files in output.dir
  dirFiles = list.files(data.dir,full.names = TRUE, recursive = FALSE)
  coreResultsFile = grep(pattern='*CORE_RESULTS*',dirFiles,value = TRUE)
  
  if(length(coreResultsFile)>0){
    # identify the folder with the CORE_RESULTS and not the .tgz
    toRemove = dirFiles != coreResultsFile
    # remove the rest of the files and folders
    unlink(dirFiles[toRemove],recursive=TRUE)
    
    # uncompress .tgz
    untar(tarfile = coreResultsFile, exdir = data.dir, tar=system("which tar", intern=TRUE))
    
    # remove .tgz
    file.remove(coreResultsFile)
  }else{
    print('No CORE_RESULTS found; probably they had already been extracted.')
  }
}

