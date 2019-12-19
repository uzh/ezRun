###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


EzAppExceRpt =
  setRefClass( "EzAppExceRpt",
               contains = "EzApp",
               methods = list(
                 initialize = function(){
                   "Initializes the application using its specific defaults."
                   runMethod <<- ezMethodExceRpt
                   name <<- "EzAppExceRpt"
                   appDefaults <<- rbind(
                     EXE_DIR = ezFrame(Type="character", DefaultValue='/usr/local/ngseq/packages/Tools/exceRpt/4.6.3', Description="Path to smallRNA makefile"),
                     
                     FASTX_CLIP_EXE=ezFrame(Type="character", DefaultValue='/usr/local/ngseq/bin/fastx_clipper', Description="Path to binary"),
                     FASTX_FILTER_EXE=ezFrame(Type="character", DefaultValue='/usr/local/ngseq/bin/fastq_quality_filter', Description="Path to binary"),
                     BOWTIE2_EXE=ezFrame(Type="character", DefaultValue='/usr/local/ngseq/packages/Aligner/Bowtie2/2.3.2/bin/bowtie2', Description="Path to binary"),
                     SAMTOOLS_EXE=ezFrame(Type="character", DefaultValue='/usr/local/ngseq/packages/Tools/samtools/1.8/bin/samtools', Description="Path to binary"),
                     JAVA_EXE=ezFrame(Type="character", DefaultValue='/usr/local/ngseq/packages/Dev/jdk/8/bin/java', Description="Path to binary"),
                     FASTQC_EXE=ezFrame(Type="character", DefaultValue="'/usr/local/ngseq/packages/Dev/jdk/8/bin/java -classpath /usr/local/ngseq/packages/QC/FastQC/0.11.7:/usr/local/ngseq/packages/QC/FastQC/0.11.7/sam-1.103.jar:/usr/local/ngseq/packages/QC/FastQC/0.11.7/jbzip2-0.9.jar'", Description="Paths to java and fastqc binaries"),
                     STAR_EXE=ezFrame(Type="character", DefaultValue='/usr/local/ngseq/packages/Aligner/STAR/2.5.2b/bin/STAR', Description="Path to binary"),
                     
                     ADAPTER_SEQ = ezFrame(Type="character", DefaultValue='guessKnown', Description="'gessKnown'|'none'|<String>"),
                     SAMPLE_NAME = ezFrame(Type="character", DefaultValue='NULL', Description="<String>"),
                     refBuild = ezFrame(Type="character", DefaultValue='Homo_sapiens/UCSC/hg38', Description="'Mus_musculus/UCSC/mm10', 'Homo_sapiens/UCSC/hg38'"),
                     CALIBRATOR_LIBRARY = ezFrame(Type="character", DefaultValue='NULL', Description="<Path>, path to a directory containing bowtie2 index of calibrator oligos used for QC or normalisation."),
                     CALIBRATOR_TRIM5p = ezFrame(Type="character", DefaultValue='0', Description="<int>"),
                     CALIBRATOR_TRIM3p = ezFrame(Type="character", DefaultValue='0', Description="<int>"),
                     ENDOGENOUS_LIB_PRIORITY = ezFrame(Type="character", DefaultValue='miRNA,tRNA,piRNA,gencode,circRNA', Description="<comma,separated,list,no,spaces>"),
                     
                     TRIM_N_BASES_5p = ezFrame(Type="character", DefaultValue='0', Description="<int>"),
                     TRIM_N_BASES_3p = ezFrame(Type="character", DefaultValue='0', Description="<int>"),
                     RANDOM_BARCODE_LENGTH = ezFrame(Type="character", DefaultValue='0', Description="<int>"),
                     RANDOM_BARCODE_LOCATION = ezFrame(Type="character", DefaultValue="'-5p -3p'", Description="'-5p -3p'/'-5p'/'-3p'"),
                     KEEP_RANDOM_BARCODE_STATS = ezFrame(Type="character", DefaultValue='false', Description="'false'|'true'"),
                     DOWNSAMPLE_RNA_READS = ezFrame(Type="character", DefaultValue='NULL', Description="<int>"),
                     MAP_EXOGENOUS = ezFrame(Type="character", DefaultValue='off', Description="'off'|'miRNA'|'on'"),
                     
                     ram = ezFrame(Type="character", DefaultValue='16G', Description="Maximum amount of RAM."),
                     cores = ezFrame(Type="character", DefaultValue='4', Description="Number of threads."),
                     JAVA_RAM = ezFrame(Type="character", DefaultValue='10G', Description="RAM amount allocated for Java."),
                     REMOVE_LARGE_INTERMEDIATE_FILES = ezFrame(Type="character", DefaultValue='true', Description="'false'|'true'"),
                     
                     MIN_READ_LENGTH = ezFrame(Type="character", DefaultValue='18', Description="Filter out reads not achieving minimum length."),
                     QFILTER_MIN_QUAL = ezFrame(Type="character", DefaultValue='20', Description="Filter out reads not achieving minimum read quality."),
                     QFILTER_MIN_READ_FRAC = ezFrame(Type="character", DefaultValue='80', Description=""),
                     STAR_alignEndsType = ezFrame(Type="character", DefaultValue='Local', Description="'Local'|'EndToEnd'"),
                     STAR_outFilterMatchNmin = ezFrame(Type="character", DefaultValue='18', Description=""),
                     STAR_outFilterMatchNminOverLread = ezFrame(Type="character", DefaultValue='0.9', Description=""),
                     STAR_outFilterMismatchNmax = ezFrame(Type="character", DefaultValue='1', Description=""),
                     MAX_MISMATCHES_EXOGENOUS = ezFrame(Type="character", DefaultValue='0', Description="")
                     )
                 }
               )
               
  )

ezMethodExceRpt = function(input=NA, output=NA, param=NA){
  ## objects from dataset information
  readFile = input$getFullPaths("Read1") # input file path
  param[['ADAPTER_SEQ']] = input$getColumn('Adapter1')
  stopifnot(length(param[['ADAPTER_SEQ']]) == 1)
  
  outputDir = output$getColumn("excerpt")
  if(!dir.exists(outputDir)){dir.create(outputDir)}
  
  ## objects from sushi form
  if(param[['SAMPLE_NAME']]=='NULL'){
    param[['SAMPLE_NAME']] = gsub('.gz','',gsub('.fastq','',gsub('.*/','',readFile)))
  }
  param[['refBuild']] = gsub('.*/','',param[['refBuild']])
  
  ## create command to run exceRpt_smallRNA
  cmd = pasteCmd(param=param)

  ezSystem(cmd)
  
  ## Function to delete all files in OUTPUT_DIR that are not in CORE_RESULTS; remove unnecessary
  ## data to create report.
  if(param[['REMOVE_LARGE_INTERMEDIATE_FILES']]=='true'){keepOnlyCoreFiles(outputDir)}
}

##' paste shell command
##' 

pasteCmd=function(param,readFile,outputDir){
  cmd = paste(
        ## makefile to execute
        paste0('make -f ',param[['EXE_DIR']],'/exceRpt_smallRNA'),
         
        ## required options
        paste0('INPUT_FILE_PATH=',readFile),   
        paste0('OUTPUT_DIR=',outputDir),                   
        
        paste0('EXE_DIR=',param[['EXE_DIR']]),
        
        ## paths to our modules
        paste0('FASTX_CLIP_EXE=',param[['FASTX_CLIP_EXE']]),
        paste0('FASTX_FILTER_EXE=',param[['FASTX_FILTER_EXE']]),
        paste0('BOWTIE2_EXE=',param[['BOWTIE2_EXE']]),
        paste0('SAMTOOLS_EXE=',param[['SAMTOOLS_EXE']]),
        paste0('JAVA_EXE=',param[['JAVA_EXE']]),
        paste0('FASTQC_EXE=',param[['FASTQC_EXE']]),
        paste0('STAR_EXE=',param[['STAR_EXE']]),
        
        ## main analysis options
        paste0('ADAPTER_SEQ=',param[['ADAPTER_SEQ']]),                         
        paste0('SAMPLE_NAME=',param[['SAMPLE_NAME']]),                        
        paste0('MAIN_ORGANISM_GENOME_ID=',param[['refBuild']]),
        paste0('CALIBRATOR_LIBRARY=',param[['CALIBRATOR_LIBRARY']]),           
        paste0('CALIBRATOR_TRIM5p=',param[['CALIBRATOR_TRIM5p']]),             
        paste0('CALIBRATOR_TRIM3p=',param[['CALIBRATOR_TRIM3p']]),           
        paste0('ENDOGENOUS_LIB_PRIORITY=',param[['ENDOGENOUS_LIB_PRIORITY']]),
        
        ## additional analysis options
        paste0('TRIM_N_BASES_5p=',param[['TRIM_N_BASES_5p']]),                    
        paste0('TRIM_N_BASES_3p=',param[['TRIM_N_BASES_3p']]),                     
        paste0('RANDOM_BARCODE_LENGTH=',param[['RANDOM_BARCODE_LENGTH']]),         
        paste0('RANDOM_BARCODE_LOCATION=',param[['RANDOM_BARCODE_LOCATION']]),     
        paste0('KEEP_RANDOM_BARCODE_STATS=',param[['KEEP_RANDOM_BARCODE_STATS']]),
        paste0('DOWNSAMPLE_RNA_READS=',param[['DOWNSAMPLE_RNA_READS']]),         
        paste0('MAP_EXOGENOUS=',param[['MAP_EXOGENOUS']]),                        
        
        ## Hardware-specific options
        paste0('MAX_RAM=',paste0(param[['ram']],'G')),                                                 
        paste0('N_THREADS=',param[['cores']]),                                             
        paste0('JAVA_RAM=',param[['JAVA_RAM']]),                                               
        paste0('REMOVE_LARGE_INTERMEDIATE_FILES=',param[['REMOVE_LARGE_INTERMEDIATE_FILES']]),
        
        ## Alignment QC options
        paste0('MIN_READ_LENGTH=',param[['MIN_READ_LENGTH']]),                                   
        paste0('QFILTER_MIN_QUAL=',param[['QFILTER_MIN_QUAL']]),                                 
        paste0('QFILTER_MIN_READ_FRAC=',param[['QFILTER_MIN_READ_FRAC']]),                      
        paste0('STAR_alignEndsType=',param[['STAR_alignEndsType']]),                            
        paste0('STAR_outFilterMatchNmin=',param[['STAR_outFilterMatchNmin']]),                   
        paste0('STAR_outFilterMatchNminOverLread=',param[['STAR_outFilterMatchNminOverLread']]), 
        paste0('STAR_outFilterMismatchNmax=',param[['STAR_outFilterMismatchNmax']]),            
        paste0('MAX_MISMATCHES_EXOGENOUS=',param[['MAX_MISMATCHES_EXOGENOUS']])                     
  )
  return(cmd)
}

##' keep only core files
##' 
##' 

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
    untar(tarfile = coreResultsFile, exdir = data.dir)
    
    # remove .tgz
    file.remove(coreResultsFile)
  }else{
    print('No CORE_RESULTS found; probably they had already been extracted.')
  }
}
