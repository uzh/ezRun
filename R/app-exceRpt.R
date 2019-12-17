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
                     MAIN_ORGANISM_GENOME_ID = ezFrame(Type="character", DefaultValue='hg38', Description="'hg38'|'hg19'|'mm10'"),
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
                     REMOVE_LARGE_INTERMEDIATE_FILES = ezFrame(Type="character", DefaultValue='false', Description="'false'|'true'"),
                     
                     MIN_READ_LENGTH = ezFrame(Type="character", DefaultValue='18', Description="Filter out reads not achieving minimum length."),
                     QFILTER_MIN_QUAL = ezFrame(Type="character", DefaultValue='20', Description="Filter out reads not achieving minimum read quality."),
                     QFILTER_MIN_READ_FRAC = ezFrame(Type="character", DefaultValue='80', Description=""),
                     STAR_alignEndsType = ezFrame(Type="character", DefaultValue='Local', Description="'Local'|'EndToEnd'"),
                     STAR_outFilterMatchNmin = ezFrame(Type="character", DefaultValue='18', Description=""),
                     STAR_outFilterMatchNminOverLread = ezFrame(Type="character", DefaultValue='0.9', Description=""),
                     STAR_outFilterMismatchNmax = ezFrame(Type="character", DefaultValue='1', Description=""),
                     MAX_MISMATCHES_EXOGENOUS = ezFrame(Type="character", DefaultValue='0', Description=""))
                 }
               )
               
  )

ezMethodExceRpt = function(input=NA, output=NA, param=NA){
  
  # init variables
  readFile = input$getFullPaths("Read1") # input file path

  outputDir = output$getColumn("excerpt")
  if(!dir.exists(outputDir)){dir.create(outputDir)}
  
  ## create sample name if not specified
  if(param[['SAMPLE_NAME']]=='NULL'){
    param[['SAMPLE_NAME']] = gsub('.gz','',gsub('.fastq','',gsub('.*/','',readFile)))
  }
  
  # create command to run exceRpt_smallRNA
  cmd = pasteCmd(param=param)

  ezSystem(cmd)
  
  ## Function to delete all files in OUTPUT_DIR that are not in CORE_RESULTS; remove unnecessary
  ## data to create report.
  keepOnlyCoreFiles(outputDir)
}

##' paste shell command
##' 

pasteCmd(param){
  cmd = paste0('make -f ',param[['EXE_DIR']],'/exceRpt_smallRNA ', # makefile to execute
         
         ## required options
         'INPUT_FILE_PATH=',readFile,' ',    # <Path>
         'OUTPUT_DIR=',outputDir,' ',                    # <Path>
         
         'EXE_DIR=',param[['EXE_DIR']],' ',       # <Path>
         
         ## paths to our modules
         'FASTX_CLIP_EXE=',param[['FASTX_CLIP_EXE']],' ',
         'FASTX_FILTER_EXE=',param[['FASTX_FILTER_EXE']],' ',
         'BOWTIE2_EXE=',param[['BOWTIE2_EXE']],' ',
         'SAMTOOLS_EXE=',param[['SAMTOOLS_EXE']],' ',
         'JAVA_EXE=',param[['JAVA_EXE']],' ',
         'FASTQC_EXE=',param[['FASTQC_EXE']],' ',
         'STAR_EXE=',param[['STAR_EXE']],' ',
         
         ## main analysis options
         'ADAPTER_SEQ=',param[['ADAPTER_SEQ']],' ',                         # 'gessKnown'|'none'|<String>
         'SAMPLE_NAME=',param[['SAMPLE_NAME']],' ',                         # <String>
         'MAIN_ORGANISM_GENOME_ID=',param[['MAIN_ORGANISM_GENOME_ID']],' ', # 'hg38'|'hg19'|'mm10'
         'CALIBRATOR_LIBRARY=',param[['CALIBRATOR_LIBRARY']],' ',           # <Path>, path to a bowtie2 index of calibrator oligos used for QC or normalisation
         'CALIBRATOR_TRIM5p=',param[['CALIBRATOR_TRIM5p']],' ',             # <int>
         'CALIBRATOR_TRIM3p=',param[['CALIBRATOR_TRIM3p']],' ',             # <int>
         'ENDOGENOUS_LIB_PRIORITY=',param[['ENDOGENOUS_LIB_PRIORITY']],' ', # <comma,separated,list,no,spaces>
         
         ## additional analysis options
         'TRIM_N_BASES_5p=',param[['TRIM_N_BASES_5p']],' ',                     # <int>
         'TRIM_N_BASES_3p=',param[['TRIM_N_BASES_3p']],' ',                     # <int>
         'RANDOM_BARCODE_LENGTH=',param[['RANDOM_BARCODE_LENGTH']],' ',         # <int>
         'RANDOM_BARCODE_LOCATION=',param[['RANDOM_BARCODE_LOCATION']],' ',     # '-5p -3p'/'-5p'/'-3p'
         'KEEP_RANDOM_BARCODE_STATS=',param[['KEEP_RANDOM_BARCODE_STATS']],' ', # 'false'|'true'
         'DOWNSAMPLE_RNA_READS=',param[['DOWNSAMPLE_RNA_READS']],' ',           # <int>
         'MAP_EXOGENOUS=',param[['MAP_EXOGENOUS']],' ',                         # 'off'|'miRNA'|'on'
         
         ## Hardware-specific options
         'MAX_RAM=',paste0(param[['ram']],'G'),' ',                                                 # <String>
         'N_THREADS=',param[['cores']],' ',                                             # <int>
         'JAVA_RAM=',param[['JAVA_RAM']],' ',                                               # <String>
         'REMOVE_LARGE_INTERMEDIATE_FILES=',param[['REMOVE_LARGE_INTERMEDIATE_FILES']],' ', # 'false'|'true'
         
         ## Alignment QC options
         'MIN_READ_LENGTH=',param[['MIN_READ_LENGTH']],' ',                                   # <int>
         'QFILTER_MIN_QUAL=',param[['QFILTER_MIN_QUAL']],' ',                                 # <int>
         'QFILTER_MIN_READ_FRAC=',param[['QFILTER_MIN_READ_FRAC']],' ',                       # <double>
         'STAR_alignEndsType=',param[['STAR_alignEndsType']],' ',                             # 'Local'|'EndToEnd'
         'STAR_outFilterMatchNmin=',param[['STAR_outFilterMatchNmin']],' ',                   # <int>
         'STAR_outFilterMatchNminOverLread=',param[['STAR_outFilterMatchNminOverLread']],' ', # <double>
         'STAR_outFilterMismatchNmax=',param[['STAR_outFilterMismatchNmax']],' ',             # <int>
         'MAX_MISMATCHES_EXOGENOUS=',param[['MAX_MISMATCHES_EXOGENOUS']]                      # <int>
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
