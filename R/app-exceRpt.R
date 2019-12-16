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
                 initiallize = function()
                 {
                   "Initializes the application using its specific defaults."
                   runMethod <<- ezMethodExceRpt
                   name <<- "EzAppExceRpt"
                   appDefaults <<- ezFrame(
                     EXE_DIR   = '/export/local/miquel/exceRpt',
                     resultDir = '.',
                     
                     ADAPTER_SEQ = 'guessKnown',
                     SAMPLE_NAME = 'NULL',
                     MAIN_ORGANISM_GENOME_ID = 'hg38',
                     CALIBRATOR_LIBRARY = 'NULL',
                     CALIBRATOR_TRIM5p = 0,
                     CALIBRATOR_TRIM3p = 0,
                     ENDOGENOUS_LIB_PRIORITY = 'miRNA,tRNA,piRNA,gencode,circRNA',
                     
                     TRIM_N_BASES_5p = '0',
                     TRIM_N_BASES_3p = '0',
                     RANDOM_BARCODE_LENGTH = 0,
                     RANDOM_BARCODE_LOCATION = "'-5p -3p'",
                     KEEP_RANDOM_BARCODE_STATS = 'false',
                     DOWNSAMPLE_RNA_READS = 'NULL',
                     MAP_EXOGENOUS = 'off',
                     
                     MAX_RAM = '16G',
                     N_THREADS = 4,
                     JAVA_RAM = '10G',
                     REMOVE_LARGE_INTERMEDIATE_FILES = 'false',
                     
                     MIN_READ_LENGTH = 18,
                     QFILTER_MIN_QUAL = 20,
                     QFILTER_MIN_READ_FRAC = 80,
                     STAR_alignEndsType = 'Local',
                     STAR_outFilterMatchNmin = 18,
                     STAR_outFilterMatchNminOverLread = 0.9,
                     STAR_outFilterMismatchNmax = 1,
                     MAX_MISMATCHES_EXOGENOUS = 0
                   )
                 }
               )
               
  )

ezMethodExceRpt = function(input=NA, output=NA, param=NA){
  # initialize method defaults
  #readFile = input$getFullPaths("Read1") # input file path
  readFile = '/export/local/miquel/exceRpt/ExampleData/testData_human.fastq.gz'
  
  cmd = paste0('make -f ',param[['EXE_DIR']],'/exceRpt_smallRNA ', # makefile to execute
               
  ## required options
  'INPUT_FILE_PATH=',readFile,' ',    # <Path>
  'OUTPUT_DIR=',param[['resultDir']],'/output ',                    # <Path>
  
  'EXE_DIR=',param[['EXE_DIR']],' ',       # <Path>
  
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
  'MAX_RAM=',param[['MAX_RAM']],' ',                                                 # <String>
  'N_THREADS=',param[['N_THREADS']],' ',                                             # <int>
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

  ezSystem(cmd)
  
  ## Function to delete all files in OUTPUT_DIR that are not in CORE_RESULTS; remove unnecessary
  ## data to create report.
  keepOnlyCoreFiles(param[['resultDir']])
}

