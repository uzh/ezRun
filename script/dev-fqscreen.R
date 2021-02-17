EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '40'
param[['scratch']] = '100'
param[['node']] = ''
param[['queue']] = 'employee'
param[['process_mode']] = 'DATASET'
param[['samples']] = 'o22914_1_11,o22914_1_12'
param[['paired']] = 'false'
param[['name']] = 'FastqScreen_Result'
param[['nReads']] = '100000'
param[['nTopSpecies']] = '5'
param[['minAlignmentScore']] = '-20'
param[['cmdOptions']] = '-k 10 --very-sensitive'
param[['trim_front1']] = '4'
param[['trim_tail1']] = '4'
param[['cut_front']] = 'false'
param[['cut_front_window_size']] = '4'
param[['cut_front_mean_quality']] = '20'
param[['cut_tail']] = 'false'
param[['cut_tail_window_size']] = '4'
param[['cut_tail_mean_quality']] = '20'
param[['cut_right']] = 'false'
param[['cut_right_window_size']] = '4'
param[['cut_right_mean_quality']] = '20'
param[['average_qual']] = '20'
param[['max_len1']] = '0'
param[['max_len2']] = '0'
param[['poly_x_min_len']] = '10'
param[['length_required']] = '18'
param[['cmdOptionsFastp']] = ''
param[['markDuplicates']] = 'true'
param[['mail']] = 'Hubert.Rehrauer@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p22849/FastqScreen_51734_2020-12-11--07-31-44'
param[['isLastJob']] = TRUE
output = list()
output[['Name']] = 'FastqScreen_Result'
output[['Report [File]']] = 'p22849/FastqScreen_51734_2020-12-11--07-31-44/FastqScreen_Result'
output[['Html [Link]']] = 'p22849/FastqScreen_51734_2020-12-11--07-31-44/FastqScreen_Result/00index.html'
input = '/srv/gstore/projects/p22849/FastqScreen_51734_2020-12-11--07-31-44/input_dataset.tsv'



appDefaults <<- rbind(
  nTopSpecies = ezFrame(Type = "integer", DefaultValue = 10,
                        Description = "number of species to show in the plots"),
  confFile = ezFrame(Type = "character", DefaultValue = "",
                     Description = "the configuration file for fastq screen"),
  virusCheck = ezFrame(Type = "logical", DefaultValue = FALSE,
                       Description = "check for viruses in unmapped data"),
  minAlignmentScore = ezFrame(Type = "integer", DefaultValue = "-20",
                              Description = "the min alignment score for bowtie2"),
  trimAdapter = ezFrame(Type = "logical", DefaultValue = TRUE,
                        Description = "whether to search for the adapters and trim them"),
  copyReadsLocally = ezFrame(Type = "logical", DefaultValue = FALSE,
                             Description = "copy reads to scratch first")
)

param = ezParam(param, appDefaults=appDefaults)
input = EzDataset$new(file=input, dataRoot=param$dataRoot)
input = input$subset(param$samples)
output = EzDataset$new(meta=output, dataRoot=param$dataRoot)

param$paired = FALSE

# Override the virus check parameter for human data
if (grepl("^Human|^Homo", input$getColumn("Species")[1])) {
  param[["virusCheck"]] <- T
}

if (input$readType() == "bam") {
  stopifnot(input$getLength() == 1L) ## We only support one uBam now.
  input <- ezMethodBam2Fastq(
    input = input, param = param,
    OUTPUT_PER_RG = TRUE
  )
}


inputRaw <- ezMethodSubsampleFastq(input = input, param = param, n=param$nReads)
param$trimAdapter <- TRUE
inputProc <- ezMethodFastpTrim(input = inputRaw, param = param)

## map to adapters ----
rawScreenResult = getFastqScreenStats(param,
                                      confFile = FASTQSCREEN_ADAPTER_CONF,
                                      inputRaw$getFullPaths("Read1"), workDir="rawReads")
procScreenResult = getFastqScreenStats(param,
                                       confFile = FASTQSCREEN_GENOMICDNA_RIBORNA_CONF,
                                       inputProc$getFullPaths("Read1"), workDir="procReads")
if (param[["virusCheck"]]) {
  unmappedFiles = gsub(".fastq.gz$", ".tagged_filter.fastq.gz", 
                       file.path("procReads", basename(inputProc$getFullPaths("Read1"))))
  names(unmappedFiles) =inputProc$getNames()
  virusResult <- map_and_count_virus(param, unmappedFiles, workDir="virusResult")
}  
refseqResult = map_and_count_refseq(param, inputProc$getFullPaths("Read1"), workDir="refseqResult", 
                                    readCount = inputProc$getColumn("Read Count"))

rRNAstrandResult <- get_rRNA_Strandness(param, inputProc)
krakenResult <- runKraken(param, inputProc)

file.remove(inputProc$getFullPaths("Read1"))

makeRmdReport(
  output = output, param = param,
  input=input,
  rawScreenResult=rawScreenResult, procScreenResult=procScreenResult, virusResult=virusResult,
  rRNAstrandResult=rRNAstrandResult, krakenResult=krakenResult,
  rmdFile = "FastqScreen.Rmd", reportTitle = paste("Fastq Screen", param$name)
)




