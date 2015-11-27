
library(ezRun)
cd = getwd()
setwdNew("/scratch/test_fastqscreen")

param = list()
param[['cores']] = '8'
param[['ram']] = '40'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'true'
param[['name']] = 'FastqScreen_Result'
param[['nReads']] = '100000'
param[['nTopSpecies']] = '5'
param[['minAlignmentScore']] = '-20'
param[['confFile']] = 'variousSpecies_rRNA_20140901_silva119.conf'
param[['cmdOptions']] = '-k 10 --trim5 4 --trim3 4 --very-sensitive'
param[['mail']] = 'peter.schmid@ieu.uzh.ch'
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['resultDir']] = 'p1001/QC_FastqScreen_734_TEST_2015-10-16--07-52-36'
output = list()
output[['Name']] = 'FastqScreen_Result'
output[['Report [File]']] = 'p1001/QC_FastqScreen_734_TEST_2015-10-16--07-52-36/FastqScreen_Result'
output[['Html [Link]']] = 'p1001/QC_FastqScreen_734_TEST_2015-10-16--07-52-36/FastqScreen_Result/00index.html'
input = EzDataset(file=system.file("extdata/ventricles_100k/dataset.tsv", package="ezRun", mustWork = TRUE))
FASTQSCREEN="/usr/local/ngseq/opt/fastq_screen_v0.5.2/fastq_screen"
FASTQSCREEN_CONF_DIR="/usr/local/ngseq/opt/fastq_screen_v0.5.2/conf/"
BOWTIE2_DIR="/usr/local/ngseq/src/bowtie2-2.2.6"
SAMTOOLS="/usr/local/ngseq/stow/samtools-1.2/bin/samtools"

myApp = EzAppFastqScreen$new()
myApp$run(input=input$copy()$subset(1:2), output=output, param=param)

setwd(cd)
