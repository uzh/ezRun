## p2401 Single-end samples and with correlation
setwd("/srv/GT/analysis/gtan/p2401-FastqScreenRmarkdown")
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
param = list()
param[['cores']] = '8'
param[['ram']] = '40'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'false'
param[['name']] = 'FastqScreen_Result'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '4'
param[['trimRight']] = '4'
param[['minTailQuality']] = '20'
param[['minAvgQuality']] = '20'
param[['minReadLength']] = '40'
param[['nReads']] = '100000'
param[['nTopSpecies']] = '5'
param[['minAlignmentScore']] = '-20'
param[['confFile']] = 'variousSpecies_rRNA_20160826_silva123.conf'
param[['cmdOptions']] = '-k 10 --very-sensitive'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2401/FastqScreen_17776_2017-04-24--09-10-26'
output = list()
output[['Name']] = 'FastqScreen_Result'
output[['Report [File]']] = 'p2401/FastqScreen_17776_2017-04-24--09-10-26/FastqScreen_Result'
output[['Html [Link]']] = 'p2401/FastqScreen_17776_2017-04-24--09-10-26/FastqScreen_Result/00index.html'
input = '/srv/GT/analysis/gtan/p2401-FastqScreenRmarkdown/input_dataset.tsv'
EzAppFastqScreen$new()$run(input=input, output=output, param=param)

# p2000 single cell unmapped bam
library(ezRun)
setwd("/scratch/gtan/p2000-SCFastqScreen")
setEnvironments("flexbar")
setEnvironments("fastq_screen")
param = list()
param[['cores']] = '8'
param[['ram']] = '40'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'false'
param[['name']] = 'FastqScreen_Result'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '4'
param[['trimRight']] = '4'
param[['minTailQuality']] = '20'
param[['minAvgQuality']] = '20'
param[['minReadLength']] = '40'
param[['nReads']] = '100000'
param[['nTopSpecies']] = '5'
param[['minAlignmentScore']] = '-20'
param[['confFile']] = 'variousSpecies_rRNA_20160826_silva123.conf'
param[['cmdOptions']] = '-k 10 --very-sensitive'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'projects/p2000/SingleCellUnmappedBam'
output = list()
output[['Name']] = 'FastqScreen_Result'
output[['Report [File]']] = 'p2401/FastqScreen_17776_2017-04-24--09-10-26/FastqScreen_Result'
output[['Html [Link]']] = 'p2401/FastqScreen_17776_2017-04-24--09-10-26/FastqScreen_Result/00index.html'
#input = '/srv/gstore/projects/p2000/SingleCellUnmappedBam/dataset.tsv'
input = 'dataset.tsv'
#debug(ezMethodFastqScreen)
EzAppFastqScreen$new()$run(input=input, output=output, param=param)

