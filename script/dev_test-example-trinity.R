
library(ezRun)
setwdNew("/scratch/trinity_test")
param = list()
param[['cores']] = '12'
param[['ram']] = '20'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['name']] = 'Trinity_Assembly'
param[['paired']] = 'true'
param[['strandMode']] = 'sense'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '0'
param[['trimRight']] = '0'
param[['minTailQuality']] = '0'
param[['minAvgQuality']] = '0'
param[['minReadLength']] = '10'
param[['trinityOpt']] = '--min_kmer_cov 2'
param[['specialOptions']] = ''
param[['mail']] = 'peter.schmid@ieu.uzh.ch'
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['resultDir']] = 'p1001/Assemble_Trinity_5750_test_2015-09-08--12-35-47'
output = list()
output[['Name']] = 'Trinity_Assembly'
output[['Fasta [File]']] = 'p1001/Assemble_Trinity_5750_test_2015-09-08--12-35-47/Trinity_Assembly.fasta'
input = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)

myApp = EzAppTrinity$new()
myApp$run(input=input, output=output, param=param)
