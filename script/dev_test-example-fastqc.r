
library(ezRun)
setwdNew("/scratch/fastqc_test")
param = list()
param[['cores']] = '4'
param[['ram']] = '4'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'true'
param[['name']] = 'FastQC_Result'
param[['cmdOptions']] = ''
param[['mail']] = ''
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['resultDir']] = 'p1001/QC_Fastqc_5750_2015-09-15--09-38-58'
output = list()
output[['Name']] = 'FastQC_Result'
output[['Report [File]']] = 'p1001/QC_Fastqc_5750_2015-09-15--09-38-58/FastQC_Result'
output[['Html [Link]']] = 'p1001/QC_Fastqc_5750_2015-09-15--09-38-58/FastQC_Result/00index.html'
input = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)

myApp = EzAppFastqc$new()
myApp$run(input=input, output=output, param=param)
