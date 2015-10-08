.libPaths("/srv/GT/analysis/course_sushi/lib")

library(ezRun)
setwdNew("/srv/GT/analysis/hubert/tmp/")
param = list()
param[['cores']] = '8'
param[['ram']] = '20'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['name']] = 'FastQC_Result'
param[['specialOptions']] = ''
param[['mail']] = 'Hubert.Rehrauer@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1001/DUMMY'
output = list()
output[['Name']] = 'FastQC_Result'
output[['Report [File]']] = 'p1001/QC_Fastqc_Dataset/FastQC_Result'
output[['Html [Link]']] = 'p1631/QC_Fastqc_Dataset/FastQC_Result/00index.html'
input = system.file("extdata/ventricles_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
#options(error=recover)
param[['dataRoot']] = system.file(package="ezRun", mustWork=TRUE)
ezRunApp(ezAppFastQC, input=input, output=output, param=param)
