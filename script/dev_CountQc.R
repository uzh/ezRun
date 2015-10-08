
library(ezRun)
setwdNew("/srv/GT/analysis/hubert/tmp/")
param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['name']] = 'Count_QC'
param[['build']] = 'Mus_musculus/Ensembl/GRCm38/Annotation/Version-2014-02-25'
param[['featureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['normMethod']] = 'logMean'
param[['expressionName']] = 'transcriptCount'
param[['runGO']] = 'false'
param[['specialOptions']] = ''
param[['mail']] = 'Hubert.Rehrauer@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1001/QC_CountQC_5450_2015-04-10--16-11-36'
output = list()
output[['Name']] = 'Count_QC'
output[['Species']] = ''
output[['build']] = 'Mus_musculus/Ensembl/GRCm38/Annotation/Version-2014-02-25'
output[['Report [File]']] = 'p1001/QC_CountQC_5450_2015-04-10--16-11-36/Count_QC'
output[['Html [Link]']] = 'p1001/QC_CountQC_5450_2015-04-10--16-11-36/Count_QC/00index.html'
input = '/srv/gstore/projects/p1001/QC_CountQC_5450_2015-04-10--16-11-36/input_dataset.tsv'
#options(error=recover)
ezRunApp(ezAppCountQC, input=input, output=output, param=param)
