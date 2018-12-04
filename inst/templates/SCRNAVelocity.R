## p2838 smart-Seq2
setwd("/scratch/gtan/dev/SCReports-p2284")

library(ezRun)
param = list()
param[['cores']] = '4'
param[['ram']] = '8'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['name']] = 'SCReport'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'


