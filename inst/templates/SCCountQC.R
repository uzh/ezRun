setwd("/export/local/scratch/gtan/p2497-SCCountQC")

library(ezRun)
param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['name']] = 'SCCount_QC'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2497/p2497-SCCountQC'
output = list()
output[['Name']] = 'SCCount_QC'
output[['Species']] = ''
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
output[['Static Report [Link]']] = 'p2497/p2497-SCCountQC/SCCount_QC/00index.html'
output[['Report [File]']] = 'p2497/p2497-SCCountQC/SCCount_QC'
input = 'input_dataset.tsv'

input = EzDataset$new(file=input, dataRoot=param$dataRoot)
param <- ezParam(param)
output = EzDataset$new(meta=output, dataRoot=param$dataRoot)

sce <- loadSCCountDataset(input, param)
EzAppSCCountQC$new()$run(input=input, output=output, param=param)


