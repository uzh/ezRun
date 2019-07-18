
#p2860 integration of wt samples
setwd("/scratch/gtan/dev/SCScranIntegration-p2860")
library(ezRun)
param = list()
param[['cores']] = '1'
param[['ram']] = '30'
param[['scratch']] = '50'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['samples']] = 'wt_3_M,wt_4_F'
param[['name']] = 'SCScranIntegration'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = 'protein_coding'
param[['scProtocol']] = '10X'
param[['batchCorrection']] = 'MNN'
param[['runPseudoTime']] = 'false'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['Rversion']] = 'Dev/R/3.6.0'
param[['dataRoot']] = '/srv/gstore/projects'
param <- ezParam(param)



