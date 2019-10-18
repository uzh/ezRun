## p2886 # o6162 Sample1
setwd("/scratch/gtan/dev/p2886-Reza-SCSeurat")

library(ezRun)
param = list()
param[['cores']] = '1'
param[['ram']] = '16'
param[['scratch']] = '50'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['samples']] = 'Sample_1,Sample_2'
param[['name']] = 'SCSeurat'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = ''
param[['scProtocol']] = '10x'
param[['minGenesPerCell']] = '1000'
param[['maxGenesPerCell']] = '7000'
param[['maxMitoPercent']] = '25'
param[['vars.to.regress']] = ''
param[['pcs']] = '30'
param[['knownMarkers']] = 'mylists=SOX2,NES,PAX6;'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2284/SCSeurat_28929_2018-10-05--22-01-26'
output = list()
output[['Name']] = 'Sample_1'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2284/SCSeurat_28929_2018-10-05--22-01-26/Sample_1_17_08_2018_SCReport/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-176.uzh.ch/shiny/fgcz_exploreSCSeurat_app/?data=p2284/SCReport_28929_2018-10-05--22-01-26/AVM_17_08_2018_SCReport/SCSeurat-butvkxmjqqvu.rds'
output[['Report [File]']] = 'p2284/SCSeurat_28929_2018-10-05--22-01-26/Sample_1_17_08_2018_SCSeurat'
input = list()
input[['Name']] = 'Sample_1'
input[['ResultDir']] = 'p2886/CellRangerCount_40261_2019-10-10--15-57-41/Sample_1'
input[['Report']] = 'p2886/CellRangerCount_40261_2019-10-10--15-57-41/Sample_1/web_summary.html'
input[['Species']] = 'Homo sapiens (human)'
input[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
input[['CountMatrix']] = 'p2886/CellRangerCount_40261_2019-10-10--15-57-41/Sample_1/filtered_feature_bc_matrix'
input[['refFeatureFile']] = 'genes.gtf'
input[['featureLevel']] = 'gene'

# input = EzDataset$new(meta=input, dataRoot=param$dataRoot)
# param <- ezParam(param)
# sce <- loadSCCountDataset(input, param)
# debug(ezMethodSCSeurat)
EzAppSCSeurat$new()$run(input=input, output=output, param=param)

