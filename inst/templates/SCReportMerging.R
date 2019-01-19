# p2529 merge 
setwd("/scratch/gtan/dev/SCReportMerging-p2529")
library(ezRun)
param = list()
param[['cores']] = '4'
param[['ram']] = '30'
param[['scratch']] = '50'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['samples']] = c('FAPS4_all', 'FAPS30_all_E1')
param[['name']] = 'SCReport'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = ''
param[['resolution']] = '0.3'
param[['batchCorrection']] = 'CCA'
param[['cc']] = '20'
param[['chosenClusters1']] =''
param[['chosenClusters2']] ='0,1,2,3,4,5,6,8,9'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2529/SCReport_32851_FAPS4_all_2019-01-08--10-47-27'
output = list()
output[['Name']] = 'FAPS4_all'
output[['Species']] = 'Mus musculus (house mouse)'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2529/SCReport_32851_FAPS4_all_2019-01-08--10-47-27/FAPS4_all_SCReport/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreSingleCell_app/?data=p2529/SCReport_32851_FAPS4_all_2019-01-08--10-47-27/FAPS4_all_SCReport/SCReport-emituitwiiek.rds'
output[['Report [File]']] = 'p2529/SCReport_32851_FAPS4_all_2019-01-08--10-47-27/FAPS4_all_SCReportMerging'
output[['ResultDir [Link]']] = 'p2529/CellRangerCount_25012_2019-01-07--16-20-19/FAPS4_all'
input <- "/scratch/gtan/dev/SCReportMerging-p2529/SCReport_FAPS4_FAPS30.tsv"

# input = EzDataset$new(file=input, dataRoot=param$dataRoot)
# output <- EzDataset$new(meta=output, dataRoot=param$dataRoot)
# param <- ezParam(param)
#debug(ezMethodSCReportMerging)
EzAppSCReportMerging$new()$run(input=input, output=output, param=param)
