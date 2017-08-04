setwd("/home/gtan/analysis/gtan/p2438-RNABamStats")
Sys.setenv("PATH"=paste("/usr/local/ngseq/packages/Dev/Python2/2.7.13/bin", "/usr/local/ngseq/packages/Tools/samtools/1.5/bin", Sys.getenv("PATH"), sep=":"))
library(ezRun)
param = list()
param[['cores']] = '16'
param[['ram']] = '50'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'false'
param[['name']] = 'RNA_BAM_Statistics'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
param[['refFeatureFile']] = 'genes.gtf'
param[['strandMode']] = 'antisense'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2438/RNABamStats_18683_2017-06-12--16-29-29'
param[['splitByChrom']] = 'false'
output = list()
output[['Name']] = 'RNA_BAM_Statistics'
output[['Report [File]']] = 'p2438/RNABamStats_18683_2017-06-12--16-29-29/RNA_BAM_Statistics'
output[['Html [Link]']] = 'p2438/RNABamStats_18683_2017-06-12--16-29-29/RNA_BAM_Statistics/00index.html'
output[['Species']] = ''
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
output[['refFeatureFile']] = 'genes.gtf'
input = 'input_dataset.tsv'
debug("getTargetTypeCounts")
EzAppRnaBamStats$new()$run(input=input, output=output, param=param)

# debug
input = EzDataset$new(file=input, dataRoot=param$dataRoot)
output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
param = ezParam(param)
