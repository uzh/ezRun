setwd("/srv/GT/analysis/gtan/debug/p1536-CountQC")
library(ezRun)
param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['name']] = 'Count_QC'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['normMethod']] = 'logMean'
param[['runGO']] = 'true'
param[['backgroundExpression']] = '10'
param[['transcriptTypes']] = ''
param[['specialOptions']] = ''
param[['expressionName']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1536/CountQC_20650_2017-09-01--10-05-09'
output = list()
output[['Name']] = 'Count_QC'
output[['Species']] = ''
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
output[['Static Report [Link]']] = 'p1536/CountQC_20650_2017-09-01--10-05-09/Count_QC/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreCountQC_app/?data=p1536/CountQC_20650_2017-09-01--10-05-09/Count_QC/counts-zlpjvingquyl-EzResult.RData'
output[['Report [File]']] = 'p1536/CountQC_20650_2017-09-01--10-05-09/Count_QC'
input = '/srv/gstore/projects/p1536/CountQC_20650_2017-09-01--10-05-09/input_dataset.tsv'

# debug
#debug(ezMethodCountQC)
EzAppCountQC$new()$run(input=input, output=output, param=param)

# p2179
setwd("/scratch/gtan/p2179-CountQC")
library(ezRun)
param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['name']] = 'Count_QC'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['normMethod']] = 'logMean'
param[['runGO']] = 'true'
param[['backgroundExpression']] = '10'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA,long_noncoding,short_noncoding,pseudogene'
param[['specialOptions']] = ''
param[['expressionName']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2179/CountQC_25738_2018-03-27--22-15-25'
output = list()
output[['Name']] = 'Count_QC'
output[['Species']] = ''
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
output[['Static Report [Link]']] = 'p2179/CountQC_25738_2018-03-27--22-15-25/Count_QC/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreCountQC_app/?data=p2179/CountQC_25738_2018-03-27--22-15-25/Count_QC/counts-pmdzoiqsglrv-EzResult.RData'
output[['Report [File]']] = 'p2179/CountQC_25738_2018-03-27--22-15-25/Count_QC'
input = '/srv/gstore/projects/p2179/CountQC_25738_2018-03-27--22-15-25/input_dataset.tsv'
debug(ezMethodCountQC)
EzAppCountQC$new()$run(input=input, output=output, param=param)
