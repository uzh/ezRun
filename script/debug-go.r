

cwd = getwd()

.libPaths("/home/pdschmid/R-libs")
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
library(ezRun)
ezSystem("rm -fr /scratch/test_go_debugging/*")
setwdNew("/scratch/test_go_debugging")
param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.PatchesSkipped/Annotation/Version-2015-06-25'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['testMethod']] = 'glm'
param[['grouping']] = 'Condition'
param[['sampleGroup']] = 'wound'
param[['runGfold']] = 'FALSE'
param[['refGroup']] = 'ctrl'
param[['normMethod']] = 'TMM'
param[['runGO']] = 'true'
param[['expressionName']] = ''
param[['specialOptions']] = ''
param[['mail']] = 'peter.schmid@ieu.uzh.ch'
param[['comparison']] = 'wound--over--ctrl'
param[['name']] = 'wound--over--ctrl'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1824/EdgeR_wound--over--control_8261_2015-11-03--18-18-32'
param[['writeScatterPlots']] = 'TRUE'
output = list()
output[['Name']] = 'wound--over--ctrl'
output[['Species']] = ''
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.PatchesSkipped/Annotation/Version-2015-06-25'
output[['Report [File]']] = 'p1824/EdgeR_wound--over--control_8261_2015-11-03--18-18-32/EdgeR--wound--over--ctrl'
output[['Html [Link]']] = 'p1824/EdgeR_wound--over--control_8261_2015-11-03--18-18-32/EdgeR--wound--over--ctrl/00index.html'
input = '/srv/gstore/projects/p1824/EdgeR_wound--over--control_8261_2015-11-03--18-18-32/input_dataset.tsv'
EzAppEdger$new()$run(input=input, output=output, param=param)

setwd(cwd)

param[['normMethod']] = 'logMean'
EzAppCountQC$new()$run(input=input, output=output, param=param)

setwd(cwd)
