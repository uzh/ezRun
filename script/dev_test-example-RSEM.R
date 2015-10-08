
options(error=recover)
library(ezRun)
setwdNew("/scratch/rsem_test")
param = list()
param[['cores']] = '8'
param[['ram']] = '30'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
param[['paired']] = 'true'
param[['strandMode']] = 'sense'
param[['refFeatureFile']] = 'genes.gtf'
param[['bowtie-e']] = '200'
param[['cmdOptions']] = ' --calc-ci '
param[['keepBam']] = 'false'
param[['trimAdapter']] = 'false'
param[['trimLeft']] = '0'
param[['trimRight']] = '0'
param[['minTailQuality']] = '0'
param[['minAvgQuality']] = '20'
param[['specialOptions']] = ''
param[['trinityFasta']] = ''
param[['mail']] = 'peter.schmid@ieu.uzh.ch'
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['resultDir']] = 'p1001/Count_RSEM_5750_2015-09-28--16-29-35'
output = list()
output[['Name']] = 'wt_1'
output[['Count [File]']] = 'p1001/Count_RSEM_5750_2015-09-28--16-29-35/wt_1.txt'
output[['Species']] = 'S. cerevisiae'
output[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
output[['featureLevel']] = 'isoform'
output[['refFeatureFile']] = 'genes.gtf'
output[['strandMode']] = 'sense'
output[['paired']] = 'true'
output[['Read Count']] = '9794'
output[['Genotype [Factor]']] = 'wt'
input = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
input = EzDataset$new(file=input)

myApp = EzAppRSEM$new()
myApp$run(input=input$copy()$subset(1), output=output, param=param)
