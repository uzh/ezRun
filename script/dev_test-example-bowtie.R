
options(error=recover)
library(ezRun)
setwdNew("/scratch/bowtie_test")
param = list()
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
param[['paired']] = 'true'
param[['cmdOptions']] = ''
param[['trimAdapter']] = 'false'
param[['trimLeft']] = '0'
param[['trimRight']] = '0'
param[['minTailQuality']] = '0'
param[['specialOptions']] = ''
param[['mail']] = 'peter.schmid@ieu.uzh.ch'
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['resultDir']] = 'p1001/Map_Bowtie_5750_2015-09-28--17-37-04'
param[['writeIgvSessionLink']] = 'false'
# or add the correct parameters in output
output = list()
output[['Name']] = 'wt_1'
output[['BAM [File]']] = 'p1001/Map_Bowtie_5750_2015-09-28--17-37-04/wt_1.bam'
output[['BAI [File]']] = 'p1001/Map_Bowtie_5750_2015-09-28--17-37-04/wt_1.bam.bai'
output[['Species']] = 'S. cerevisiae'
output[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
output[['paired']] = 'true'
output[['Read Count']] = '9794'
output[['Genotype [Factor]']] = 'wt'
input = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
input = EzDataset$new(file=input)

myApp = EzAppBowtie$new()
myApp$run(input=input$copy()$subset(1), output=output, param=param)

