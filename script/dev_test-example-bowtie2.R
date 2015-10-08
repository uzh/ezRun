
options(error=recover)
library(ezRun)
setwdNew("/scratch/bowtie2_test")
param = list()
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.PatchesSkipped/Annotation/Version-2015-07-05'
param[['paired']] = 'false'
param[['cmdOptions']] = '--no-unal'
param[['trimAdapter']] = 'false'
param[['trimLeft']] = '0'
param[['trimRight']] = '0'
param[['minTailQuality']] = '0'
param[['specialOptions']] = ''
param[['mail']] = 'peter.schmid@ieu.uzh.ch'
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['resultDir']] = 'p1001/Map_Bowtie2_5750_2015-07-05--21-52-41'
output = list()
output[['Name']] = 'wt_1'
output[['BAM [File]']] = 'p1001/Map_Bowtie2_5750_2015-07-05--21-52-41/wt_1.bam'
output[['BAI [File]']] = 'p1001/Map_Bowtie2_5750_2015-07-05--21-52-41/wt_1.bam.bai'
output[['IGV Starter [Link]']] = 'p1001/Map_Bowtie2_5750_2015-07-05--21-52-41/wt_1-igv.jnlp'
output[['Species']] = 'S. cerevisiae'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.PatchesSkipped/Annotation/Version-2015-07-05'
output[['paired']] = 'false'
output[['Read Count']] = '9794'
output[['IGV Starter [File]']] = 'p1001/Map_Bowtie2_5750_2015-07-05--21-52-41/wt_1-igv.jnlp'
output[['IGV Session [File]']] = 'p1001/Map_Bowtie2_5750_2015-07-05--21-52-41/wt_1-igv.xml'
output[['Genotype [Factor]']] = 'wt'
input = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
input = EzDataset$new(file=input)

myApp = EzAppBowtie2$new()
myApp$run(input=input$copy()$subset(1), output=output, param=param)
