setwd("/scratch/gtan/dev/quickdev")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '30'
param[['scratch']] = '180'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['samples']] = 'PID106_HSC,PID115star_HSC,PID282_HSC'
param[['paired']] = 'false'
param[['name']] = 'GATK_RnaVariants'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2378/GatkRnaHaplotyper_37795_2019-06-27--14-11-25'
output = list()
output[['Name']] = 'GATK_RnaVariants'
output[['VCF [File]']] = 'p2378/GatkRnaHaplotyper_37795_2019-06-27--14-11-25/GATK_RnaVariants-haplo.vcf.gz'
output[['TBI [File]']] = 'p2378/GatkRnaHaplotyper_37795_2019-06-27--14-11-25/GATK_RnaVariants-haplo.vcf.gz.tbi'
output[['Report [File]']] = 'p2378/GatkRnaHaplotyper_37795_2019-06-27--14-11-25/GATK_RnaVariants'
output[['Html [Link]']] = 'p2378/GatkRnaHaplotyper_37795_2019-06-27--14-11-25/GATK_RnaVariants/00index.html'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
input = '/srv/gstore/projects/p2378/GatkRnaHaplotyper_37795_2019-06-27--14-11-25/input_dataset.tsv'
debug(ezMethodGatkRnaHaplotyper)
EzAppGatkRnaHaplotyper$new()$run(input=input, output=output, param=param)

