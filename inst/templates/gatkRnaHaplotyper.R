setwd("/scratch/gtan/dev/quickdev")
library(ezRun)
param = list()
param[['cores']] = '4'
param[['ram']] = '30'
param[['scratch']] = '500'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['samples']] = 'PID106_HSC,PID115star_HSC,PID282_HSC'
param[['paired']] = 'false'
param[['name']] = 'GATK_RnaVariants_p_HU_HSC'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2378/GatkRnaHaplotyper_37795_2019-07-03--12-36-29'
output = list()
output[['Name']] = 'GATK_RnaVariants_p_HU_HSC'
output[['VCF [File]']] = 'p2378/GatkRnaHaplotyper_37795_2019-07-03--12-36-29/GATK_RnaVariants_p_HU_HSC-haplo.vcf.gz'
output[['TBI [File]']] = 'p2378/GatkRnaHaplotyper_37795_2019-07-03--12-36-29/GATK_RnaVariants_p_HU_HSC-haplo.vcf.gz.tbi'
output[['Report [File]']] = 'p2378/GatkRnaHaplotyper_37795_2019-07-03--12-36-29/GATK_RnaVariants_p_HU_HSC'
output[['Html [Link]']] = 'p2378/GatkRnaHaplotyper_37795_2019-07-03--12-36-29/GATK_RnaVariants_p_HU_HSC/00index.html'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
input = '/srv/gstore/projects/p2378/GatkRnaHaplotyper_37795_2019-07-03--12-36-29/input_dataset.tsv'
# debug(ezMethodGatkRnaHaplotyper)
EzAppGatkRnaHaplotyper$new()$run(input=input, output=output, param=param)