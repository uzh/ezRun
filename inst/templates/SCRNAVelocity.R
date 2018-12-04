## p2838 smart-Seq2
setwd("/scratch/gtan/dev/RNAVelocity-p2838")

library(ezRun)
param = list()
param[['cores']] = '12'
param[['ram']] = '24'
param[['scratch']] = '150'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['name']] = 'SCRNAVelocity'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
param[['refFeatureFile']] = 'genes.gtf'
param[['scProtocol']] = 'smart-Seq2'
param[['markersToCheck']] = 'Calb1,Pvalb,Aqp2,Slc12a3,Trpv5'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'

input = list()
input[['Name']] = 'NCC_Trpv5_Tomato_sc_A01'
input[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
input[['refFeatureFile']] = 'genes.gtf'
input[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreSingleCell_app/?data=p2838/SCReport_31393_2018-11-18--17-12-56/SCReport/SCReport-cinggfxwqcul.rds'
input[['Report [File]']] = 'p2838/SCReport_31393_2018-11-18--17-12-56/SCReport'
input[['BAM']] = 'p2838/SCCountsApp_28469_2018-11-14--08-58-02/NCC_Trpv5_Tomato_sc_A01.bam'

output = list()
output[['Name']] = 'NCC_Trpv5_Tomato_sc_A01'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2838/SCRNAVelocity_31393_2018-11-18--17-12-56/SCRNAVelocity/00index.html'
output[['Report [File]']] = 'p2838/SCRNAVelocity_31393_2018-11-18--17-12-56/SCRNAVelocity'

# input = EzDataset$new(meta=input, dataRoot=param$dataRoot)
# param <- ezParam(param)

debug(ezMethodSCRNAVelocity)
EzAppSCRNAVelocity$new()$run(input=input, output=output, param=param)

