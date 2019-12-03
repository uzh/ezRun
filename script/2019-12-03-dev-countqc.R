

EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
library(ezRun)

setwdNew("/scratch/dev-countqc")



param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['process_mode']] = 'DATASET'
#param[['samples']] = 'WT_GM_S1,WT_GM_S2,WT_GM_S3,WT_GM_IL10_S5,WT_GM_IL10_S7,WT_GM_IL10_S8,WT_IL10_S9,WT_IL10_S11,WT_IL10_S12,WT_ctrl_S13,WT_ctrl_S15,WT_ctrl_S16,IRF5_GM_S18,IRF5_GM_S19,IRF5_GM_S20,IRF5_ctrl_S21,IRF5_ctrl_S23,IRF5_ctrl_S24,Csf2r_GM_S25,Csf2r_GM_S26,Csf2r_GM_S27,Csf2r_ctrl_S28,Csf2r_ctrl_S29,Csf2r_ctrl_S30'
param[['name']] = 'Count_QC'
param[['refBuild']] = 'Mus_musculus/GENCODE/GRCm38.p6/Annotation/Release_M23-2019-11-05'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['normMethod']] = 'logMean'
param[['runGO']] = 'false'
param[['backgroundExpression']] = '10'
param[['transcriptTypes']] = 'protein_coding'
param[['specialOptions']] = ''
param[['expressionName']] = ''
param[['mail']] = 'Hubert.Rehrauer@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2223/CountQC_41657_2019-11-21--15-58-29'
param[['isLastJob']] = TRUE
output = list()
output[['Name']] = 'Count_QC'
output[['Species']] = 'Mus musculus (house mouse)'
output[['refBuild']] = 'Mus_musculus/GENCODE/GRCm38.p6/Annotation/Release_M23-2019-11-05'
output[['Static Report [Link]']] = 'p2223/CountQC_41657_2019-11-21--15-58-29/Count_QC/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreCountQC_app/?data=p2223/CountQC_41657_2019-11-21--15-58-29/Count_QC/counts-myhzcyvdsofs-EzResult.RData'
output[['Report [File]']] = 'p2223/CountQC_41657_2019-11-21--15-58-29/Count_QC'
input = '/srv/gstore/projects/p2223/CountQC_41657_2019-11-21--15-58-29/input_dataset.tsv'
EzAppCountQC$new()$run(input=input, output=output, param=param)
