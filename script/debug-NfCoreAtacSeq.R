#detach("package:ezRun", unload=TRUE)

library(ezRun, lib.loc = "/home/marconotaro/R/x86_64-pc-linux-gnu-library/4.5/")
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
source(EZ_GLOBAL_VARIABLES)

library(ezRun)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

setwdNew('/scratch/mnotaro/test-nfcore/atacseq')

param = list()
param[['cores']] = '8'
param[['ram']] = '10'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'dataset'
param[['refBuild']] = 'Mus_musculus/GENCODE/GRCm39/Annotation/Release_M31-2023-01-30' ## root: /srv/GT/reference
param[['specialOptions']] = ''
param[['mail']] = 'marco.notaro@fgcz.uzh.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['peakStyle']] = 'broad'
param[['varStabilizationMethod']] = 'vst'
param[['grouping']] = 'Condition'
param[['name']] <- 'NfCoreAtacSeq'
param[['runTwoGroupAnalysis']] <- 'true'
# param[['resultDir']] <- 'p36614/o35568_NfCoreAtacSeq_2025-06-18--13-21-07'
  
input = '/srv/gstore/projects/p36614/o35568_NovaSeq_240807_X112_5M/dataset.tsv'
# input = '/srv/gstore/projects/p35006/o37938_NovaSeq_250325_X242/dataset.tsv' ## dataset-nocond.tsv

output = list()
output[['Name']] = 'o37938_NfCoreCutAndRun'
output[['Result [File]']] = 'p36614/o35568_NfCoreAtacSeq_2025-06-18--13-21-07/'

debug(ezMethodNfCoreAtacSeq)

EzAppNfCoreAtacSeq$new()$run(input=input, output=output, param=param)


