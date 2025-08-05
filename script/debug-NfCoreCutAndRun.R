#detach("package:ezRun", unload=TRUE)

library(ezRun, lib.loc = "/home/marconotaro/R/x86_64-pc-linux-gnu-library/4.5/")
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
source(EZ_GLOBAL_VARIABLES)

library(ezRun)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

setwdNew('/scratch/mnotaro/test-nfcore/cutandrun')

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
param[['peakCaller']] = 'macs2' ## can be seacr, macs2
param[['spikeinGenome']] = 'K12-MG1655' ## spike-in genome, defaulting to E. coli K12-MG1655 (def), for yeast set to R64-1-1, for fruit fly BDGP6
param[['normalization']] = 'Spikein' ## "spikein" (def), "RPKM", "CPM", "BPM", "None" 
param[['name']] <- 'NfCoreCutAndRun'
param[['grouping']] = 'Condition'
param[['control']] = 'Control'

## move to ruby
buildName = 'GRCm39' ## param$ezRef@refBuildName
blackPath = '/scratch/mnotaro/test-nfcore/.nextflow/assets/nf-core/cutandrun/assets/blacklists/' ## move to a global path
param[['blacklist']] <- case_when(
  buildName == 'GRCm39' ~ paste0(blackPath, buildName, '-blacklist.bed'),
  buildName == 'GRCm38' ~ paste0(blackPath, buildName, '-blacklist.bed'),
  buildName == 'GRCh38' ~ paste0(blackPath, buildName, '-blacklist.bed'),
  buildName == 'GRCh37' ~ paste0(blackPath, buildName, '-blacklist.bed'),
  TRUE ~ ""
)

input = '/srv/gstore/projects/p36614/o35568_NovaSeq_240807_X112_5M/dataset_cutandrun.tsv'  ## a control column was added

output = list()
output[['Name']] = 'NfCoreCutAndRun'
# output[['Result [File]']] = 'NfCoreCutAndRun'

debug(ezMethodNfCoreCutAndRun)

EzAppNfCoreCutAndRun$new()$run(input=input, output=output, param=param)


