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
param[['refBuild']] = 'Homo_sapiens/GENCODE/GRCh38.p14/Annotation/Release_48-2025-07-03' ## root: /srv/GT/reference
param[['specialOptions']] = ''
param[['mail']] = 'marco.notaro@fgcz.uzh.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['peakStyle']] = 'broad'
param[['peakCaller']] = 'macs2' ## can be seacr, macs2
param[['spikeinGenome']] = 'K12-MG1655' ## spike-in genome, defaulting to E. coli K12-MG1655 (def), for yeast set to R64-1-1, for fruit fly BDGP6
param[['normalization']] = 'Spikein' ## "spikein" (def), "RPKM", "CPM", "BPM", "None" 
param[['name']] <- 'NfCoreCutAndRun'
param[['grouping']] = 'Condition'
param[['controlColumn']] = 'Control'

# input = '/srv/gstore/projects/p36614/o35568_NovaSeq_240807_X112_5M/dataset_cutandrun.tsv'  ## a control column was added
input = '/srv/gstore/projects/p35864/o36419_NovaSeq_241023_X157/dataset_cutandrun.tsv'

output = list()
output[['Name']] = 'NfCoreCutAndRun'
# output[['Result [File]']] = 'NfCoreCutAndRun'

debug(ezMethodNfCoreCutAndRun)

EzAppNfCoreCutAndRun$new()$run(input=input, output=output, param=param)
