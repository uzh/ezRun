
setwd("/scratch/gtan/p2578-atacENCODE")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '32'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['Species']] = 'Homo sapiens (human)'
param[['paired']] = 'true'
param[['name']] = 'asthmatic_P'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2578/AtacENCODE_18564_2017-06-07--14-06-33'
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
output = list()
output[['Name']] = 'AtacENCODE_Result'
output[['Report [File]']] = 'p2578/AtacENCODE_18564_2017-06-07--14-06-33/asthmatic_P'
output[['Html [Link]']] = 'p2578/AtacENCODE_18564_2017-06-07--14-06-33/asthmatic_P/out/asthmatic_P_report.html'
input = 'HiSeq2500_ATAC_o3757_asthmaticP.tsv'
EzAppAtacENCODE$new()$run(input=input, output=output, param=param)


## atacBamFilter
library(ezRun)
setwd("/scratch/gtan/p2578-atacENCODE/atacBamFilter")
localBamFile <- "/srv/gstore/projects/p2578/Bowtie2_22155_2017-11-06--22-58-25/A8901US.bam"
param <- list()
param$paired <- TRUE
param$cores <- 8L
input <- EzDataset(file="dataset.tsv", dataRoot="/srv/gstore/projects")
system.time(output <- atacBamProcess(input=input, output=NA, param=param))
