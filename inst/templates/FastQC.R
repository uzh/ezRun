
## p1997 Paired-end samples
setwd("/srv/GT/analysis/gtan/p1997-FastQCRmarkdown")
param = list()
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'true'
param[['name']] = 'FastQC_Result'
param[['cmdOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1997/Fastqc_17423_2017-04-06--11-46-10'
output = list()
output[['Name']] = 'FastQC_Result'
output[['Report [File]']] = 'p1997/Fastqc_17423_2017-04-06--11-46-10/FastQC_Result'
output[['Html [Link]']] = 'p1997/Fastqc_17423_2017-04-06--11-46-10/FastQC_Result/00index.html'
input = '/srv/gstore/projects/p1997/Fastqc_17423_2017-04-06--11-46-10/input_dataset.tsv'
EzAppFastqc$new()$run(input=input, output=output, param=param)

#rmarkdown::render(input="/home/gtan/Repos/FGCZ/ezRun/inst/templates/FastQC.Rmd",
#                  output_dir=".", output_file="00index.html")

## p2401 Single-end samples and with correlation
setwd("/srv/GT/analysis/gtan/p2401-FastQCRmarkdown")
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'false'
param[['name']] = 'FastQC_Result'
param[['cmdOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2401/Fastqc_17776_2017-04-24--09-09-19'
output = list()
output[['Name']] = 'FastQC_Result'
output[['Report [File]']] = 'p2401/Fastqc_17776_2017-04-24--09-09-19/FastQC_Result'
output[['Html [Link]']] = 'p2401/Fastqc_17776_2017-04-24--09-09-19/FastQC_Result/00index.html'
input = '/srv/gstore/projects/p2401/Fastqc_17776_2017-04-24--09-09-19/input_dataset5.tsv'
EzAppFastqc$new()$run(input=input, output=output, param=param)
rmarkdown::render(input="/Users/gtan/Repos/FGCZ/ezRun/inst/templates/FastQC.Rmd", envir = new.env(),
                  output_dir="FastQC_Result", output_file="00index.html")

## p2438 with plate position
setwd("/srv/GT/analysis/gtan/p2438-FastQC")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'false'
param[['name']] = 'FastQC_Result'
param[['cmdOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2438/Fastqc_18564_2017-06-07--14-06-33'
output = list()
output[['Name']] = 'FastQC_Result'
output[['Report [File]']] = 'p2438/Fastqc_18564_2017-06-07--14-06-33/FastQC_Result'
output[['Html [Link]']] = 'p2438/Fastqc_18564_2017-06-07--14-06-33/FastQC_Result/00index.html'
input = '/srv/gstore/projects/p2438/Fastqc_18564_2017-06-07--14-06-33/input_dataset.tsv'
EzAppFastqc$new()$run(input=input, output=output, param=param)

# p2000 single cell bam
setwd("/scratch/gtan/p2000-SCFastQC")
library(ezRun)
setEnvironments("fastqc")
param = list()
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'false'
param[['name']] = 'FastQC_Result'
param[['cmdOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'projects/p2000/SingleCellUnmappedBam'
output = list()
output[['Name']] = 'FastQC_Result'
output[['Report [File]']] = 'p2438/Fastqc_18564_2017-06-07--14-06-33/FastQC_Result'
output[['Html [Link]']] = 'p2438/Fastqc_18564_2017-06-07--14-06-33/FastQC_Result/00index.html'
input = '/srv/gstore/projects/p2000/SingleCellUnmappedBam/dataset.tsv'
#input = 'dataset.tsv'
debug(ezMethodFastQC)
EzAppFastqc$new()$run(input=input, output=output, param=param)

# p3221 RNA-seq
setwd("/scratch/gtan/dev/quickdev")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['samples']] = ''
param[['paired']] = 'false'
param[['perLibrary']] = 'true'
param[['name']] = 'FastQC_Result'
param[['cmdOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p3221/Fastqc_40422_2019-10-17--09-45-22'
output = list()
output[['Name']] = 'FastQC_Result'
output[['Report [File]']] = 'p3221/Fastqc_40422_2019-10-17--09-45-22/FastQC_Result'
output[['Html [Link]']] = 'p3221/Fastqc_40422_2019-10-17--09-45-22/FastQC_Result/00index.html'
output[['MultiQC [Link]']] = 'p3221/Fastqc_40422_2019-10-17--09-45-22/FastQC_Result/multiqc_report.html'
input = '/srv/gstore/projects/p3221/Fastqc_40422_2019-10-17--09-45-22/input_dataset.tsv'
EzAppFastqc$new()$run(input=input, output=output, param=param)