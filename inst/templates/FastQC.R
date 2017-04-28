setwd("/home/gtan/analysis/p1997-FastQCRmarkdown")
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
param = list()
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'true'
param[['name']] = 'FastQC_Result'
param[['cmdOptions']] = ''
param[['mail']] = ''
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
