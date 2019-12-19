library(ezRun)
# Report
## input
input = '/home/miquel/dataset.tsv'
## output
output = list()
output[['Name']] = 'Count_QC'
output[['Species']] = ''
output[['refBuild']] = 'Homo_sapiens/UCSC/hg38'
output[['Static Report [Link]']] = 'p3093/Excerpt_38298_2019-12-19--09-55-33/Excerpt_Report/00index.html'
output[['Report [File]']] = 'p3093/Excerpt_38298_2019-12-19--09-55-33/Excerpt_Report'
## param
param=list()
param$name='Excerpt_Report'
param$process_mode='DATASET'
param$dataRoot='/srv/gstore/projects'

## run report app
EzAppExceRptReport$new()$run(input=input, output=output, param=param)