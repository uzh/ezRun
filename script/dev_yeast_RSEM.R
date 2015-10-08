
.libPaths("/srv/GT/analysis/course_sushi/lib")
library(ezRun)
setwdNew("/srv/GT/analysis/hubert/tmp/")
param = list()
param[['cores']] = '8'
param[['ram']] = '20'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
param[['paired']] = 'true'
param[['cmdOptions']] = '' # '--calc-ci'
param[['trimAdapter']] = 'false'
param[['trimLeft']] = '0'
param[['trimRight']] = '0'
param[['minTailQuality']] = '0'
param[['specialOptions']] = ''
param[['mail']] = 'Hubert.Rehrauer@fgcz.ethz.ch'
#param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1001/DUMMY'
#param = ezParam(param)
#param$ezRef
#ezParam(param)$paired
input = EzDataset(file=system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE))

#options(error=recover)
param[['dataRoot']] = system.file(package="ezRun", mustWork=TRUE)
outputMeta = input$meta[ , !input$columnHasTag("File")]
outputMeta[["Count [File]"]] = paste(param$resultDir, "/", rownames(outputMeta), ".txt", sep="")
outputMeta[["BAM [File]"]] = paste(param$resultDir, "/", rownames(outputMeta), ".bam", sep="")
outputMeta[["BAI [File]"]] = paste(param$resultDir, "/", rownames(outputMeta), ".bai", sep="")
outputMeta[["featureLevel"]] = "isoform"
output = EzDataset(meta=outputMeta)

# i = 1
# for (i in 1:nrow(input$meta)){
#   ezRunApp(ezAppRSEM, input=input[i, ], output=output[i, ], param=param)
# }


## run edgeR
param[['process_mode']] = 'DATASET'
param$cores = "1"
#param$refFeatureFile = "genes.gtf"
#param$strandMode = "both"
param$name = "mycomp"
param$expressionName = "transcriptCount"
param$dataRoot = getwd()
param$grouping = "Genotype"
param$sampleGroup = "mut"
param$refGroup = "wt"
param$comparisn = "mut--over-wt"
input = output
input$meta[["Count [File]"]] = basename(input$meta[["Count [File]"]])
#input$meta[["BAI [File]"]] = basename(input$meta[["BAI [File]"]])
output = list()
output[['Name']] = 'outputName'
output[['Report [File]']] = 'twoGroupsReportName'
output[['Html [Link]']] = 'twoGroupsReportName/00index.html'
output[['Species']] = ''
output[['refBuild']] = param$refBuild
output = EzDataset(meta=output)
options(error=recover)
ezRunApp(ezAppEdger, input=input, output=output, param=param)


