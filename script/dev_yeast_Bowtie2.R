
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
param[['cmdOptions']] = '--no-unal'
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
outputMeta[["BAM [File]"]] = paste(param$resultDir, "/", rownames(outputMeta), ".bam", sep="")
outputMeta[["BAI [File]"]] = paste(param$resultDir, "/", rownames(outputMeta), ".bai", sep="")
outputMeta[["IGV Starter [File]"]] = paste(param$resultDir, "/", rownames(outputMeta), "-igv.jnlp", sep="")
outputMeta[["IGV Starter [Link]"]] = paste(param$resultDir, "/", rownames(outputMeta), "-igv.jnlp", sep="")
outputMeta[["IGV Session [File]"]] = paste(param$resultDir, "/", rownames(outputMeta), "-igv.xml", sep="")

output = EzDataset(meta=outputMeta)
i = 1
for (i in 1:nrow(input$meta)){
   ezRunApp(ezAppBowtie2, input=input[i, ], output=output[i, ], param=param)
}


## run bamstats
param$cores = "8"
param$refFeatureFile = "genes.gtf"
param$strandMode = "both"
param$name = "RNA_BAM_Stats"
param$dataRoot = getwd()
input = output
input$meta[["BAM [File]"]] = basename(input$meta[["BAM [File]"]])
input$meta[["BAI [File]"]] = basename(input$meta[["BAI [File]"]])
output = list()
output[['Name']] = 'RNA_BAM_Statistics'
output[['Report [File]']] = 'p1001/DUMMY/RNA_BAM_Statistics'
output[['Html [Link]']] = 'p1001/DUMMY/RNA_BAM_Statistics/00index.html'
output[['Species']] = ''
output[['build']] = param$build
output[['featureFile']] = 'genes.gtf'
output = EzDataset(meta=output)
options(error=recover)
ezRunApp(ezAppRnaBamStats, input=input, output=output, param=param)



