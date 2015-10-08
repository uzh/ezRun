
options(error=recover)
library(ezRun)
setwdNew("/scratch/bampreview_test")
param = list()
param[['cores']] = '8'
param[['ram']] = '50'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'false'
param[['name']] = 'BAM_Preview'
param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
param[['refFeatureFile']] = 'genes.gtf'
param[['strandMode']] = 'both'
param[['subsampleReads']] = 2
param[['mapMethod']] = 'STAR'
param[['mapOptions']] = ''
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '1'
param[['trimRight']] = '0'
param[['minTailQuality']] = '0'
param[['specialOptions']] = ''
param[['mail']] = 'peter.schmid@ieu.uzh.ch'
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['resultDir']] = 'p1001/QC_BAMPreview_736_2015-09-03--09-23-09'
output = list()
output[['Name']] = 'BAM_Preview'
output[['Report [File]']] = 'p1001/QC_BAMPreview_736_2015-09-03--09-23-09/BAM_Preview'
output[['Html [Link]']] = 'p1001/QC_BAMPreview_736_2015-09-03--09-23-09/BAM_Preview/00index.html'
output[['Species']] = ''
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.PatchesSkipped/Annotation/Version-2015-06-25'
output[['refFeatureFile']] = 'genes.gtf'
input = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
# input = EzDataset$new(file=input)
# input = ezMethodSubsampleReads(input=input, param=param)


myApp = EzAppBamPreview$new()
myApp$run(input=input, output=output, param=param)
