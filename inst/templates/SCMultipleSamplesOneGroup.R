library(ezRun)
source("~/workspaceR/dayme-scripts/sushi_scripts_mod/app-SCMultipleSamplesOneGroup.R")  #remove later
setwd("/scratch/daymegr") #remove later
param = list()
param[['cores']] = '4'
param[['ram']] = '50'
param[['scratch']] = '50'
param[['node']] = 'fgcz-c-065'
param[['process_mode']] = 'DATASET'
param[['samples']] = 'temporal_lobe_20190219,TLE'
param[['name']] = 'SCReportMerging'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['refFeatureFile']] = 'genes.gtf'
param[['npcs']] = '20'
param[['resolution']] = '0.6'
param[['batchCorrection']] = 'TRUE'
param[['chosenClusters']] = ''
param[['all2allMarkers']] = 'false'
param[['specialOptions']] = ''
param[['mail']] = 'daymegr@fgcz.uzh.ch'
param[['Rversion']] = 'Dev/R/3.6.0'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2979/SCReportMerging_37761_2019-06-25--15-34-23'
output = list()
output[['Name']] = 'SCReportMerging'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2979/SCReportMerging_37761_2019-06-25--15-34-23/SCReportMerging/00index.html'
output[['Report [File]']] = 'p2979/SCReportMerging_37761_2019-06-25--15-34-23/SCReportMerging'
input = '/srv/gstore/projects/p2979/SCReportMerging_37761_2019-06-25--15-34-23/input_dataset.tsv'
#debug(ezMethodSCReportMerging)
EzAppSCMultipleSamplesOneGroup$new()$run(input=input, output=output, param=param)

input = EzDataset$new(file=input, dataRoot=param$dataRoot)
output <- EzDataset$new(meta=output, dataRoot=param$dataRoot)

appDefaults = rbind(scProtocol=ezFrame(npcs=ezFrame(Type="numeric", 
                                                    DefaultValue=20, 
                                                    Description="The maximal dimensions to use for reduction"),
                                       resolution=ezFrame(Type="numeric", 
                                                          DefaultValue=0.6,
                                                          Description="Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities."),
                                       batchCorrection=ezFrame(Type="character", 
                                                               DefaultValue="CCA",
                                                               Description="Which batch correction method to use? None or CCA"),
                                       chosenClusters=ezFrame(Type="charList",
                                                              DefaultValue="",
                                                              Description="The clusters to choose from each sample.In the format of sample1=cluster1,cluster2;sample2=cluster1,cluster2."),
                                       all2allMarkers=ezFrame(Type="logical",
                                                              DefaultValue=FALSE, 
                                                              Description="Run all against all cluster comparisons?"),
                                       markersToShow=ezFrame(Type="numeric", 
                                                             DefaultValue=10, 
                                                             Description="The markers to show in the heatmap of cluster marker genes"),
                                       maxSamplesSupported=ezFrame(Type="numeric", 
                                                                   DefaultValue=5, 
                                                                   Description="Maximum number of samples to compare")))
param <- ezParam(param, appDefaults = appDefaults)
debug(ezMethodSCReportMerging)
ezMethodSCReportMerging(input=input, output=output, param=param)
