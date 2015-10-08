
options(error=recover)
library(ezRun)
setwdNew("/scratch/fastqscreen_test")
param = list()
param[['cores']] = '8'
param[['ram']] = '40'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'false'
param[['name']] = 'FastqScreen_Result'
param[['nReads']] = '500000'
param[['nTopSpecies']] = '5'
param[['minAlignmentScore']] = '-20'
param[['confFile']] = 'variousSpecies_rRNA_20140901_silva119.conf'
param[['cmdOptions']] = '-k 10 --trim5 1 --trim3 4 --very-sensitive'
param[['mail']] = 'peter.schmid@ieu.uzh.ch'
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['resultDir']] = 'p1001/QC_FastqScreen_6576_FastQ_Plate_RNAs_JKG20_2015-06-30--13-15-52'
output = list()
output[['Name']] = 'FastqScreen_Result'
output[['Report [File]']] = 'p1001/QC_FastqScreen_6576_FastQ_Plate_RNAs_JKG20_2015-06-30--13-15-52/FastqScreen_Result'
output[['Html [Link]']] = 'p1001/QC_FastqScreen_6576_FastQ_Plate_RNAs_JKG20_2015-06-30--13-15-52/FastqScreen_Result/00index.html'
input = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)

myApp = EzAppFastqScreen$new()
myApp$run(input=input, output=output, param=param)
