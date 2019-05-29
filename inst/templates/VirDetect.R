Sys.setenv("PATH"=paste("/usr/local/ngseq/bin/", "/usr/local/ngseq/packages/Dev/jdk/8/bin", "/usr/local/ngseq/packages/QC/Flexbar/3.0.3/bin", "/usr/local/ngseq/packages/Aligner/Bowtie2/2.3.2/bin", "/usr/local/ngseq/packages/Tools/samtools/1.5/bin", Sys.getenv("PATH"), sep=":"))

setwd("/srv/GT/analysis/gtan/debug/p2150-VirDetect")
library(ezRun)
param = list()
param[['hostBuild']] = 'Sus_scrofa/NAGRP/Sscrofa11.1'
param[['virBuild']] = 'Virome/refseq/vertebrate_20170831'
param[['paired']] = 'true'
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['process_mode']] = 'SAMPLE'
param[['cmdOptionsHost']] = '--very-sensitive'
param[['cmdOptions']] = '-a --very-sensitive --no-mixed --no-discordant -X 1000'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '5'
param[['trimRight']] = '0'
param[['minTailQuality']] = '10'
param[['minAvgQuality']] = '20'
param[['minReadLength']] = '50'
param[['minReadCount']] = '10'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
output = list()
output[['Name']] = '4-2'
output[['Species']] = 'Sus scrofa'
output[['OutReport']] = '4-2.html'

input = list()
input[['Name']] = '4-2'
input[['Condition']] = ''
input[['Read1']] = 'p2150/NextSeq500_20170831_NS121_o3707/20170831.A-4-2_R1.fastq.gz'
input[['Read2']] = 'p2150/NextSeq500_20170831_NS121_o3707/20170831.A-4-2_R2.fastq.gz'
input[['Species']] = 'Sus scrofa'
input[['FragmentSize']] = '0'
input[['SampleConc']] = '50'
input[['Tube']] = 'p2150_3707/23'
input[['Extract Id']] = 'bfe_61147'
input[['Index']] = 'CTGAAGCT-CAGGACGT'
input[['PlatePosition']] = 'ready-made by the user_'
input[['LibConc_100_800bp']] = '0'
input[['LibConc_qPCR']] = '0'
input[['Adapter1']] = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
input[['Adapter2']] = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
input[['strandMode']] = 'NA'
input[['LibraryPrepKit']] = 'NEB Next Ultra DNA'
input[['EnrichmentMethod']] = 'None'
input[['InputAmount']] = '200'
input[['Read Count']] = '360280'
#undebug(ezMethodVirDetect)
EzAppVirDetect$new()$run(input=input, output=output, param=param)
