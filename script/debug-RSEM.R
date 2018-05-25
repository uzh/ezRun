setwd("/scratch/gtan/debug/p1688-RSEM")
Sys.setenv("PATH"=paste("/usr/local/ngseq/packages/Dev/jdk/8/bin:/usr/local/ngseq/packages/QC/Flexbar/3.0.3/bin:/usr/local/ngseq/packages/Aligner/RSEM/1.3.0/bin:/usr/local/ngseq/packages/Aligner/STAR/2.5.3a/bin:/usr/local/ngseq/packages/Aligner/Bowtie2/2.3.2/bin:/usr/local/ngseq/packages/Aligner/Bowtie/1.1.2/bin:/usr/local/ngseq/packages/Tools/samtools/1.3.1/bin", Sys.getenv("PATH"), sep=":"))
Sys.setenv("Trimmomatic_jar"="/usr/local/ngseq/packages/QC/Trimmomatic/0.36/trimmomatic-0.36.jar")

library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '30'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
param[['paired']] = 'true'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['bowtie-e']] = '200'
param[['cmdOptions']] = ' --calc-ci --sort-bam-by-read-name'
param[['keepBam']] = 'false'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '0'
param[['trimRight']] = '0'
param[['minTailQuality']] = '0'
param[['minAvgQuality']] = '20'
param[['specialOptions']] = ''
param[['transcriptFasta']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1688/RSEM_20178_2017-08-16--13-24-26'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA,long_noncoding,short_noncoding,pseudogene'
output = list()
output[['Name']] = 's28A24Veh'
output[['Count [File]']] = 'p1688/RSEM_20178_2017-08-16--13-24-26/s28A24Veh.txt'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
output[['featureLevel']] = 'isoform'
output[['refFeatureFile']] = 'genes.gtf'
output[['strandMode']] = 'both'
output[['paired']] = 'true'
output[['Read Count']] = '6575082'
output[['PreprocessingLog [File]']] = 'p1688/RSEM_20178_2017-08-16--13-24-26/s28A24Veh_preprocessing.log'
output[['Conditon [Factor]']] = 'asthmatic'
output[['Treatment [Factor]']] = 'Vehicle'
input = list()
input[['Name']] = 's28A24Veh'
input[['Conditon']] = 'asthmatic'
input[['Treatment']] = 'Vehicle'
input[['Read1']] = 'p1688/GSE61141-fastq/GSE61141-fastq/s28A24Veh_1.fastq.gz'
input[['Read2']] = 'p1688/GSE61141-fastq/GSE61141-fastq/s28A24Veh_2.fastq.gz'
input[['Species']] = 'Homo sapiens (human)'
input[['strandMode']] = 'unstranded'
input[['Read Count']] = '6575082'
#debug(ezMethodRSEM)
EzAppRSEM$new()$run(input=input, output=output, param=param)

