
setwd("/scratch/gtan/p2578-atacENCODE")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '32'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['Species']] = 'Homo sapiens (human)'
param[['paired']] = 'true'
param[['name']] = 'asthmatic_P'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2578/AtacENCODE_18564_2017-06-07--14-06-33'
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
output = list()
output[['Name']] = 'AtacENCODE_Result'
output[['Report [File]']] = 'p2578/AtacENCODE_18564_2017-06-07--14-06-33/asthmatic_P'
output[['Html [Link]']] = 'p2578/AtacENCODE_18564_2017-06-07--14-06-33/asthmatic_P/out/asthmatic_P_report.html'
input = 'HiSeq2500_ATAC_o3757_asthmaticP.tsv'
EzAppAtacENCODE$new()$run(input=input, output=output, param=param)

## atacBamFilter and MACS2
library(ezRun)
setwd("/scratch/gtan/p2578-atacENCODE/atacBamFilter")
param = list()
param[['cores']] = '8'
param[['ram']] = '60'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
param[['paired']] = 'true'
param[['useControl']] = 'false'
param[['shiftATAC']] = 'false'
param[['refFeatureFile']] = 'genes.gtf'
param[['cmdOptions']] = '--nomodel --extsize 147 -g hs --bw 200'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2578/MACS2_22192_2018-01-24--14-31-29'
param[['mode']] = 'ATAC-seq'

input = list()
input[['Name']] = 'A9502US'
input[['BAM']] = 'p2578/Bowtie2_22155_2017-11-06--22-58-25/A9502US.bam'
input[['BAI']] = 'p2578/Bowtie2_22155_2017-11-06--22-58-25/A9502US.bam.bai'
input[['IGV Starter']] = 'p2578/Bowtie2_22155_2017-11-06--22-58-25/A9502US-igv.jnlp'
input[['Species']] = 'Homo sapiens (human)'
input[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
input[['paired']] = 'true'
input[['Read Count']] = '37857601'
input[['IGV Session']] = 'p2578/Bowtie2_22155_2017-11-06--22-58-25/A9502US-igv.xml'
input[['PreprocessingLog']] = 'p2578/Bowtie2_22155_2017-11-06--22-58-25/A9502US_preprocessing.log'
input[['Bowtie2Log']] = 'p2578/Bowtie2_22155_2017-11-06--22-58-25/A9502US_bowtie2.log'
input[['Condition']] = ''
input[['FragmentSize']] = '0'
input[['SampleConc']] = '85'
input[['Tube']] = 'p2578_3757/3'
input[['Index']] = 'AGGCAGAA'
input[['PlatePosition']] = 'usermade_caf_'
input[['LibConc_100_800bp']] = '0'
input[['LibConc_qPCR']] = '0'
input[['InputAmount']] = '0'
input[['Sample Id']] = 'bfs_164950'

output = list()
output[['Name']] = 'A9502US'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
output[['refFeatureFile']] = 'genes.gtf'
output[['paired']] = 'true'
output[['CalledPeaks [File]']] = 'p2578/MACS2_22192_2018-01-24--14-31-29/A9502US_peaks.xls'
output[['BED [File]']] = 'p2578/MACS2_22192_2018-01-24--14-31-29/A9502US_peaks.bed'
output[['PeakSequences [File]']] = 'p2578/MACS2_22192_2018-01-24--14-31-29/A9502US_peaks.fa'
output[['BigWigFile [File]']] = 'p2578/MACS2_22192_2018-01-24--14-31-29/A9502US_processed.bw'
output[['BAM [File]']] = 'p2578/Bowtie2_22155_2017-11-06--22-58-25/A9502US_processed.bam'
output[['BAI [File]']] = 'p2578/Bowtie2_22155_2017-11-06--22-58-25/A9502US_processed.bam.bai'
output[['Condition [Factor]']] = ''
output[['Sample Id [B-Fabric]']] = 'bfs_164950'
output[['FragmentSize [Characteristic]']] = '0'
output[['SampleConc [Characteristic]']] = '85'
output[['Tube [Characteristic]']] = 'p2578_3757/3'
output[['Index [Characteristic]']] = 'AGGCAGAA'
output[['PlatePosition [Characteristic]']] = 'usermade_caf_'
output[['LibConc_100_800bp [Characteristic]']] = '0'
output[['LibConc_qPCR [Characteristic]']] = '0'
output[['InputAmount [Characteristic]']] = '0'

## atacBamProcess
input = EzDataset$new(meta=input, dataRoot=param$dataRoot)
output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
param = ezParam(param)
system.time(output <- atacBamProcess(input=input, output=NA, param=param))

## annotatePeaks
annotatePeaks(peakFile="A9502US_peaks.xls", param=param)

## bam2bw
bam2bw(file="A9502US_processed.bam", destination="A9502US_processed_deepTpools.bw",
       paired=TRUE, method="deepTools", cores=8)

## MACS2
## module load Dev/Python2
## module load Tools/MACS2
#debug(atacBamProcess)
EzAppMacs2$new()$run(input=input, output=output, param=param)

## HOMER
### makeTagDirectory App
library(ezRun)
setwd("/scratch/gtan/p2578-atacENCODE/HOMER")
param = list()
param[['cores']] = '4'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['refBuild']] = 'Homo_sapiens/GENCODE/GRCh38.p10/Annotation/Release_27-2018-02-01'
param[['paired']] = 'true'
param[['refFeatureFile']] = 'genes.gtf'
param[['grouping']] = 'Condition'
param[['sampleGroup']] = 'asthmatic_P'
param[['refGroup']] = 'asthmatic_US'
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2578/HomerDiffPeaks_2018'
### special params for HOMER
param[['refBuildHOMER']] = 'hg38'
param[['repFoldChange']] = '2'
param[['repFDR']] = '0.1'
param[['balanced']] = 'true'
param[['style']] = 'histone' # factor, tss, groseq, dnase, super, mC

input <- "input_dataset.tsv"
output = list()
output[['Name']] = 'P--over--US'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/GENCODE/GRCh38.p10/Annotation/Release_27-2018-02-01'
output[['Report [File]']] = 'p2578/HomerDiffPeaks_2018/HomerDiffPeaks_asthmatic_P--over--US'
output[['DiffPeaks [File]']] = 'p2578/HomerDiffPeaks_2018/HomerDiffPeaks_asthmatic_P--over--US/diffPeaks.txt'

# input = EzDataset$new(file=input, dataRoot=param$dataRoot)
# param <- ezParam(param)
# output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
EzAppHomerDiffPeaks$new()$run(input=input, output=output, param=param)

## ataqv
library(ezRun)
setEnvironments("ataqv")
setwd("/export/local/scratch/gtan/p2578-atacENCODE/ataqv")
param = list()
param[['cores']] = '8'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['refBuild']] = 'Homo_sapiens/GENCODE/GRCh38.p10/Annotation/Release_27-2018-02-01'
param[['paired']] = 'true'
param[['refFeatureFile']] = 'genes.gtf'
param[['grouping']] = 'Condition'
param[['sampleGroup']] = 'asthmatic_P'
param[['refGroup']] = 'asthmatic_US'
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2578/ataqvApp-2018'
param[['name']] = 'ataqv'
### special params for HOMER
param[['refBuildHOMER']] = 'hg38'
param[['repFoldChange']] = '2'
param[['repFDR']] = '0.1'
param[['balanced']] = 'true'
param[['style']] = 'histone' # factor, tss, groseq, dnase, super, mC

input <- "input_dataset.tsv"
output = list()
output[['Name']] = 'ataqv'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/GENCODE/GRCh38.p10/Annotation/Release_27-2018-02-01'
output[['Report [File]']] = 'p2578/ataqvApp-2018/ataqv'
output[['Static Report [Link]']] = 'p2578/ataqvApp-2018/ataqv/index.html'

#input = EzDataset$new(file=input, dataRoot=param$dataRoot)
#param <- ezParam(param)
#output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
EzAppAtaqv$new()$run(input=input, output=output, param=param)

