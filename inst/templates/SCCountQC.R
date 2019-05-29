## smart-Seq2 single cell
setwd("/export/local/scratch/gtan/dev/p2497-SCCountQC")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '30'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['name']] = 'SCCount_QC'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2497/SCCountQC_24763_2018-03-01--10-00-22'
output = list()
output[['Name']] = '20171222.A-SiCSeq_SCs_P5'
output[['Species']] = 'Mus musculus (house mouse)'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
output[['Static Report [Link]']] = 'p2497/SCCountQC_24763_2018-03-01--10-00-22/20171222.A-SiCSeq_SCs_P5_SCCountQC/00index.html'
output[['Report [File]']] = 'p2497/SCCountQC_24763_2018-03-01--10-00-22/20171222.A-SiCSeq_SCs_P5_SCCountQC'
input = list()
input[['Name']] = '20171222.A-SiCSeq_SCs_P5'
input[['Species']] = 'Mus musculus (house mouse)'
input[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
input[['paired']] = 'false'
input[['refFeatureFile']] = 'genes.gtf'
input[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
input[['CellDataset']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5-dataset.tsv'
input[['CountMatrix']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5-counts.mtx'
input[['Stats']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5-stats.txt'
input[['CellCyclePhase']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5-CellCyclePhase.txt'
input[['BAM']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5.bam'
input[['BAI']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5.bam.bai'
input[['PreprocessingLog']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5_preprocessing.log'
input[['STARLog']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5_STAR.log'
input[['featureLevel']] = 'gene'

#input = EzDataset$new(meta=input, dataRoot=param$dataRoot)
#param <- ezParam(param)
#output = EzDataset$new(meta=output, dataRoot=param$dataRoot)

#sce <- loadSCCountDataset(input, param)
debug(ezMethodSCCountQC)
#debug(txEndBias)
#debug(getRpkm)
#debug(loadSCCountDataset)
EzAppSCCountQC$new()$run(input=input, output=output, param=param)

## 10x single cell
# https://fgcz-sushi.uzh.ch/data_set/p2497/26599
setwd("/export/local/scratch/gtan/dev/p2497-SCCountQC")
param = list()
param[['cores']] = '8'
param[['ram']] = '30'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['name']] = 'SCCount_QC'
param[['reference']] = '10X_Ref_Mouse_GRCm38.p5_20180305_Release_91'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2497/SCCountQC_24763_2018-03-01--10-00-22'
input = list()
input[['Name']] = 'Mouse_p60_1'
input[['Species']] = 'Mus musculus (house mouse)'
input[['reference']] = '10X_Ref_Mouse_GRCm38.p5_20180305_Release_91'
input[['BAM']] = 'p2497/CellRangerCount_26591_2018-05-11--14-23-14/Mouse_p60_1/possorted_genome_bam.bam'
input[['BAI']] = 'p2497/CellRangerCount_26591_2018-05-11--14-23-14/Mouse_p60_1/possorted_genome_bam.bam.bai'
input[['CountMatrix']] = 'p2497/CellRangerCount_26591_2018-05-11--14-23-14/Mouse_p60_1/outs/filtered_gene_bc_matrices/10X_Ref_Mouse_GRCm38.p5_20180305_Release_91/matrix.mtx'
input[['featureLevel']] = 'gene'

scProtocol <- ifelse("STARLog" %in% input$colNames, "smart-Seq2", "10X")
param$scProtocol <- scProtocol

input = EzDataset$new(meta=input, dataRoot=param$dataRoot)
param <- ezParam(param)
sce <- loadSCCountDataset(input, param)

EzAppSCCountQC$new()$run(input=input, output=output, param=param)
