## p2284 10x
setwd("/scratch/gtan/dev/SCReports-p2284")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '100'
param[['scratch']] = '200'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['name']] = 'SCCount_QC'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = ''
param[['min_genes']] = '500'
param[['max_genes']] = '3000'
param[['min_counts']] = '50000'
param[['pcs']] = '14'
param[['markersToCheck']] = 'Zhengs lists=Psat1,Phgdh,Psph,Dnmt1,Dnmt3a,Hdac9,Vwf,Selp,Kdr,Ramp3,Slc6a6,Car4,Car8,Ankrd37,Rgcc,Podoxl,Tmem176a,Tmem176b,Lrg1,Esm1,Dll4,CD36,Lamb1,Apold1,Cxcl12,Lpl,Mb,Slc26a10,Hba-a1,Hba-a2B cells;=Cd79a,Ly6d,Cd79b,H2-DMb2,H2-Ob,Fcmr,Ccr7,Bank1,Cd55,Ms4a1;'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2284/SCCountQC_28929_2018-08-22--19-38-16'
output = list()
output[['Name']] = 'AVM_17_08_2018'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2284/SCCountQC_28929_2018-08-22--19-38-16/AVM_17_08_2018_SCCountQC/00index.html'
output[['Report [File]']] = 'p2284/SCCountQC_28929_2018-08-22--19-38-16/AVM_17_08_2018_SCCountQC'
input = list()
input[['Name']] = 'AVM_17_08_2018'
input[['ResultDir']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018'
input[['Report']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018/outs/web_summary.html'
input[['BAM']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018/outs/possorted_genome_bam.bam'
input[['BAI']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018/outs/possorted_genome_bam.bam.bai'
input[['Species']] = 'Homo sapiens (human)'
input[['CountMatrix']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018/outs/filtered_gene_bc_matrices'
input[['Sample Id']] = 'bfs_183566'
input[['Order Id']] = '4718'
input[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
input[['refFeatureFile']] = 'genes.gtf'
input[['featureLevel']] = 'gene'

# input = EzDataset$new(meta=input, dataRoot=param$dataRoot)
# param <- ezParam(param)
# param$scProtocol <- ifelse("STARLog" %in% input$colNames, "smart-Seq2", "10x")
# sce <- loadSCCountDataset(input, param)
#debug(ezMethodSCReport)
EzAppSCReport$new()$run(input=input, output=output, param=param)
