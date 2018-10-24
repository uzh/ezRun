## p2284 10x
setwd("/scratch/gtan/dev/SCReports-p2284")

library(ezRun)
param = list()
param[['cores']] = '4'
param[['ram']] = '8'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['name']] = 'SCReport'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = ''
param[['minReadsPerGene']] = '3'
param[['minGenesPerCell']] = '500'
param[['maxGenesPerCell']] = '3000'
param[['minReadsPerCell']] = '50000'
param[['pcs']] = '14'
param[['markersToCheck']] = 'Zhengs lists=Psat1,Phgdh,Psph,Dnmt1,Dnmt3a,Hdac9,Vwf,Selp,Kdr,Ramp3,Slc6a6,Car4,Car8,Ankrd37,Rgcc,Podoxl,Tmem176a,Tmem176b,Lrg1,Esm1,Dll4,CD36,Lamb1,Apold1,Cxcl12,Lpl,Mb,Slc26a10,Hba-a1,Hba-a2B cells;=Cd79a,Ly6d,Cd79b,H2-DMb2,H2-Ob,Fcmr,Ccr7,Bank1,Cd55,Ms4a1;'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2284/SCReport_28929_2018-10-05--22-01-26'
output = list()
output[['Name']] = 'AVM_17_08_2018'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2284/SCReport_28929_2018-10-05--22-01-26/AVM_17_08_2018_SCReport/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-176.uzh.ch/shiny/fgcz_exploreSC_app/?data=p2284/SCReport_28929_2018-10-05--22-01-26/AVM_17_08_2018_SCReport/SCReport-butvkxmjqqvu.rds'
output[['Report [File]']] = 'p2284/SCReport_28929_2018-10-05--22-01-26/AVM_17_08_2018_SCReport'
input = list()
input[['Name']] = 'AVM_17_08_2018'
input[['ResultDir']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018'
input[['Report']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018/outs/web_summary.html'
input[['BAM']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018/outs/possorted_genome_bam.bam'
input[['BAI']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018/outs/possorted_genome_bam.bam.bai'
input[['Species']] = 'Homo sapiens (human)'
input[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
input[['CountMatrix']] = 'p2284/CellRangerCount_28912_NS214_o4718_2018-08-20--11-24-01/AVM_17_08_2018/outs/filtered_gene_bc_matrices'
input[['Sample Id']] = 'bfs_183566'
input[['Order Id']] = '4718'
input[['refFeatureFile']] = 'genes.gtf'
input[['featureLevel']] = 'gene'

# input = EzDataset$new(meta=input, dataRoot=param$dataRoot)
# param <- ezParam(param)
# param$scProtocol <- ifelse("STARLog" %in% input$colNames, "smart-Seq2", "10x")
# sce <- loadSCCountDataset(input, param)
# debug(seuratPreProcess)
EzAppSCReport$new()$run(input=input, output=output, param=param)
