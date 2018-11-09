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


## p2284 o4718 10X
setwd("/scratch/gtan/dev/SCReports-p2284")
library(ezRun)
param = list()
param[['cores']] = '4'
param[['ram']] = '8'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['samples']] = 'AVM_17_08_2018'
param[['name']] = 'SCReport'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = ''
param[['minGenesPerCell']] = '500'
param[['maxGenesPerCell']] = '3000'
param[['minReadsPerCell']] = '50000'
param[['pcs']] = '14'
param[['pcGenes']] = ''
param[['x.low.cutoff']] = '0.1'
param[['x.high.cutoff']] = '8'
param[['y.cutoff']] = '1'
param[['resolution']] = '0.6'
param[['markersToCheck']] = 'Zheng List=Psat1, Phgdh, Psph, Dnmt1, Dnmt3a, Hdac9, Vwf, Selp, Kdr, Ramp3, Slc6a6, Car4, Car8, Ankrd37, Rgcc, Podoxl, Tmem176a, Tmem176b, Lrg1, Esm1, Dll4, CD36, Lamb1, Apold1, Cxcl12, Lpl, Mb, Slc26a10, Hba-a1, Hba-a2;B cell=Cd79a,Ly6d,Cd79b,H2-DMb2,H2-Ob,Fcmr,Ccr7,Bank1,Cd55,Ms4a1;DC-like=Lgals3,Napsa,Plbd1,Plbd2,Ccr2,Rnase6,Plac8,Ifitm6,Naaa,Ear2,H2afy,Cd209a;Endothelial Cells=Ly6c1,Egfl7,Gpihbp1,Cdh5,Mgll,Slc9a3r2,Emcn,Kdr,Pecam1,Rgcc,Cdh5,Kdr,Flt1;Granulocytes=S100a8,S100a9,Slpi,Csf3r,Lmnb1,Retnlg,Clec4d,Hp,Hdc, S100a9;Macrophages=Dab2,Adgre1,Mgl2,Mrc1,Hpgd,P2ry6,C3ar1,F13a1,Maf,Ms4a7,Fcgr1;NK Cells=Ccl5,Nkg7,Klrk1,Klrd1,Ncr1,Ctsw,Klrb1c,Gzma,Gzmb;Pericytes=Kcnj8,Vtn,Colec11,Steap4,Abcc9,Myo1b,Cog7,P2ry14,Heyl,Gnb4,Pdgfrb;Schwann Cells=Plp1,Kcna1,Cnp,S100b,Gfra3,Gpr3,Nrn1,Aspa,Cd59a,Stmn1;Smooth Muscle Cells=PTagln,Mustn1,Myh11,Mylk,Pcp4l1,Sncg,Lmod1,Des,Pln,Nrp2,Acta2;T Cells=Cd3g,CD3d,Lat,Cd3e,Skap1,Il7r,Lef1,CD247,Tcf7,Itk;Fibroblasts=Lamc1,Pcsk6,Pdgfra,Entpd2,Dpep1,Adamts5,Medag,Ms4a4d,Lamb1, Tcf21,Dkk3,Tbx29,Wif1,Frzb,Meox1,Prg4,Abl3bp,Pdgfr,Mdk,Gstm5,Igfbp6,Col1a1,Tcf4;Muscle Stem Cells=Pax 7,Myf5,Myod,Myog;FAPs=Pdgfra,Sca1;'
param[['runPseudoTime']] = 'true'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2284/SCReport_28929_2018-11-06--00-11-59'
output = list()
output[['Name']] = 'AVM_17_08_2018'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2284/SCReport_28929_2018-11-06--00-11-59/AVM_17_08_2018_SCReport/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-176.uzh.ch/shiny/fgcz_exploreSingleCell_app/?data=p2284/SCReport_28929_2018-11-06--00-11-59/AVM_17_08_2018_SCReport/SCReport-nrtwjwznzizf.rds'
output[['Report [File]']] = 'p2284/SCReport_28929_2018-11-06--00-11-59/AVM_17_08_2018_SCReport'
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
#debug(seuratPreProcess)
EzAppSCReport$new()$run(input=input, output=output, param=param)

## p2284 o4866 10X multi plates/runs
setwd("/scratch/gtan/dev/SCReports-p2284")
library(ezRun)
param = list()
param[['cores']] = '4'
param[['ram']] = '8'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['samples']] = 'AVM_sample_26092018'
param[['name']] = 'SCReport'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = ''
param[['minGenesPerCell']] = '500'
param[['maxGenesPerCell']] = '3000'
param[['minReadsPerCell']] = '50000'
param[['pcs']] = '20'
param[['pcGenes']] = ''
param[['x.low.cutoff']] = '0.0125'
param[['x.high.cutoff']] = '3'
param[['y.cutoff']] = '0.5'
param[['vars.to.regress']] = 'nUMI,perc_mito'
param[['resolution']] = '0.4'
param[['markersToShow']] = '10'
param[['markersToCheck']] = ''
param[['runPseudoTime']] = 'true'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2284/SCReport_30486_AVM_sample_26092018_2018-11-07--15-36-46'
output = list()
output[['Name']] = 'AVM_sample_26092018'
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2284/SCReport_30486_AVM_sample_26092018_2018-11-07--15-36-46/AVM_sample_26092018_SCReport/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreSingleCell_app/?data=p2284/SCReport_30486_AVM_sample_26092018_2018-11-07--15-36-46/AVM_sample_26092018_SCReport/SCReport-avzlqrgakgyf.rds'
output[['Report [File]']] = 'p2284/SCReport_30486_AVM_sample_26092018_2018-11-07--15-36-46/AVM_sample_26092018_SCReport'
input = "/srv/gstore/projects/p2284/CellRangerCount_30477_NOV28_o4866_2018-10-11--09-10-56/dataset.tsv"

# input <- EzDataset$new(file=input, dataRoot=param$dataRoot)
# param <- ezParam(param)
# param$scProtocol <- "10x"
# debug(loadSCCountDataset)
# sce <- loadSCCountDataset(input, param)
# debug(seuratPreProcess)
EzAppSCReport$new()$run(input=input, output=output, param=param)
