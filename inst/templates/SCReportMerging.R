# p2529 merge 
setwd("/scratch/gtan/dev/SCReportMerging-p2529")
library(ezRun)
param = list()
param[['cores']] = '4'
param[['ram']] = '30'
param[['scratch']] = '50'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['samples']] = c('FAPS4_all', 'FAPS30_all_E1')
param[['name']] = 'SCReport'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = ''
param[['resolution']] = '0.6'
param[['chosenClusters1']] =''
param[['chosenClusters2']] ='0,1,2,3,4,5,6,8,9'
param[['markersToShow']] = '10'
param[['markersToCheck']] = 'Zheng List=Psat1, Phgdh, Psph, Dnmt1, Dnmt3a, Hdac9, Vwf, Selp, Kdr, Ramp3, Slc6a6, Car4, Car8, Ankrd37, Rgcc, Podoxl, Tmem176a, Tmem176b, Lrg1, Esm1, Dll4, CD36, Lamb1, Apold1, Cxcl12, Lpl, Mb, Slc26a10, Hba-a1, Hba-a2;B cell=Cd79a,Ly6d,Cd79b,H2-DMb2,H2-Ob,Fcmr,Ccr7,Bank1,Cd55,Ms4a1;DC-like=Lgals3,Napsa,Plbd1,Plbd2,Ccr2,Rnase6,Plac8,Ifitm6,Naaa,Ear2,H2afy,Cd209a;Endothelial Cells=Ly6c1,Egfl7,Gpihbp1,Cdh5,Mgll,Slc9a3r2,Emcn,Kdr,Pecam1,Rgcc,Cdh5,Kdr,Flt1;Granulocytes=S100a8,S100a9,Slpi,Csf3r,Lmnb1,Retnlg,Clec4d,Hp,Hdc, S100a9;Macrophages=Dab2,Adgre1,Mgl2,Mrc1,Hpgd,P2ry6,C3ar1,F13a1,Maf,Ms4a7,Fcgr1;NK Cells=Ccl5,Nkg7,Klrk1,Klrd1,Ncr1,Ctsw,Klrb1c,Gzma,Gzmb;Pericytes=Kcnj8,Vtn,Colec11,Steap4,Abcc9,Myo1b,Cog7,P2ry14,Heyl,Gnb4,Pdgfrb;Schwann Cells=Plp1,Kcna1,Cnp,S100b,Gfra3,Gpr3,Nrn1,Aspa,Cd59a,Stmn1;Smooth Muscle Cells=PTagln,Mustn1,Myh11,Mylk,Pcp4l1,Sncg,Lmod1,Des,Pln,Nrp2,Acta2;T Cells=Cd3g,CD3d,Lat,Cd3e,Skap1,Il7r,Lef1,CD247,Tcf7,Itk;Fibroblasts=Lamc1,Pcsk6,Pdgfra,Entpd2,Dpep1,Adamts5,Medag,Ms4a4d,Lamb1, Tcf21,Dkk3,Tbx29,Wif1,Frzb,Meox1,Prg4,Abl3bp,Pdgfr,Mdk,Gstm5,Igfbp6,Col1a1,Tcf4;Muscle Stem Cells=Pax 7,Myf5,Myod,Myog;FAPs=Pdgfra,Sca1,Ly6a;'
param[['runPseudoTime']] = 'true'
param[['all2allMarkers']] = 'true'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2529/SCReport_32851_FAPS4_all_2019-01-08--10-47-27'
output = list()
output[['Name']] = 'FAPS4_all'
output[['Species']] = 'Mus musculus (house mouse)'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2529/SCReport_32851_FAPS4_all_2019-01-08--10-47-27/FAPS4_all_SCReport/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreSingleCell_app/?data=p2529/SCReport_32851_FAPS4_all_2019-01-08--10-47-27/FAPS4_all_SCReport/SCReport-emituitwiiek.rds'
output[['Report [File]']] = 'p2529/SCReport_32851_FAPS4_all_2019-01-08--10-47-27/FAPS4_all_SCReportMerging'
output[['ResultDir [Link]']] = 'p2529/CellRangerCount_25012_2019-01-07--16-20-19/FAPS4_all'
input <- "/scratch/gtan/dev/SCReportMerging-p2529/SCReport_FAPS4_FAPS30.tsv"

# input = EzDataset$new(file=input, dataRoot=param$dataRoot)
# output <- EzDataset$new(meta=output, dataRoot=param$dataRoot)
# param <- ezParam(param)
debug(ezMethodSCReportMerging)
EzAppSCReportMerging$new()$run(input=input, output=output, param=param)
