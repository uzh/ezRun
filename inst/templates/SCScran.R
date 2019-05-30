# p2860 10X
setwd("/scratch/gtan/dev/SCScran-p2860")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '30'
param[['scratch']] = '50'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['samples']] = 'wt_4_F'
param[['name']] = 'SCScran'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['transcriptTypes']] = 'protein_coding'
param[['scProtocol']] = '10X'
param[['snnK']] = '15'
param[['visMethod']] = 'tSNE'
param[['knownMarkers']] = ''
param[['runPseudoTime']] = 'true'
param[['all2allMarkers']] = 'true'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2860/SCScran_2019-02-21--15-27-39'
output = list()
output[['Name']] = 'wt_4_F'
output[['Species']] = 'Mus musculus (house mouse)'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
output[['refFeatureFile']] = 'genes.gtf'
output[['Static Report [Link]']] = 'p2860/SCScran_2019-02-21--15-27-39/wt_4_F_SCran/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreScran_app/?data=p2860/2860/SCScran_2019-02-21--15-27-39/wt_4_F_SCran/sce-mwixochmusro.rds'
output[['Report [File]']] = 'p2860/SCScran_2019-02-21--15-27-39/wt_4_F_SCran'
output[['ResultDir [Link]']] = 'p2860/CellRangerCount_33994_2019-05-22--11-02-46/wt_4_F'
input = list()
input[['Name']] = 'wt_4_F'
input[['ResultDir']] = 'p2860/CellRangerCount_33994_2019-05-22--11-02-46/wt_4_F'
input[['Report']] = 'p2860/CellRangerCount_33994_2019-05-22--11-02-46/wt_4_F/web_summary.html'
input[['BAM']] = 'p2860/CellRangerCount_33994_2019-05-22--11-02-46/wt_4_F/possorted_genome_bam.bam'
input[['BAI']] = 'p2860/CellRangerCount_33994_2019-05-22--11-02-46/wt_4_F/possorted_genome_bam.bam.bai'
input[['Species']] = 'Mus musculus (house mouse)'
input[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
input[['refFeatureFile']] = 'genes.gtf'
input[['featureLevel']] = 'gene'
input[['CountMatrix']] = 'p2860/CellRangerCount_33994_2019-05-22--11-02-46/wt_4_F/filtered_feature_bc_matrix'
debug(ezMethodSCScran)
EzAppSCScran$new()$run(input=input, output=output, param=param)


