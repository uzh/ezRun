#!/bin/bash
#$ -hold_jid 537404,537405,537406,537407
set -e
set -o pipefail
umask 0002

#### SET THE STAGE
SCRATCH_DIR=/scratch/CountQC_42329_2019-12-17--15-00-27_temp$$
GSTORE_DIR=/srv/gstore/projects
INPUT_DATASET=/srv/gstore/projects/p1001/CountQC_42329_2019-12-17--15-00-27/input_dataset.tsv
LAST_JOB=TRUE
echo "Job runs on `hostname`"
echo "at $SCRATCH_DIR"
mkdir $SCRATCH_DIR || exit 1
cd $SCRATCH_DIR || exit 1
source /usr/local/ngseq/etc/lmod_profile
module add Dev/R/3.6.1

#### NOW THE ACTUAL JOBS STARTS
R --vanilla --slave<<  EOT
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
library(ezRun)
param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['node']] = 'fgcz-c-042,fgcz-c-046,fgcz-c-048,fgcz-c-056,fgcz-c-063,fgcz-c-065,fgcz-h-004,fgcz-h-006,fgcz-h-007,fgcz-h-008,fgcz-h-009,fgcz-h-010,fgcz-h-011,fgcz-h-012,fgcz-h-013,fgcz-h-014,fgcz-h-015,fgcz-h-016,fgcz-h-017,fgcz-h-018'
param[['process_mode']] = 'DATASET'
param[['samples']] = 'wt_1,wt_2,mut_1,mut_2'
param[['name']] = 'Count_QC'
param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/R64/Annotation/Release_98-2019-12-03'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['normMethod']] = 'logMean'
param[['runGO']] = 'true'
param[['backgroundExpression']] = '10'
param[['transcriptTypes']] = 'protein_coding'
param[['specialOptions']] = ''
param[['expressionName']] = ''
param[['mail']] = ''
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1001/CountQC_42329_2019-12-17--15-00-27'
param[['isLastJob']] = TRUE
output = list()
output[['Name']] = 'Count_QC'
output[['Species']] = 'S. cerevisiae'
output[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/R64/Annotation/Release_98-2019-12-03'
output[['Static Report [Link]']] = 'p1001/CountQC_42329_2019-12-17--15-00-27/Count_QC/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreCountQC_app/?data=p1001/CountQC_42329_2019-12-17--15-00-27/Count_QC/counts-qfzzemtkhcgo-EzResult.RData'
output[['Report [File]']] = 'p1001/CountQC_42329_2019-12-17--15-00-27/Count_QC'
input = '/srv/gstore/projects/p1001/CountQC_42329_2019-12-17--15-00-27/input_dataset.tsv'
EzAppCountQC\$new()\$run(input=input, output=output, param=param)
EOT


#### JOB IS DONE WE PUT THINGS IN PLACE AND CLEAN AUP
g-req -w copy Count_QC /srv/gstore/projects/p1001/CountQC_42329_2019-12-17--15-00-27
cd /scratch
rm -rf /scratch/CountQC_42329_2019-12-17--15-00-27_temp$$ || exit 1

