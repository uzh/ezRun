#!/bin/bash

set -e
set -o pipefail

source /usr/local/ngseq/etc/lmod_profile
module add Tools/htslib Tools/samtools Tools/bcftools Tools/sambamba
module add QC/Flexbar QC/Trimmomatic
module add Aligner/STAR Aligner/Bowtie Aligner/BWA Aligner/Bowtie2
module add QC/FastQC QC/FastQScreen Tools/exceRpt

SCRIPT=$(readlink -f "$0")
SCRIPTPATH=$(dirname "$SCRIPT")
cd $SCRIPTPATH/..

R --vanilla --slave<< EOT
Sys.setenv(RUN_LONG_TEST=TRUE)
EZ_GLOBAL_VARIABLES <<- '/usr/local/ngseq/opt/EZ_GLOBAL_VARIABLES.txt'
library(ezRun)
.libPaths(c(".", .libPaths()))
.libPaths()
devtools::test()
EOT
