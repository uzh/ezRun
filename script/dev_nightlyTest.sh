#!/bin/bash
set -eu  ## also fail if there is an unset variable
set -o pipefail
set -x 
set -o history -o histexpand

#mailaddress="gxtx_data_mngt@fgcz.ethz.ch"
mailaddress="Hubert.Rehrauer@fgcz.ethz.ch"


## run with
##  /srv/GT/analysis/course_sushi_testCases/nightlyTest.sh > /srv/GT/analysis/course_sushi_testCases/nightlyTest.out 2> /srv/GT/analysis/course_sushi_testCases/nightlyTest.err


dieWithMail() {
        echo -e "The failed commandline is:\n $1" | mail -s "ERROR: sushi test run script failed" $mailaddress
        exit 1
}

sushiCmd="/usr/local/ngseq/opt/Ruby_Gems/ruby/1.9.1//bin/sushi_fabric -I /srv/GT/analysis/course_sushi"
readDataset="/srv/GT/analysis/course_sushi/public/gstore/projects/p1000/yeast_10k/dataset.tsv"
alignedDataset="/srv/GT/analysis/course_sushi/public/gstore/projects/p1000/yeast10k_STAR_622_2015-05-04--19-28-50/dataset.tsv"
countDataset1="/srv/GT/analysis/course_sushi/public/gstore/projects/p1000/yeast_10k_RSEM_2015-05-04--19-15-08/dataset.tsv"
countDataset2="/srv/GT/analysis/course_sushi/public/gstore/projects/p1000/yeast_10k_CountOverlaps_638_2015-05-04--19-34-02/dataset.tsv"
cd /srv/GT/analysis/course_sushi


$sushiCmd --class FastqcApp -p 1000 --run -d $readDataset --next_dataset_name testcase_fastqc || dieWithMail " !! !:* "


