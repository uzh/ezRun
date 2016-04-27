
# clean up everyting, start new
rm(list = ls())
# flag whether or not to debug 
bDebug <- FALSE

# define global variables
EZ_GLOBAL_VARIABLES <<- '/home/petervr/myRepo/ezRun/inst/extdata/EZ_GLOBAL_VARIABLES.txt'
sPackWd <- "~/myRepo/ezRun"

### # setting working directory to local ezRun directory
setwd(sPackWd)
sRunMode <- "dev"
# load local version of ezRun or installed version
if (sRunMode == "dev") {
  devtools::load_all()
} else {
  library(ezRun)
}

# helper function

# change the working directory to where the results will be placed
sCurWd <- getwd()
sResultDir <- '/scratch/PVR_test/bamFiles_ExonCounting_DEXSeqResults'
### # create empty results dir
setwdNew(sResultDir)

# start constructing parameter
sRefFeatureFile <- "/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.PatchesSkipped/Annotation/Version-2015-06-25/Genes/genes.gtf"
refBuild = "Mus_musculus/Ensembl/GRCm38.PatchesSkipped/Annotation/Version-2015-06-25"
param <- ezParam(list(refBuild=refBuild))
# check whether feature file can be found
stopifnot(file.exists(param[['ezRef']]@refFeatureFile))

param[['cores']] = '1'
param[['ram']] = '16'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['mail']] = 'peter.vonrohr@gmail.com'
param[['dataRoot']] = ''
param[['grouping']] = 'Condition'
param[['sampleGroup']] = 'Colon'
param[['refGroup']] = 'SI'

# continue here  with testing both plots
param[['disp_plot']] = 'dispersion_estimate_plot'
param[['ma_plot']] = 'ma_plot'

# specify extension of counts files
param[['countfile_ext']] = 'count'
param[['countfile_path']] = sResultDir

# specify name of gff file
param[['gff_file']] = 'genes.gff'

# specify the name that will appear in the report
param[['name']] = 'Exon Usage Comparison in MM Colon Vs SI'

# specify output information, at least output[['Name']]  must be specified
output = list()
output[['Name']] = 'bamFiles_ExonCounting_DEXSeqResults'

# file with input meta data
input = '/scratch/PVR_test/bamFiles_ExonCounting/complete_dataset_local.tsv'


### # run the ezApp analysis
EzAppDEXSeqAnalysis$new()$run(input = input, output = output, param = param)

# reset wd
setwd(sCurWd)
