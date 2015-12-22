context("differential expression apps with example data")

## this tests do take long therefore we only run them if the environment variable RUN_LONG_TEST is set to TRUE
# Sys.setenv(RUN_LONG_TEST=TRUE)

cwd = getwd()

skipLong = function(){
  if (Sys.getenv("RUN_LONG_TEST") == "TRUE"){
    return()
  } else {
    skip("not running lengthy tests")
  }
}

yeastCommonDiffExprParam = function(){
  param = list()
  param[['cores']] = '1'
  param[['ram']] = '2'
  param[['scratch']] = '10'
  param[['node']] = ''
  param[['process_mode']] = 'DATASET'
  param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
  param[['refFeatureFile']] = 'genes.gtf'
  param[['featureLevel']] = 'gene'
  param[['grouping']] = 'Genotype'
  param[['sampleGroup']] = 'mut'
  param[['refGroup']] = 'wt'
  param[['runGO']] = 'true'
  param[['expressionName']] = ''
  param[['specialOptions']] = ''
  param[['mail']] = ''
  param[['comparison']] = 'mut--over--wt'
  param[['name']] = 'mut--over--wt'
  param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
  param[['resultDir']] = 'p1001/Count_Result'
  return(param)
}

test_that("deseq2_withgo", {
  skipLong()
  ezSystem("rm -fr /scratch/test_deseq2_withgo/*")
  setwdNew("/scratch/test_deseq2_withgo")
  input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE))
  output = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts_deseq2/dataset.tsv", package="ezRun", mustWork = TRUE))
  param = yeastCommonDiffExprParam()
  myApp = EzAppDeseq2$new()
  myApp$run(input=input, output=output, param=param)
  setwd(cwd)
})

test_that("edger_withgo", {
  skipLong()
  ezSystem("rm -fr /scratch/test_edger_withgo/*")
  setwdNew("/scratch/test_edger_withgo")
  input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE))
  output = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts_edger/dataset.tsv", package="ezRun", mustWork = TRUE))
  param = yeastCommonDiffExprParam()
  param[['normMethod']] = 'TMM'
  myApp = EzAppEdger$new()
  myApp$run(input=input, output=output, param=param)
  setwd(cwd)
})

test_that("count_QC", {
  skipLong()
  ezSystem("rm -fr /scratch/test_count_QC/*")
  setwdNew("/scratch/test_count_QC")
  input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE))
  output = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts_edger/dataset.tsv", package="ezRun", mustWork = TRUE))
  param = yeastCommonDiffExprParam()
  param[['name']] = 'Count_QC'
  param[['normMethod']] = 'logMean'
  myApp = EzAppCountQC$new()
  myApp$run(input=input, output=output, param=param)
  setwd(cwd)
})

# test_that("edger_multi", {
#   skipLong()
#   setwdNew("/scratch/test_edger_multi")
#   input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE))
#   output = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts_edger/dataset.tsv", package="ezRun", mustWork = TRUE))
#   param = yeastCommonDiffExprParam()
#   param[['name']] = 'Edger_Multi'
#   param[['normMethod']] = 'TMM'
#   param[['mail']] = ''
#   myApp = EzAppEdgerMulti$new()
#   myApp$run(input=input, output=output, param=param)
#   setwd(cwd)
# })
