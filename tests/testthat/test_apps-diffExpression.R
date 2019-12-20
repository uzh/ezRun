context("differential expression apps with example data")

## this tests do take long therefore we only run them if the environment variable RUN_LONG_TEST is set to TRUE
# Sys.setenv(RUN_LONG_TEST=TRUE)

cwd = getwd()

testScratchDir = "/srv/GT/analysis/ezRunTestScratch"
testScratchDir = "~/tmp/ezRunTestScratch"


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
  param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/R64/Annotation/Release_98-2019-12-03'
  param[['refFeatureFile']] = 'genes.gtf'
  param[['featureLevel']] = 'gene'
  param[['grouping']] = 'Genotype'
  param[['grouping2']] = ''
  param[['sampleGroup']] = 'mut'
  param[['refGroup']] = 'wt'
  param[['runGO']] = 'false'
  param[['backgroundExpression']] = '10'
  param[['expressionName']] = ''
  param[['transcriptTypes']] = 'protein_coding'
  param[['specialOptions']] = ''
  param[['mail']] = ''
  param[['comparison']] = 'mut--over--wt'
  param[['name']] = 'mut--over--wt'
  param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
  param[['resultDir']] = 'p1001/Count_Result'
  param$linkHtmlLibDir = '' ## disables linking of libs
  return(param)
}

test_that("deseq2_withgo", {
  skipLong()
  ezSystem("rm -fr /scratch/test_deseq2_withgo/*")
  setwdNew("/scratch/test_deseq2_withgo")
  param = yeastCommonDiffExprParam()
  input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE),
                        dataRoot=param$dataRoot)
  output = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts_deseq2/dataset.tsv", package="ezRun", mustWork = TRUE),
                         dataRoot=param$dataRoot)
  param$runGO = FALSE
  myApp = EzAppDeseq2$new()
  myApp$run(input=input, output=output, param=param)
  setwd(cwd)
})

test_that("edger_withgo", {
  skipLong()
  ezSystem("rm -fr /scratch/test_edger_withgo/*")
  setwdNew("/scratch/test_edger_withgo")
  param = yeastCommonDiffExprParam()
  input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE),
                        dataRoot=param$dataRoot)
  output = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts_edger/dataset.tsv", package="ezRun", mustWork = TRUE),
                         dataRoot=param$dataRoot)
  param[['normMethod']] = 'TMM'
  myApp = EzAppEdger$new()
  myApp$run(input=input, output=output, param=param)
  setwd(cwd)
})

test_that("count_QC", {
  skipLong()
  testDir = file.path(testScratchDir, "countqc")
  unlink(testDir)
  setwdNew(testDir)
  param = yeastCommonDiffExprParam()
  #param$refAnnotationFile = file=system.file("extdata/genes_annotation.txt", package="ezRun", mustWork = TRUE)
  input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE),
                        dataRoot=param$dataRoot)
  output = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_countqc/dataset.tsv", package="ezRun", mustWork = TRUE),
                         dataRoot=param$dataRoot)
  param[['name']] = 'Count_QC'
  param[['normMethod']] = 'logMean'
  myApp = EzAppCountQC$new()
  myApp$run(input=input, output=output, param=param)
  setwd(cwd)
})

# test_that("edger_multi", {
#   skipLong()
#   setwdNew("/scratch/test_edger_multi")
#   input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE), dataRoot=param$dataRoot)
#   output = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts_edger/dataset.tsv", package="ezRun", mustWork = TRUE),
# dataRoot=param$dataRoot)
#   param = yeastCommonDiffExprParam()
#   param[['name']] = 'Edger_Multi'
#   param[['normMethod']] = 'TMM'
#   param[['mail']] = ''
#   myApp = EzAppEdgerMulti$new()
#   myApp$run(input=input, output=output, param=param)
#   setwd(cwd)
# })
