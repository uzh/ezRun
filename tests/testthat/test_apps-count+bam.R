context("Counting and bam apps with example data")

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

yeastCommonCountParam = function(){
  param = list()
  param[['cores']] = '1'
  param[['ram']] = '10'
  param[['scratch']] = '10'
  param[['node']] = ''
  param[['process_mode']] = 'SAMPLE'
  param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/R64/Annotation/Release_98-2019-12-03'
  param[['paired']] = 'true'
  param[['strandMode']] = 'sense'
  param[['refFeatureFile']] = 'genes.gtf'
  param[['featureLevel']] = 'gene'
  param[['specialOptions']] = ''
  param[['mail']] = ''
  param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
  param[['resultDir']] = 'p1001/Count_Result'
  return(param)
}

test_that("Count_RSEM", {
  skipLong()
  ezSystem("rm -fr /scratch/test_rsem/*")
  setwdNew("/scratch/test_rsem")
  param = yeastCommonCountParam()
  input = EzDataset$new(file=system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE),
                        dataRoot=param$dataRoot)
  output = list()
  output[['Name']] = 'wt_1'
  output[['Count [File]']] = 'p1001/Count_RSEM/wt_1.txt'
  output[['Species']] = 'S. cerevisiae'
  output[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
  output[['featureLevel']] = 'isoform'
  output[['refFeatureFile']] = 'genes.gtf'
  output[['strandMode']] = 'sense'
  output[['paired']] = 'true'
  output[['Read Count']] = '9794'
  output[['Genotype [Factor]']] = 'wt'
  param[['trimAdapter']] = 'false'
  param[['trimLeft']] = '0'
  param[['trimRight']] = '0'
  param[['minTailQuality']] = '0'
  param[['bowtie-e']] = '200'
  param[['cmdOptions']] = ' --calc-ci '
  param[['keepBam']] = 'false'
  param[['minAvgQuality']] = '10'
  param[['transcriptFasta']] = ''
  myApp = EzAppRSEM$new()
  myApp$run(input=input$copy()$subset(1), output=EzDataset$new(metaNew=output,dataRoot=param$dataRoot), param=param)
  setwd(cwd)
})

test_that("Count_FeatureCounts", {
  skipLong()
  ezSystem("rm -fr /scratch/test_featureCounts/*")
  setwdNew("/scratch/test_featureCounts")
  param = yeastCommonCountParam()
  input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR/dataset.tsv", package="ezRun", mustWork = TRUE),
                        dataRoot=param$dataRoot)
  output = EzDataset$new(file=system.file("extdata/yeast_10k_STAR_counts/dataset.tsv", package="ezRun", mustWork = TRUE),
                         dataRoot=param$dataRoot)
  param[['gtfFeatureType']] = 'exon'
  param[['allowMultiOverlap']] = 'true'
  param[['countPrimaryAlignmentsOnly']] = 'true'
  param[['minFeatureOverlap']] = '10'
  param[['minMapQuality']] = '10'
  param[['keepMultiHits']] = 'true'
  myApp = EzAppFeatureCounts$new()
  myApp$run(input=input$copy()$subset(1), output=output$copy()$subset(1), param=param)
  setwd(cwd)
})

test_that("RNA_Bamstats", {
  skipLong()
  ezSystem("rm -fr /scratch/test_RNA_Bamstats/*")
  setwdNew("/scratch/test_RNA_Bamstats")
  param = yeastCommonCountParam()
  input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR/dataset.tsv", package="ezRun", mustWork = TRUE),
                        dataRoot=param$dataRoot)
  output = list()
  output[['Name']] = 'RNA_BAM_Statistics'
  output[['Report [File]']] = 'p1001/QC_RNABamStats_8617_2015-11-17--10-15-39/RNA_BAM_Statistics'
  output[['Html [Link]']] = 'p1001/QC_RNABamStats_8617_2015-11-17--10-15-39/RNA_BAM_Statistics/00index.html'
  output[['Species']] = ''
  output[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
  output[['refFeatureFile']] = 'genes.gtf'
  param[['process_mode']] = 'DATASET'
  param[['name']] = 'RNA_BAM_Statistics'
  myApp = EzAppRnaBamStats$new()
  myApp$run(input=input, output=output, param=param)
  setwd(cwd)
})

test_that("TEQC", {
  skipLong()
  ezSystem("rm -fr /scratch/test_TEQC/*")
  setwdNew("/scratch/test_TEQC")
  param = yeastCommonCountParam()
  input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR/dataset.tsv", package="ezRun", mustWork=TRUE),
                        dataRoot=param$dataRoot)
  output = list()
  output[['Name']] = 'TEQC_Result'
  output[['Report [File]']] = 'p1001/QC_Teqc_5579_2015-05-04--13-41-58/TEQC_Result'
  output[['Html [Link]']] = 'p1001/QC_Teqc_5579_2015-05-04--13-41-58/TEQC_Result/00index.html'
  param[['cores']] = 1
  param[['paired']] = "false"
  param[['process_mode']] = 'DATASET'
  param[['name']] = 'TEQC_Result'
  param[['designFile']] = system.file("extdata/genes.bed", package="ezRun", mustWork=TRUE)
  param[['covUniformityPlot']] = 'true'
  param[['covTargetLengthPlot']] = 'true'
  param[['duplicatesPlot']] = 'true'
  param[['cmdOptions']] = ''
  myApp = EzAppTeqc$new()
  myApp$run(input=input, output=output, param=param)
  setwd(cwd)
})
