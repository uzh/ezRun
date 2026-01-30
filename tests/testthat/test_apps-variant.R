context("Variant apps with example data")

## this tests do take long therefore we only run them if the environment variable RUN_LONG_TEST is set to TRUE
# Sys.setenv(RUN_LONG_TEST=TRUE)

cwd = getwd()

skipLong = function() {
  if (Sys.getenv("RUN_LONG_TEST") == "TRUE") {
    return()
  } else {
    skip("not running lengthy tests")
  }
}

yeastCommonVariantParam = function() {
  param = list()
  param[['cores']] = '8'
  param[['ram']] = '10'
  param[['scratch']] = '10'
  param[['node']] = ''
  param[['process_mode']] = 'DATASET'
  param[[
    'refBuild'
  ]] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
  param[['paired']] = 'true'
  param[['strandMode']] = 'sense'
  param[['refFeatureFile']] = 'genes.gtf'
  param[['mail']] = ''
  param[['dataRoot']] = system.file(package = "ezRun", mustWork = TRUE)
  param[['resultDir']] = 'p1001/Count_Result'
  return(param)
}

test_that("Variant_Mpileup", {
  skipLong()
  ezSystem("rm -fr /scratch/test_mpileup/*")
  setwdNew("/scratch/test_mpileup")
  param = yeastCommonVariantParam()
  input = EzDataset$new(
    file = system.file(
      "extdata/yeast_10k_STAR/dataset.tsv",
      package = "ezRun",
      mustWork = TRUE
    ),
    dataRoot = param$dataRoot
  )
  output = list()
  output[['Name']] = 'Mpileup_Variants'
  output[[
    'VCF [File]'
  ]] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants.vcf.gz'
  output[[
    'TBI [File]'
  ]] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants.vcf.gz.tbi'
  output[[
    'IGV Starter [Link]'
  ]] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants-igv.jnlp'
  output[[
    'Report [File]'
  ]] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants'
  output[[
    'Html [Link]'
  ]] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants/00index.html'
  output[['Species']] = ''
  output[[
    'refBuild'
  ]] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
  output[[
    'IGV Starter [File]'
  ]] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants-igv.jnlp'
  output[[
    'IGV Session [File]'
  ]] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants-igv.xml'
  output = EzDataset$new(meta = output, dataRoot = param$dataRoot)
  param[['name']] = 'Mpileup_Variants'
  param[['region']] = ''
  param[[
    'mpileupOptions'
  ]] = '--skip-indels --output-tags DP,DV,DPR,INFO/DPR,DP4,SP'
  param[['callOptions']] = '--multiallelic-caller --keep-alts --variants-only'
  param[['filterOptions']] = '--include "MIN(DP)>5"'
  param[['specialOptions']] = ''
  param[['dataRoot']] = system.file(package = "ezRun", mustWork = TRUE)
  param[[
    'resultDir'
  ]] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54'
  myApp = EzAppMpileup$new()
  myApp$run(
    input = input$copy()$subset(1),
    output = output$copy()$subset(1),
    param = param
  )
  setwd(cwd)
})
