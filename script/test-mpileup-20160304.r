
require(ezRun)
setwdNew("/scratch/test_mpileup")

param = list()
param[['cores']] = '8'
param[['ram']] = '10'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
param[['paired']] = 'true'
param[['strandMode']] = 'sense' ## not needed for DNA variants
param[['refFeatureFile']] = 'genes.gtf' ## not needed for DNA variants
param[['mail']] = ''
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['specialOptions']] = ''

param[['name']] = 'Mpileup_Variants'
param[['region']] = ''
param[['mpileupOptions']] = '--skip-indels --output-tags DP,DV,DPR,INFO/DPR,DP4,SP'
param[['callOptions']] = '--multiallelic-caller --keep-alts --variants-only'
param[['filterOptions']] = '--include "MIN(DP)>5"'
param[['resultDir']] = 'Variants'


input = EzDataset$new(file=system.file("extdata/yeast_10k_STAR/dataset.tsv", package="ezRun", mustWork = TRUE))

output = list()
output[['Name']] = 'Mpileup_Variants'
output[['VCF [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants.vcf.gz'
output[['TBI [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants.vcf.gz.tbi'
output[['IGV Starter [Link]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants-igv.jnlp'
output[['Report [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants'
output[['Html [Link]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants/00index.html'
output[['Species']] = ''
output[['refBuild']] = param$refBuild
output[['IGV Starter [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants-igv.jnlp'
output[['IGV Session [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants-igv.xml'
output = EzDataset$new(meta=output)

myApp = EzAppMpileup$new()
## run mpileup on a single file
myApp$run(input=input$copy()$subset(1), output=output$copy()$subset(1), param=param)

## run mpileup on all files together
myApp$run(input=input, output=output, param=param)
