## input, output and param of yet to repair apps

### EzAppMpileup
param = list()
param[['cores']] = '8'
param[['ram']] = '30'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['paired']] = 'true'
param[['name']] = 'Mpileup_Variants'
param[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
param[['region']] = ''
param[['mpileupOptions']] = '--skip-indels --output-tags DP,DV,DPR,INFO/DPR,DP4,SP'
param[['callOptions']] = '--multiallelic-caller --keep-alts --variants-only'
param[['filterOptions']] = '--include "MIN(DP)>5"'
param[['specialOptions']] = ''
param[['mail']] = ''
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54'
output = list()
output[['Name']] = 'Mpileup_Variants'
output[['VCF [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants.vcf.gz'
output[['TBI [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants.vcf.gz.tbi'
output[['IGV Starter [Link]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants-igv.jnlp'
output[['Report [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants'
output[['Html [Link]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants/00index.html'
output[['Species']] = ''
output[['refBuild']] = 'Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Version-2013-03-18'
output[['IGV Starter [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants-igv.jnlp'
output[['IGV Session [File]']] = 'p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/Mpileup_Variants-igv.xml'
input = '/srv/gstore/projects/p1001/Variant_Analysis_samtoolsmpileup_8617_2015-11-12--12-19-54/input_dataset.tsv'


