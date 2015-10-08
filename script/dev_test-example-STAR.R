
options(error=recover)
library(ezRun)
setwdNew("/scratch/STAR_test")
param = list()
param[['cores']] = '8'
param[['ram']] = '40'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Version-2013-03-07'
param[['paired']] = 'true'
param[['strandMode']] = 'sense'
param[['refFeatureFile']] = 'genes.gtf'
param[['cmdOptions']] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
param[['getChimericJunctions']] = 'false'
param[['trimAdapter']] = 'false'
param[['trimLeft']] = '0'
param[['trimRight']] = '0'
param[['minTailQuality']] = '0'
param[['specialOptions']] = ''
param[['mail']] = 'peter.schmid@ieu.uzh.ch'
param[['dataRoot']] = system.file(package="ezRun", mustWork = TRUE)
param[['resultDir']] = 'p1001/Map_STAR_5750_2015-09-14--14-08-53'
output = list()
output[['Name']] = 'wt_1'
output[['BAM [File]']] = 'p1001/Map_STAR_5750_2015-09-14--14-08-53/wt_1.bam'
output[['BAI [File]']] = 'p1001/Map_STAR_5750_2015-09-14--14-08-53/wt_1.bam.bai'
output[['IGV Starter [Link]']] = 'p1001/Map_STAR_5750_2015-09-14--14-08-53/wt_1-igv.jnlp'
output[['Species']] = 'S. cerevisiae'
output[['refBuild']] = 'Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Version-2013-03-07'
output[['paired']] = 'true'
output[['refFeatureFile']] = 'genes.gtf'
output[['strandMode']] = 'sense'
output[['Read Count']] = '9794'
output[['IGV Starter [File]']] = 'p1001/Map_STAR_5750_2015-09-14--14-08-53/wt_1-igv.jnlp'
output[['IGV Session [File]']] = 'p1001/Map_STAR_5750_2015-09-14--14-08-53/wt_1-igv.xml'
output[['Genotype [Factor]']] = 'wt'
input = system.file("extdata/yeast_10k/dataset.tsv", package="ezRun", mustWork = TRUE)
input = EzDataset$new(file=input)

myApp = EzAppSTAR$new()
myApp$run(input=input$copy()$subset(1), output=output, param=param)


