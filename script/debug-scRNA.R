# R manual approach
setwd("/srv/local/scratch/gtan/smartSeq2")

refFn <- "/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Sequence/WholeGenomeFasta/genome.fa"

fastqFn <- "test_R1.fastq.gz"
bamFn <- "test_R1"
fastqI1Fn <- "test_I1.fastq.gz"
fastqI2Fn <- "test_I2.fastq.gz"

fastqFn <- "20170822.B-Pool1_EBO6_11_A2B5_R1.fastq.gz"
bamFn <- "20170822.B-Pool1_EBO6_11_A2B5_R1"
fastqI1Fn <- "20170822.B-Pool1_EBO6_11_A2B5_I1.fastq.gz"
fastqI2Fn <- "20170822.B-Pool1_EBO6_11_A2B5_I2.fastq.gz"

system.time(fastq2bam(fastqFn, refFn, bamFn, fastqI1Fn=fastqI1Fn, fastqI2Fn=fastqI2Fn))
#user    system   elapsed on GT/analysis
#7641.624   484.496 13885.729
#user    system   elapsed  on scratch
#7957.988   399.820 12582.806

# Read group approach
setwd("/scratch/gtan/scRNA")
fastqFns <- list.files(path="/srv/gstore/projects/p2288/HiSeq2500_20171011_RR99_o3511/dmx",
                       pattern="\\.fastq\\.gz$", full.names=TRUE)
fastqs2bam(fastqFns, bamFn="20171011.A-C1_HT_24H.bam")

setEnvironments("star")
setEnvironments("flexbar")
setEnvironments("trimmomatic")
setEnvironments("python2")
setEnvironments("samtools")
setwd("/scratch/gtan/scRNA")
param = list()
param[['cores']] = '8'
param[['ram']] = '40'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
param[['paired']] = 'false'
param[['strandMode']] = 'sense'
param[['refFeatureFile']] = 'genes.gtf'
param[['cmdOptions']] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
param[['getChimericJunctions']] = 'false'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '5'
param[['trimRight']] = '0'
param[['minTailQuality']] = '10'
param[['minAvgQuality']] = '10'
param[['minReadLength']] = '20'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '.'
param[['resultDir']] = 'p2438/STAR_18564_2017-06-12--13-46-30'
param[['twopassMode']] = "true" #false
output = list()
output[['Name']] = '20171011.A-C1_HT_24H'
output[['BAM [File]']] = '20171011.A-C1_HT_24H.bam'
output[['BAI [File]']] = '20171011.A-C1_HT_24H.bam.bai'
output[['IGV Starter [Link]']] = '20171011.A-C1_HT_24H-igv.jnlp'
output[['Species']] = 'Mus musculus (mouse)'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
output[['paired']] = 'false'
output[['refFeatureFile']] = 'genes.gtf'
output[['strandMode']] = 'sense'
output[['Read Count']] = '39534111'
output[['IGV Starter [File]']] = '20171011.A-C1_HT_24H-igv.jnlp'
output[['IGV Session [File]']] = '20171011.A-C1_HT_24H-igv.xml'
output[['PreprocessingLog [File]']] = '20171011.A-C1_HT_24H_preprocessing.log'
output[['STARLog [File]']] = '20171011.A-C1_HT_24H_STAR.log'
input = list()
input[['Name']] = '20171011.A-C1_HT_24H'
input[['Condition']] = ''
input[['Read1']] = '20171011.A-C1_HT_24H_unmapped.bam'
input[['Species']] = 'Mus musculus (mouse)'
input[['FragmentSize']] = '0'
input[['SampleConc']] = '9.1'
input[['Tube']] = 'p2438_3365/24'
input[['Extract Id']] = 'bfe_53715'
input[['Index']] = 'GAATTCGT-CAGGACGT'
input[['PlatePosition']] = '1_G6'
input[['LibConc_100_800bp']] = '1.75'
input[['LibConc_qPCR']] = '0'
input[['Adapter1']] = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
input[['strandMode']] = 'both'
input[['LibraryPrepKit']] = 'TruSeq RNA stranded'
input[['EnrichmentMethod']] = 'PolyA, TruSeq RNA Kit'
input[['InputAmount']] = '100'
input[['Read Count']] = '39534111'
debug(ezMethodSTAR)
EzAppSTAR$new()$run(input=input, output=output, param=param)

# p2277 HiSeq4000_20170822_RUN376_o3645
# fastqs2bam: is supposed to run on dmx server
# For each lane, it produces one unmapped Bam file and a dataset.tsv
setwd("/scratch/gtan/scRNA-p2277/HiSeq4000_20170822_RUN376_o3645")
library(ezRun)
input <- EzDataset(file="/srv/gstore/projects/p2277/HiSeq4000_20170822_RUN376_o3645_A2B5/A2B5-cells-dataset.tsv", dataRoot="/srv/gstore/projects")
fastqs2bam(fastqFns=input$getFullPaths("Read1"), readGroupNames=input$getNames(),
           bamFn="A2B5_unmapped.bam")
countMeta = input$meta[ , !input$columnHasTag("File")]
ezWrite.table(countMeta, file="A2B5-dataset.tsv", head='Name')

input <- EzDataset(file="/srv/gstore/projects/p2277/HiSeq4000_20170822_RUN376_o3645_Pool2/Auto-cells-dataset.tsv", dataRoot="/srv/gstore/projects")
fastqs2bam(fastqFns=input$getFullPaths("Read1"), readGroupNames=input$getNames(), 
           bamFn="Auto_unmapped.bam")
countMeta = input$meta[ , !input$columnHasTag("File")]
ezWrite.table(countMeta, file="Auto-dataset.tsv", head='Name')

input <- EzDataset(file="/srv/gstore/projects/p2277/HiSeq4000_20170822_RUN376_o3645_Pool3/NEG-cells-dataset.tsv", dataRoot="/srv/gstore/projects")
fastqs2bam(fastqFns=input$getFullPaths("Read1"), readGroupNames=input$getNames(), 
           bamFn="NEG_unmapped.bam")
countMeta = input$meta[ , !input$columnHasTag("File")]
ezWrite.table(countMeta, file="NEG-dataset.tsv", head='Name')

## STAR
setwd("/scratch/gtan/scRNA-p2277/STAR")
setEnvironments("star")
setEnvironments("flexbar")
setEnvironments("trimmomatic")
setEnvironments("python2")
setEnvironments("samtools")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '40'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['cmdOptions']] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
param[['getChimericJunctions']] = 'false'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '5'
param[['trimRight']] = '0'
param[['minTailQuality']] = '10'
param[['minAvgQuality']] = '10'
param[['minReadLength']] = '20'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/scratch/gtan'
param[['resultDir']] = 'scRNA-p2277/STAR'
output = list()
output[['Name']] = 'A2B5'
output[['BAM [File]']] = 'scRNA-p2277/STAR/A2B5.bam'
output[['BAI [File]']] = 'scRNA-p2277/STAR/A2B5.bam.bai'
output[['IGV Starter [Link]']] = 'scRNA-p2277/STAR/A2B5-igv.jnlp'
output[['Species']] = 'Mus musculus'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
output[['paired']] = 'false'
output[['refFeatureFile']] = 'genes.gtf'
output[['strandMode']] = 'antisense'
#output[['Read Count']] = '32231106'
output[['IGV Starter [File]']] = 'scRNA-p2277/STAR/A2B5-igv.jnlp'
output[['IGV Session [File]']] = 'scRNA-p2277/STAR/A2B5-igv.xml'
output[['PreprocessingLog [File]']] = 'scRNA-p2277/STAR/A2B5_preprocessing.log'
output[['STARLog [File]']] = 'scRNA-p2277/STAR/A2B5_STAR.log'
input = list()
input[['Name']] = 'A2B5'
input[['Read1']] = 'scRNA-p2277/HiSeq4000_20170822_RUN376_o3645/A2B5_unmapped_part.bam'
input[['Species']] = 'Mus musculus'
input[['strandMode']] = 'both'
#debug(ezMethodSingleCellSTAR)
EzAppSingleCellSTAR$new()$run(input=input, output=output, param=param)
input[['Name']] = 'Auto'
input[['Read1']] = 'scRNA-p2277/HiSeq4000_20170822_RUN376_o3645/Auto_unmapped_part.bam'
output[['BAM [File]']] = 'scRNA-p2277/STAR/Auto.bam'
output[['BAI [File]']] = 'scRNA-p2277/STAR/Auto.bam.bai'
output[['IGV Starter [Link]']] = 'scRNA-p2277/STAR/Auto-igv.jnlp'
output[['IGV Starter [File]']] = 'scRNA-p2277/STAR/Auto-igv.jnlp'
output[['IGV Session [File]']] = 'scRNA-p2277/STAR/Auto-igv.xml'
output[['PreprocessingLog [File]']] = 'scRNA-p2277/STAR/Auto_preprocessing.log'
output[['STARLog [File]']] = 'scRNA-p2277/STAR/Auto_STAR.log'
EzAppSingleCellSTAR$new()$run(input=input, output=output, param=param)
input[['Name']] = 'NEG'
input[['Read1']] = 'scRNA-p2277/HiSeq4000_20170822_RUN376_o3645/NEG_unmapped_part.bam'
output[['BAM [File]']] = 'scRNA-p2277/STAR/NEG.bam'
output[['BAI [File]']] = 'scRNA-p2277/STAR/NEG.bam.bai'
output[['IGV Starter [Link]']] = 'scRNA-p2277/STAR/NEG-igv.jnlp'
output[['IGV Starter [File]']] = 'scRNA-p2277/STAR/NEG-igv.jnlp'
output[['IGV Session [File]']] = 'scRNA-p2277/STAR/NEG-igv.xml'
output[['PreprocessingLog [File]']] = 'scRNA-p2277/STAR/NEG_preprocessing.log'
output[['STARLog [File]']] = 'scRNA-p2277/STAR/NEG_STAR.log'
EzAppSingleCellSTAR$new()$run(input=input, output=output, param=param)


## FeatureCounts
setwd("/scratch/gtan/scRNA-p2277/FeatureCounts")
library(ezRun)
param = list()
param[['cores']] = '8'
param[['ram']] = '20'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['gtfFeatureType']] = 'exon'
param[['allowMultiOverlap']] = 'true'
param[['countPrimaryAlignmentsOnly']] = 'true'
param[['minFeatureOverlap']] = '10'
param[['minMapQuality']] = '10'
param[['keepMultiHits']] = 'true'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA,long_noncoding,short_noncoding,pseudogene'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/scratch/gtan'
param[['resultDir']] = 'scRNA-p2277/FeatureCounts'
output = list()
output[['Name']] = 'A2B5'
output[['Count [File]']] = 'scRNA-p2277/FeatureCounts/A2B5.mtx'
output[['Stats [File]']] = 'scRNA-p2277/FeatureCounts/A2B5-stats.txt'
output[['CellCyclePhase [File]']] = 'scRNA-p2277/FeatureCounts/A2B5-CellCyclePhase.txt'
output[['Species']] = 'Mus musculus'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
output[['featureLevel']] = 'gene'
output[['refFeatureFile']] = 'genes.gtf'
output[['strandMode']] = 'both'
output[['paired']] = 'false'
output[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA,long_noncoding,short_noncoding,pseudogene'
input = list()
input[['Name']] = 'A2B5'
input[['BAM']] = 'scRNA-p2277/STAR/A2B5.bam'
input[['BAI']] = 'scRNA-p2277/STAR/A2B5.bam.bai'
input[['Species']] = 'Mus musculus'
input[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
input[['paired']] = 'false'
input[['refFeatureFile']] = 'genes.gtf'
input[['strandMode']] = 'both'
#debug(ezMethodSingleCellFeatureCounts)
EzAppSingleCellFeatureCounts$new()$run(input=input, output=output, param=param)

input[['Name']] = 'Auto'
input[['BAM']] = 'scRNA-p2277/STAR/Auto.bam'
input[['BAI']] = 'scRNA-p2277/STAR/Auto.bam.bai'
output[['Name']] = 'Auto'
output[['Count [File]']] = 'scRNA-p2277/FeatureCounts/Auto.mtx'
output[['Stats [File]']] = 'scRNA-p2277/FeatureCounts/Auto-stats.txt'
output[['CellCyclePhase [File]']] = 'scRNA-p2277/FeatureCounts/Auto-CellCyclePhase.txt'
EzAppSingleCellFeatureCounts$new()$run(input=input, output=output, param=param)

input[['Name']] = 'NEG'
input[['BAM']] = 'scRNA-p2277/STAR/NEG.bam'
input[['BAI']] = 'scRNA-p2277/STAR/NEG.bam.bai'
output[['Name']] = 'NEG'
output[['Count [File]']] = 'scRNA-p2277/FeatureCounts/NEG.mtx'
output[['Stats [File]']] = 'scRNA-p2277/FeatureCounts/NEG-stats.txt'
output[['CellCyclePhase [File]']] = 'scRNA-p2277/FeatureCounts/NEG-CellCyclePhase.txt'
EzAppSingleCellFeatureCounts$new()$run(input=input, output=output, param=param)

## EzAppSCCounts: smart-seq2 with smaller one uBam
setwd("/scratch/gtan/dev/scRNA-p2277")
library(ezRun)
setEnvironments("star")
setEnvironments("flexbar")
setEnvironments("trimmomatic")
setEnvironments("python2")
setEnvironments("samtools")
setEnvironments("picard")
param = list()
param[['cores']] = '8'
param[['ram']] = '40'
param[['scratch']] = '100'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['cmdOptions']] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
param[['getChimericJunctions']] = 'false'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '5'
param[['trimRight']] = '0'
param[['minTailQuality']] = '10'
param[['minAvgQuality']] = '10'
param[['minReadLength']] = '20'
param[['controlSeqs']] = 'TdTomato_55420622,Cre_NC_005856.1'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/scratch/gtan/dev'
param[['resultDir']] = 'scRNA-p2277/scCount'
output = list()
output[['Species']] = 'Mus musculus'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
output[['paired']] = 'false'
output[['refFeatureFile']] = 'genes.gtf'
output[['strandMode']] = 'antisense'
input = list()
input[['Species']] = 'Mus musculus'
input[['strandMode']] = 'both'
## A2B5
input[['Name']] = 'A2B5'
input[['Read1']] = 'scRNA-p2277/HiSeq4000_20170822_RUN376_o3645/A2B5_unmapped_part.bam'
input[['CellDataset']] = 'scRNA-p2277/HiSeq4000_20170822_RUN376_o3645/A2B5-dataset.tsv'
output[['Name']] = 'A2B5'
output[['ResultDir']] = 'scRNA-p2277/scCount'
output[['BAM [File]']] = 'scRNA-p2277/scCount/A2B5.bam'
output[['BAI [File]']] = 'scRNA-p2277/scCount/A2B5.bam.bai'
output[['PreprocessingLog [File]']] = 'scRNA-p2277/scCount/A2B5_preprocessing.log'
output[['STARLog [File]']] = 'scRNA-p2277/scCount/A2B5_STAR.log'
output[['CountMatrix [File]']] = 'scRNA-p2277/scCount/A2B5-counts.mtx'
output[['Stats [File]']] = 'scRNA-p2277/scCount/A2B5-stats.txt'
output[['CellCyclePhase [File]']] = 'scRNA-p2277/scCount/A2B5-CellCyclePhase.txt'
output[['CellDataset [File]']] = 'scRNA-p2277/scCount/A2B5-dataset.txt'
# debug(ezMethodSingleCellSTAR)
# debug(ezMethodSingleCellFeatureCounts)
EzAppSCCounts$new()$run(input=input, output=output, param=param)
## Auto
input[['Name']] = 'Auto'
input[['Read1']] = 'scRNA-p2277/HiSeq4000_20170822_RUN376_o3645/Auto_unmapped_part.bam'
input[['CellDataset']] = 'scRNA-p2277/HiSeq4000_20170822_RUN376_o3645/Auto-dataset.tsv'
output[['Name']] = 'Auto'
output[['BAM [File]']] = 'scRNA-p2277/scCount/Auto.bam'
output[['BAI [File]']] = 'scRNA-p2277/scCount/Auto.bam.bai'
output[['PreprocessingLog [File]']] = 'scRNA-p2277/scCount/Auto_preprocessing.log'
output[['STARLog [File]']] = 'scRNA-p2277/scCount/Auto_STAR.log'
output[['CountMatrix [File]']] = 'scRNA-p2277/scCount/Auto-counts.txt'
output[['Stats [File]']] = 'scRNA-p2277/scCount/Auto-stats.txt'
output[['CellCyclePhase [File]']] = 'scRNA-p2277/scCount/Auto-CellCyclePhase.txt'
output[['CellDataset [File]']] = 'scRNA-p2277/scCount/Auto-dataset.txt'
EzAppSCCounts$new()$run(input=input, output=output, param=param)
## NEG
input[['Name']] = 'NEG'
input[['Read1']] = 'scRNA-p2277/HiSeq4000_20170822_RUN376_o3645/NEG_unmapped_part.bam'
input[['CellDataset']] = 'scRNA-p2277/HiSeq4000_20170822_RUN376_o3645/NEG-dataset.tsv'
output[['Name']] = 'NEG'
output[['BAM [File]']] = 'scRNA-p2277/scCount/NEG.bam'
output[['BAI [File]']] = 'scRNA-p2277/scCount/NEG.bam.bai'
output[['PreprocessingLog [File]']] = 'scRNA-p2277/scCount/NEG_preprocessing.log'
output[['STARLog [File]']] = 'scRNA-p2277/scCount/NEG_STAR.log'
output[['CountMatrix [File]']] = 'scRNA-p2277/scCount/NEG-counts.txt'
output[['Stats [File]']] = 'scRNA-p2277/scCount/NEG-stats.txt'
output[['CellCyclePhase [File]']] = 'scRNA-p2277/scCount/NEG-CellCyclePhase.txt'
output[['CellDataset [File]']] = 'scRNA-p2277/scCount/NEG-dataset.txt'
EzAppSCCounts$new()$run(input=input, output=output, param=param)

# p2497 SCCounts from uBam
setwd("/export/local/scratch/gtan/dev/p2497-SCCounts")
library(ezRun)
setEnvironments("star")
setEnvironments("flexbar")
setEnvironments("trimmomatic")
setEnvironments("python2")
setEnvironments("samtools")
param = list()
param[['cores']] = '8'
param[['ram']] = '50'
param[['scratch']] = '500'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['mapMethod']] = 'STAR'
param[['mapOptions']] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
param[['getChimericJunctions']] = 'false'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '4'
param[['trimRight']] = '0'
param[['minTailQuality']] = '10'
param[['minAvgQuality']] = '20'
param[['minReadLength']] = '20'
param[['featureLevel']] = 'gene'
param[['gtfFeatureType']] = 'exon'
param[['allowMultiOverlap']] = 'true'
param[['countPrimaryAlignmentsOnly']] = 'true'
param[['minFeatureOverlap']] = '10'
param[['minMapQuality']] = '1'
param[['keepMultiHits']] = 'true'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42'
output = list()
output[['Name']] = '20171222.A-SiCSeq_SCs_P5'
output[['Species']] = 'Mus musculus (house mouse)'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
output[['paired']] = 'false'
output[['refFeatureFile']] = 'genes.gtf'
output[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
output[['CellDataset [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5-dataset.tsv'
output[['CountMatrix [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5-counts.mtx'
output[['Stats [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5-stats.txt'
output[['CellCyclePhase [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5-CellCyclePhase.txt'
output[['BAM [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5.bam'
output[['BAI [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5.bam.bai'
output[['PreprocessingLog [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5_preprocessing.log'
output[['STARLog [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/20171222.A-SiCSeq_SCs_P5_STAR.log'
input = list()
input[['Name']] = '20171222.A-SiCSeq_SCs_P5'
#input[['Read1']] = 'p2497/HiSeq4000_20171222_RUN420_o3705_fixedRG/20171222.A-SiCSeq_SCs_P5_unmapped.bam'
input[['Read1']] = 'p2497/HiSeq4000_20171222_RUN420_o3705_fixedRG/20171222.A-SiCSeq_SCs_P5_subset.bam'
input[['Read Count']] = '312016990'
input[['Species']] = 'Mus musculus (house mouse)'
input[['CellDataset']] = 'p2497/HiSeq4000_20171222_RUN420_o3705_fixedRG/fq_dataset.tsv'
#debug(ezMethodSingleCellSTAR)
#debug(ezMethodSingleCellFeatureCounts)
EzAppSCCounts$new()$run(input=input, output=output, param=param)

# p2497 SCCounts from uBam, second example
library(ezRun)
setwd("/export/local/scratch/gtan/dev/p2497-SCCounts")
setEnvironments("star")
setEnvironments("flexbar")
setEnvironments("trimmomatic")
setEnvironments("python2")
setEnvironments("samtools")
param = list()
param[['cores']] = '8'
param[['ram']] = '50'
param[['scratch']] = '500'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
param[['paired']] = 'false'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['spikeInSet']] = ''
param[['mapMethod']] = 'STAR'
param[['mapOptions']] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000  --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
param[['getChimericJunctions']] = 'false'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '0'
param[['trimRight']] = '0'
param[['minTailQuality']] = '10'
param[['minAvgQuality']] = '10'
param[['minReadLength']] = '20'
param[['featureLevel']] = 'gene'
param[['gtfFeatureType']] = 'exon'
param[['allowMultiOverlap']] = 'true'
param[['countPrimaryAlignmentsOnly']] = 'true'
param[['minFeatureOverlap']] = '10'
param[['minMapQuality']] = '1'
param[['keepMultiHits']] = 'true'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2497/SCCountsApp_26718_2018-05-23--11-47-01'
output = list()
output[['Name']] = 'P60_2'
output[['Species']] = ''
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26'
output[['paired']] = 'false'
output[['featureLevel']] = 'gene'
output[['refFeatureFile']] = 'genes.gtf'
output[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
output[['CellDataset [File]']] = 'p2497/SCCountsApp_26718_2018-05-23--11-47-01/P60_2-dataset.tsv'
output[['CountMatrix [File]']] = 'p2497/SCCountsApp_26718_2018-05-23--11-47-01/P60_2-counts.mtx'
output[['Stats [File]']] = 'p2497/SCCountsApp_26718_2018-05-23--11-47-01/P60_2-stats.txt'
output[['CellCyclePhase [File]']] = 'p2497/SCCountsApp_26718_2018-05-23--11-47-01/P60_2-CellCyclePhase.txt'
output[['BAM [File]']] = 'p2497/SCCountsApp_26718_2018-05-23--11-47-01/P60_2.bam'
output[['BAI [File]']] = 'p2497/SCCountsApp_26718_2018-05-23--11-47-01/P60_2.bam.bai'
output[['PreprocessingLog [File]']] = 'p2497/SCCountsApp_26718_2018-05-23--11-47-01/P60_2_preprocessing.log'
output[['STARLog [File]']] = 'p2497/SCCountsApp_26718_2018-05-23--11-47-01/P60_2_STAR.log'
input = list()
input[['Name']] = 'P60_2'
#input[['Read1']] = 'p2497/HiSeq4000_20180518_RUN456_o4446/P60_2_unmapped.bam'
# samtools view -s 0.01 -b /srv/gstore/projects/p2497/HiSeq4000_20180518_RUN456_o4446/P60_2_unmapped.bam > P60_2_subset.bam
input[['Read1']] = 'p2497/HiSeq4000_20180518_RUN456_o4446/P60_2_subset.bam'
input[['Read Count']] = '318102190'
input[['Species']] = ''
input[['CellDataset']] = 'p2497/HiSeq4000_20180518_RUN456_o4446/scFastq_dataset.tsv'
#debug(ezMethodSingleCellFeatureCounts)
EzAppSCCounts$new()$run(input=input, output=output, param=param)


# p2214 SCCounts from fastqs
library(ezRun)
setEnvironments("star")
setEnvironments("flexbar")
setEnvironments("trimmomatic")
setEnvironments("python2")
setEnvironments("samtools")
setwd("/export/local/scratch/gtan/p2214-SCCounts")
param = list()
param[['cores']] = '8'
param[['ram']] = '50'
param[['scratch']] = '500'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['paired']] = 'true'
param[['strandMode']] = 'both'
param[['refFeatureFile']] = 'genes.gtf'
param[['mapMethod']] = 'STAR'
param[['mapOptions']] = '--outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif'
param[['getChimericJunctions']] = 'false'
param[['trimAdapter']] = 'true'
param[['trimLeft']] = '4'
param[['trimRight']] = '0'
param[['minTailQuality']] = '10'
param[['minAvgQuality']] = '20'
param[['minReadLength']] = '20'
param[['featureLevel']] = 'gene'
param[['gtfFeatureType']] = 'exon'
param[['allowMultiOverlap']] = 'true'
param[['countPrimaryAlignmentsOnly']] = 'true'
param[['minFeatureOverlap']] = '10'
param[['minMapQuality']] = '1'
param[['keepMultiHits']] = 'true'
param[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
param[['specialOptions']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2214/SCCountsApp_24762_2018-02-20--21-32-42'
output = list()
output[['Name']] = 'plate_2 '
output[['Species']] = 'Homoe sapiens'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
output[['paired']] = 'true'
output[['refFeatureFile']] = 'genes.gtf'
output[['transcriptTypes']] = 'protein_coding,rRNA,tRNA,Mt_rRNA,Mt_tRNA'
output[['CellDataset [File]']] = 'p2214/SCCountsApp_24762_2018-02-20--21-32-42/plate_2-dataset.tsv'
output[['CountMatrix [File]']] = 'p2214/SCCountsApp_24762_2018-02-20--21-32-42/plate_2-counts.txt'
output[['Stats [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/plate_2-stats.txt'
output[['CellCyclePhase [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/plate_2-CellCyclePhase.txt'
output[['BAM [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/plate_2.bam'
output[['BAI [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/plate_2.bam.bai'
output[['PreprocessingLog [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/plate_2_preprocessing.log'
output[['STARLog [File]']] = 'p2497/SCCountsApp_24762_2018-02-20--21-32-42/plate_2_STAR.log'
input = list()
input[['Name']] = 'plate_2'
#input[['Read1']] = 'p2214/HiSeq4000_20171222_RUN420_o3705_fixedRG/20171222.A-SiCSeq_SCs_P5_unmapped.bam'
#input[['Read1']] = 'p2497/HiSeq4000_20171222_RUN420_o3705_fixedRG/20171222.A-SiCSeq_SCs_P5_subset.bam'
#input[['Read Count']] = '312016990'
input[['Species']] = 'Homoe sapiens'
input[['CellDataset']] = 'p2214/HiSeq4000_20180309_RUN428_o4166_GeTesting3Samples/dataset.tsv'
#input[['CellDataset']] = 'p2214/HiSeq4000_20180309_RUN428_o4166/dataset.tsv'
debug(ezMethodSingleCellSTAR)
#debug(ezMethodSCCounts)
EzAppSCCounts$new()$run(input=input, output=output, param=param)

