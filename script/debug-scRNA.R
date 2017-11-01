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

## Mapping pipeline
# module add Dev/jdk/8 Aligner/STAR/2.5.3a Tools/samtools/1.5 Aligner/BWA/0.7.15 Aligner/Bowtie/1.2.1.1 Aligner/Bowtie2/2.3.2 Aligner/TopHat/2.1.1 QC/Trimmomatic/0.36 QC/Flexbar/3.0.3 Tools/Picard/2.9.0 Dev/Python2/2.7.13
# java -jar $Picard_jar SamToFastq  I=20170822.B-Pool1_EBO6_11_A2B5_R1.bam  FASTQ=20170822.B-Pool1_EBO6_11_A2B5_R1.fastq
# cp /srv/GT/databases/contaminants/allIllumina-forTrimmomatic-20160202.fa adapters.fa
# java -jar $Trimmomatic_jar SE -threads 8 -phred33 20170822.B-Pool1_EBO6_11_A2B5_R1.fastq trimmed-R1.fastq ILLUMINACLIP:adapters.fa:1:20:7:5:true SLIDINGWINDOW:4:10 AVGQUAL:10 MINLEN:30 > trimmomatic.out 2> trimmomatic.err
# cat trimmomatic.err >> merged_RG_preprocessing.log
# flexbar --threads 8 -r trimmed-R1.fastq -u 20 --pre-trim-left 2 --pre-trim-right 0 --min-read-length 30 --target flexbar > flexbar.out 2> flexbar.err
# cat flexbar.out >> merged_RG_preprocessing.log
# mv flexbar.fastq merged_RG_trimmed.fastq
# STAR --genomeDir /srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31/Genes/genes_STARIndex --sjdbOverhang 150 --readFilesIn merged_RG_trimmed.fastq --twopassMode None --runThreadN 16 --genomeLoad LoadAndKeep --outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 100000 --alignMatesGapMax 100000 --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted > Aligned.out.bam

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
#java -jar $Picard_jar SamToFastq  I=20170822.B-Pool1_EBO6_11_A2B5_R1.bam  FASTQ=20170822.B-Pool1_EBO6_11_A2B5_R1.fastq

# Mapping with merged bam
## module add Dev/jdk/8 Aligner/STAR/2.5.3a Tools/samtools/1.5 Aligner/BWA/0.7.15 Aligner/Bowtie/1.2.1.1 Aligner/Bowtie2/2.3.2 Aligner/TopHat/2.1.1 QC/Trimmomatic/0.36 QC/Flexbar/3.0.3 Tools/Picard/2.9.0 Dev/Python2/2.7.13
## not used: samtools fastq -0 merged_RG.fastq.gz merged_RG.bam
## java -jar $Picard_jar SamToFastq  I=merged_RG.bam  FASTQ=merged_RG.fastq
## cp /srv/GT/databases/contaminants/allIllumina-forTrimmomatic-20160202.fa adapters.fa
## java -jar $Trimmomatic_jar SE -threads 8 -phred33 merged_RG.fastq trimmed-R1.fastq ILLUMINACLIP:adapters.fa:1:20:7:5:true SLIDINGWINDOW:4:10 AVGQUAL:10 MINLEN:30 > trimmomatic.out 2> trimmomatic.err
## cat trimmomatic.err >> merged_RG_preprocessing.log
## flexbar --threads 8 -r trimmed-R1.fastq -u 20 --pre-trim-left 2 --pre-trim-right 0 --min-read-length 30 --target flexbar > flexbar.out 2> flexbar.err
## cat flexbar.out >> merged_RG_preprocessing.log
## mv flexbar.fastq merged_RG_trimmed.fastq
## STAR --genomeDir /srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31/Genes/genes_STARIndex --sjdbOverhang 150 --readFilesIn merged_RG_trimmed.fastq --twopassMode None --runThreadN 16 --genomeLoad LoadAndKeep --outFilterType BySJout --outFilterMatchNmin 30 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.05 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 --alignIntronMax 100000 --alignMatesGapMax 100000 --outFilterMultimapNmax 50 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --outSAMstrandField intronMotif --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted > Aligned.out.bam
## not used: samtools sort -l 9 -m 612M -@ 8 Aligned.out.bam -o merged_STAR.bam
## not used: samtools index merged_STAR.bam
## not used: infer_experiment.py -r /srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31/Genes/genes.bed -i merged_STAR.bam -s 1000000

## java -jar $Picard_jar MergeBamAlignment ALIGNED=Aligned.out.bam UNMAPPED=merged_RG.bam O=merge_alignments.bam R=/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Sequence/WholeGenomeFasta/genome.fa
#samtools sort -@ 8 -o merge_alignments_sorted.bam merge_alignments.bam
countResult = Rsubread::featureCounts("20171011.A-C1_HT_24H.bam", annot.inbuilt=NULL,
                                      annot.ext="/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31/Genes/genes.gtf", isGTFAnnotationFile=TRUE,
                                      GTF.featureType="exon",
                                      GTF.attrType="gene_id",
                                      useMetaFeatures=TRUE,
                                      allowMultiOverlap=TRUE, isPairedEnd=FALSE, 
                                      requireBothEndsMapped=FALSE,
                                      checkFragLength=FALSE,minFragLength=50,maxFragLength=600,
                                      nthreads=4,
                                      strandSpecific=0,
                                      minMQS=1,
                                      readExtension5=0,readExtension3=0,read2pos=NULL,
                                      minOverlap=10,
                                      ignoreDup=NA,
                                      splitOnly=FALSE,
                                      countMultiMappingReads=TRUE,
                                      fraction=FALSE,
                                      primaryOnly=TRUE,
                                      countChimericFragments=TRUE,
                                      chrAliases=NULL,reportReads=NULL,
                                      byReadGroup=TRUE)


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

## Scater
library(ezRun)
param = list()
param[['cores']] = '1'
param[['ram']] = '8'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'SAMPLE'
param[['name']] = 'Scater_QC'
param[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['specialOptions']] = ''
param[['mail']] = ''
param[['dataRoot']] = '/scratch/gtan'
param[['resultDir']] = 'scRNA-p2277/Scater'
output = list()
output[['Name']] = 'Scater_QC'
output[['Species']] = 'Mus musculus'
output[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
output[['Static Report [Link]']] = 'p2277/ScaterApp_20774_2017-09-05--07-25-28/Scater_QC/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_scater_app/?data=p2277/ScaterApp_20774_2017-09-05--07-25-28/Scater_QC/counts-cbbohwdtniyq-EzResult.RData'
output[['Report [File]']] = 'p2277/ScaterApp_20774_2017-09-05--07-25-28/Scater_QC'
input = list()
input[['Name']] = 'NEG'
input[['CountDataset']] = 'p2277/SingleCellCountsApp_20773_2017-09-04--17-25-44/NEG-dataset.tsv'
input[['CountMatrix']] = 'p2277/SingleCellCountsApp_20773_2017-09-04--17-25-44/NEG-counts.txt'
input[['CountFolder']] = 'p2277/SingleCellCountsApp_20773_2017-09-04--17-25-44/NEG-Counts'
input[['CellCyclePhase']] = 'p2277/SingleCellCountsApp_20773_2017-09-04--17-25-44/NEG-CellCyclePhase.txt'
input[['Species']] = 'Mus musculus'
input[['refBuild']] = 'Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31'
input[['paired']] = 'false'
input[['refFeatureFile']] = 'genes.gtf'
EzAppScater$new()$run(input=input, output=output, param=param)
