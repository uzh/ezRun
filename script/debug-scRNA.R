setwd("/home/gtan/analysis/gtan/scRNA/smartSeq2")

refFn <- "/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Sequence/WholeGenomeFasta/genome.fa"

fastqFn <- "test_R1.fastq"
bamFn <- "test_R1"
fastqI1Fn <- "test_I1.fastq"
fastqI2Fn <- "test_I2.fastq"


fastqFn <- "20170822.B-Pool1_EBO6_11_A2B5_R1.fastq.gz"
bamFn <- "20170822.B-Pool1_EBO6_11_A2B5_R1"
fastqI1Fn <- "20170822.B-Pool1_EBO6_11_A2B5_I1.fastq.gz"
fastqI2Fn <- "20170822.B-Pool1_EBO6_11_A2B5_I2.fastq.gz"

system.time(fastq2bam(fastqFn, refFn, bamFn, fastqI1Fn=fastqI1Fn, fastqI2Fn=fastqI2Fn))
#user    system   elapsed 
#7641.624   484.496 13885.729


## Read group approach
setwd("/scratch/gtan/smartSeq2/splitFastq")
fastqFns <- list.files(path=".", pattern="\\.fastq\\.gz$", full.names = TRUE)
bamFn <- "merged_RG.bam"

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

countResult = Rsubread::featureCounts("merge_alignments.bam", annot.inbuilt=NULL,
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

# subread-1.5.3-Linux-x86_64/bin/featureCounts -a /srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_89-2017-05-31/Genes/genes.gtf -o featureCounts.txt --byReadGroup merge_alignments.bam