

require(ezRun)
require(ShortRead)

setwdNew("/srv/GT/analysis/hubert/yeast_10k")


FLUX <<- "export FLUX_MEM=\"20G\"; /usr/local/ngseq/src/flux-simulator-1.2.1/bin/flux-simulator"
my.system(paste(FLUX, "--printParameters"))
my.system(paste(FLUX, "--help"))
parLines = c("GEN_DIR"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Chromosomes",
             "PAIRED_END"="TRUE",
             "READ_LENGTH"="36",
             "UNIQUE_IDS"="true",
             "READ_NUMBER"="10000",
             #"REF_FILE_NAME"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf",
             "REF_FILE_NAME"="/srv/GT/analysis/hubert/yeast_10k/clean.gtf")
# grep -v _codon genes_sorted.gtf > clean.gtf
my.write.table(parLines, file="wt.par", col.names = FALSE, append = FALSE)
my.system(paste(FLUX, "-x", "-p", "wt.par"))

parLines = c("GEN_DIR"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Chromosomes",
             "PAIRED_END"="TRUE",
             "READ_LENGTH"="36",
             "UNIQUE_IDS"="true",
             "READ_NUMBER"="10000",
             #"REF_FILE_NAME"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf",
             "REF_FILE_NAME"="/srv/GT/analysis/hubert/yeast_10k/clean.gtf",
             "PRO_FILE_NAME"="wt.pro",
             "FASTA"="true")
my.write.table(parLines, file="wt_1.par", col.names = FALSE, append = FALSE)
my.system(paste(FLUX, "-s -l", "-p", "wt_1.par"))



pro = my.read.table("wt.pro", header=FALSE, row.names=NULL)
sum(pro$V5)
pro[1:50, "V6"] = pro[1:50, "V6"] *2
pro[51:100, "V6"] = round(pro[51:100, "V6"] /2)
pro$V5 = pro$V6 / sum(pro$V6)
sum(pro$V6)
my.write.table(pro, file="mut.pro", col.names=FALSE, row.names=FALSE)


parLines = c("GEN_DIR"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Chromosomes",
             "PAIRED_END"="true",
             "READ_LENGTH"="36",
             "ERR_FILE"="35",
             "UNIQUE_IDS"="true",
             "READ_NUMBER"="20000",
             #"REF_FILE_NAME"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf",
             "REF_FILE_NAME"="/srv/GT/analysis/hubert/yeast_10k/clean.gtf",
             "PRO_FILE_NAME"="mut.pro",
             "FASTA"="true")
my.write.table(parLines, file="mut_1.par", col.names = FALSE, append = FALSE)
my.system(paste(FLUX, "--threads 8", "-s -l", "-p", "mut_1.par"))


parLines = c("GEN_DIR"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Chromosomes",
             "PAIRED_END"="true",
             "READ_LENGTH"="36",
             "ERR_FILE"="35",
             "UNIQUE_IDS"="true",
             "READ_NUMBER"="20000",
             #"REF_FILE_NAME"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf",
             "REF_FILE_NAME"="/srv/GT/analysis/hubert/yeast_10k/clean.gtf",
             "PRO_FILE_NAME"="mut.pro",
             "FASTA"="true")
my.write.table(parLines, file="mut_2.par", col.names = FALSE, append = FALSE)
my.system(paste(FLUX, "--threads 8", "-s -l", "-p", "mut_2.par"))



# parLines = c("GEN_DIR"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Chromosomes",
#              "PAIRED_END"="true",
#              "READ_LENGTH"="36",
#              "ERR_FILE"="35",
#              "UNIQUE_IDS"="true",
#              "READ_NUMBER"="20000",
#              #"REF_FILE_NAME"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf",
#              "REF_FILE_NAME"="/srv/GT/analysis/hubert/yeast_10k/clean.gtf",
#              "PRO_FILE_NAME"="wt.pro",
#              "FASTA"="true")
# my.write.table(parLines, file="wt_1.par", col.names = FALSE, append = FALSE)
# my.system(paste(FLUX, "--threads 8", "-s -l", "-p", "wt_1.par"))


parLines = c("GEN_DIR"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/Chromosomes",
             "PAIRED_END"="true",
             "READ_LENGTH"="36",
             "ERR_FILE"="35",
             "UNIQUE_IDS"="true",
             "READ_NUMBER"="20000",
             #"REF_FILE_NAME"="/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Annotation/Genes/genes.gtf",
             "REF_FILE_NAME"="/srv/GT/analysis/hubert/yeast_10k/clean.gtf",
             "PRO_FILE_NAME"="wt.pro",
             "FASTA"="true")
my.write.table(parLines, file="wt_2.par", col.names = FALSE, append = FALSE)
my.system(paste(FLUX, "--threads 8", "-s -l", "-p", "wt_2.par"))




## build the dataset
dir.create("yeast_10k")
require(ShortRead)
ds = data.frame(row.names=c("wt_1", "wt_2", "mut_1", "mut_2"))
ds$"Genotype [Factor]" = c("wt", "wt", "mut", "mut")
ds$"Read1 [File]" = ""
ds$"Read2 [File]" = ""
ds$"Read Count" = 0
sm = "mut_2"
for (sm in rownames(ds)){
  reads = readFastq(paste(sm, ".fastq", sep=""))
  r1 = reads[grepl("/1$", id(reads))]
  r2 = reads[grepl("/2$", id(reads))]
  use = width(r1) == 36 & width(r2) == 36
  r1@id = BStringSet(sub("/1", "", id(r1)))
  outFile = paste("yeast_10k/", sm, "_R1.fastq", sep="")
  writeFastq(r1[use], file = outFile, compress=FALSE)
  my.system(paste("pigz -p 4 --best", outFile))
  ds[sm, "Read1 [File]"] = paste("extdata", "/", outFile, ".gz", sep="")
  r2@id = BStringSet(sub("/2", "", id(r2)))
  outFile = paste("yeast_10k/", sm, "_R2.fastq", sep="")
  writeFastq(r2[use], file = outFile, compress=FALSE)
  my.system(paste("pigz -p 4 --best", outFile))
  ds[sm, "Read2 [File]"] = paste("extdata", "/", outFile, ".gz", sep="")
  ds[sm, "Read Count"] = sum(use)
}
ds$"Adapter1" = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ds$"Adapter2" = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
ds$strandMode = "sense"
ds$Species = "S. cerevisiae"
ds$"Enrichment Kit" = "poly-A"
my.write.table(ds, file="yeast_10k/dataset.tsv", head="Name")



ds = my.read.table("inst/extdata/yeast_10k/dataset.tsv")
ds$"Adapter1" = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ds$"Adapter2" = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
my.write.table(ds, file="inst/extdata/yeast_10k/dataset.tsv", head="Name")

