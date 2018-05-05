###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodVirDetect = function(input=NA, output=NA, param=NA,
                             htmlFile="00index.html"){
  require(Rsamtools)
  sampleName = input$getNames() ##first parameter pass to rmarkdown::render
#  stopifnot((param$paired))
#  start_path = getwd()
  setwdNew(sampleName)
  
  ## trim reads
  trimmedInput = ezMethodTrim(input = input, param = param)

  ## align trimmed reads to human genome, get read pairs, for which both reads were unmapped
  defOpt = paste("-p", ezThreads())
  readGroupOpt = paste0("--rg-id ", sampleName," --rg SM:", sampleName,
                        " --rg LB:RGLB_", sampleName," --rg PL:illumina",
                        " --rg PU:RGPU_", sampleName)
  cmd = paste("bowtie2", param$cmdOptionsHost, defOpt, readGroupOpt,
              "-x", "/srv/GT/reference/Homo_sapiens/Ensembl/GRCh38.p10/Sequence/BOWTIE2Index/genome")
  if(param$paired){
    cmd = paste(cmd, "-1", trimmedInput$getColumn("Read1"),
              "-2", trimmedInput$getColumn("Read2"))
  } else {
    cmd = paste(cmd, "-U", trimmedInput$getColumn("Read1"))
  }
  cmd = paste(cmd, "2>", "bowtie2.log", "|", "samtools", "view -S -b -", 
              " > human.bam")
  ezSystem(cmd)
  if(param$paired){
    cmd = "samtools view -b -f 12 -F 256 human.bam > human.both_unmapped.bam"
  } else {
  cmd = "samtools view -b -f 4 -F 256 human.bam > human.both_unmapped.bam"
  }
  ezSystem(cmd)
  sortBam("human.both_unmapped.bam", "human.both_unmapped.sorted", byQname=TRUE,
          maxMemory=10240)
  #cmd = "samtools sort -n -m 10G human.both_unmapped.bam human.both_unmapped.sorted"
  #ezSystem(cmd)
  if(param$paired){
    cmd = "bedtools bamtofastq -i human.both_unmapped.sorted.bam -fq tr_human_removed_R1.fastq -fq2 tr_human_removed_R2.fastq"
  } else {
    cmd = "bedtools bamtofastq -i human.both_unmapped.sorted.bam -fq tr_human_removed_R1.fastq"
  }
  ezSystem(cmd) 
  
  ## align filtered reads to the selected host genome, get read pairs, in which both reads were unmapped
  paramHost <- param
  paramHost$refBuild <- param$hostBuild
  paramHost$ezRef <- EzRef(paramHost)
  host = getBowtie2Reference(paramHost)
  cmd = paste("bowtie2", param$cmdOptionsHost, defOpt, readGroupOpt,"-x", host)
  if(param$paired){
    cmd = paste(cmd, "-1 tr_human_removed_R1.fastq",
              "-2 tr_human_removed_R2.fastq",
              "2>>", "bowtie2.log", "|", "samtools", "view -S -b -", 
              " > host.bam")
  } else {
    cmd = paste(cmd, "-U tr_human_removed_R1.fastq",
                "2>>", "bowtie2.log", "|", "samtools", "view -S -b -", 
                " > host.bam")
  }
  ezSystem(cmd)
  if(param$paired){
    cmd = "samtools view -b -f 12 -F 256 host.bam > host.both_unmapped.bam"
  } else {
    cmd = "samtools view -b -f 4 -F 256 host.bam > host.both_unmapped.bam"
  }
  ezSystem(cmd)
  sortBam("host.both_unmapped.bam", "host.both_unmapped.sorted", byQname=TRUE,
          maxMemory=10240)
  #cmd = "samtools sort -n -m 10G host.both_unmapped.bam host.both_unmapped.sorted"
  #ezSystem(cmd)
  if(param$paired){
    cmd = "bedtools bamtofastq -i host.both_unmapped.sorted.bam -fq tr_host_removed_R1.fastq -fq2 tr_host_removed_R2.fastq"
  } else {
    cmd = "bedtools bamtofastq -i host.both_unmapped.sorted.bam -fq tr_host_removed_R1.fastq"
  }
  ezSystem(cmd)
  ## align filtered reads to the viral reference database, get sorted bam file and idex, output idxstats into a text file
  paramVirom <- param
  paramVirom$refBuild <- param$virBuild
  paramVirom$ezRef <- EzRef(paramVirom)
  vir = getBowtie2Reference(paramVirom)
  cmd = paste("bowtie2", param$cmdOptions, defOpt, readGroupOpt,
              "-x", vir)
  if(param$paired){
    cmd = paste(cmd, "-1 tr_host_removed_R1.fastq", 
              "-2 tr_host_removed_R2.fastq",
              "2>>", "bowtie2.log", "|", "samtools", "view -S -b -", 
              " > virome.bam")
  } else {
    cmd = paste(cmd, "-U tr_host_removed_R1.fastq",
                "2>>", "bowtie2.log", "|", "samtools", "view -S -b -", 
                " > virome.bam")
 }
  ezSystem(cmd)
  ezSortIndexBam("virome.bam", "virome.sorted.bam", ram=param$ram, 
                 removeBam=TRUE, cores=ezThreads())
  cmd = "samtools idxstats virome.sorted.bam > virome.idxstats.txt"
  ezSystem(cmd)
  
  bamFile <- "virome.sorted.bam"
  
  ## collect summary statistics and save in a summary table, collect per base coverage of each mapped viral genomes and save in individual csv files
  idx<-read.table("virome.idxstats.txt", header=FALSE, stringsAsFactors=FALSE)
  sub<-idx[idx$V3>0, ]
  csvFile = sub(".fa$", ".csv", paramVirom$ezRef["refFastaFile"])
  names<-read.csv(csvFile, quote="", stringsAsFactors=FALSE, header=FALSE)
  sub<-merge(sub, names, by="V1")
  if (nrow(sub)!=0){
  	for(i in 1:nrow(sub)) {
        	chr=sub[i, 1]
        	len=sub[i, 2]
        	common_name=sub[i, 6]
        	temp.df<-data.frame(c1=c(chr), c2=c("0"), c3=c(len), c4=c(common_name))
        	bed.file<-paste0(chr, ".bed")
        	csv.file<-paste0(chr, ".csv")
        	write.table(temp.df, file=bed.file, quote=FALSE, col.names=FALSE, 
        	            row.names=FALSE, sep="\t")
		system(paste0("samtools view -b ", bamFile, " ",
                              chr, " > ", chr, ".bam"))
                system(paste0("samtools index ", chr, ".bam"))
        	system(paste0("bedtools coverage -a ", bed.file, 
        	              " -b ", chr, ".bam", " -d > ", csv.file))
        	cov<-read.table(csv.file, header=FALSE, sep="\t", quote="",
        	                stringsAsFactors=FALSE)
        	sub[i,8]<-sum(cov$V6!=0)
        	sub[i,9]<-sum(cov$V6!=0)/len*100
        	sub[i,10]<-sum(cov$V6)/len
  	}
	sub<-sub[order(sub$V10, decreasing=TRUE), ]
  	out<-sub[, c(1,6,7,2,3,8,9,10)]
  	colnames(out)<-c("ID", "CommonName", "Family", "Len", "mappedReads", 
  	                 "mappedBases", "genomeCov_pect", "aveDepth")
  	write.table(out, file="summary_table.tsv", col.names=TRUE, row.names=FALSE, 
  	            quote=FALSE, sep="\t")
  }
  ##delete intermediate result files, this folder will be copied back to gstore
  ezSystem("rm -f *.bed")
  ezSystem("rm -f virome.bam")
  ezSystem("rm -f *host*")
  ezSystem("rm -f *human*")
  ezSystem("rm -f *.gz")
  ezSystem("rm -f *.fastq")
  
  ##html file  
  #setwd(start_path)
  htmlFile = output$getColumn("OutReport")
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "VirDetect.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  params = list(sample=sampleName,
                minReadCount=param$minReadCount)
  rmarkdown::render(input="VirDetect.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
  return("Success")
}

##' @template app-template
##' @templateVar method ezMethodSpades()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppVirDetect <-
  setRefClass("EzAppVirDetect",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodVirDetect
                  name <<- "EzAppVirDetect"
                  appDefaults <<- rbind(minReadCount = ezFrame(Type="integer", DefaultValue="9", Description="use for reporting only viral genomes with mapped reads higher than"))
                }
              )
)
