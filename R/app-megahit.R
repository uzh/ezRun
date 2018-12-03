###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMegahit = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  ### de novo metagenome assemby with Megahit, gene prediction with prodigal and annotation with diamond
  
  library(plyr)
  library(dplyr)
  
  sampleName = input$getNames()
  file1PathInDatset <- input$getFullPaths("Read1")
  fastqName1 <- paste0(sampleName,".R1.fastq")
  fastqName1Trimmed  <- paste0(sampleName,".R1.trimmed.fastq") 
  cpCmd1 <- paste0("gunzip -c ", file1PathInDatset, "  > ", fastqName1)
  ezSystem(cpCmd1)  
  if(param$paired){
    file2PathInDatset <- input$getFullPaths("Read2")
    fastqName2 <- paste0(sampleName,".R2.fastq")
    fastqName2Trimmed  <- paste0(sampleName,".R2.trimmed.fastq")
    cpCmd2 <- paste0("gunzip -c ", file2PathInDatset, "  > ", fastqName2)
    ezSystem(cpCmd2)
    inputString <- paste("\"-1", fastqName1Trimmed, "-2", fastqName2Trimmed,"\"")
    pairedString="YES"
  } else {
    inputString <- paste("-r",fastqName1Trimmed)
    pairedString="NO"
  }
  megahitTemplScript <- MEGAHIT_TEMPLATE_SCRIPT
  megahitToBeExec <- paste0("megahit.",sampleName,".sh")
  ##update template
  updateTemplateScriptCmd <- paste0("sed -e s/\"SAMPLE_NAME\"/", sampleName, "/g",
                                    " -e s/\"INPUT_FILE_STRING\"/", inputString, "/g ",
                                    " -e s/\"KMER_LIST\"/", param$megahitKmerList, "/g ",
                                    " -e s/\"ARE_READ_PAIRED\"/", pairedString, "/g ",
                                    megahitTemplScript, " >",
                                    megahitToBeExec)
  ezSystem(updateTemplateScriptCmd)
  ## run script
  megahitCmd <- paste("bash",megahitToBeExec)
  ezSystem(megahitCmd)
  ### megahit complete
  
  ## format output  OTUcount and OTUtax files for phyloseq
  ##OTU count
  covFile <- ezRead.table("mockCommShotgun.cov.txt")
  DF1 <- t(data.frame(t(covFile)["Median_fold",], stringsAsFactors = F))
  OTUsCountTable <- data.frame(Group="test", DF1, stringsAsFactors = F)


  ##OTU tax
  OTUsTaxAbun <- data.frame(t(subset(OTUsCountTable,select=-c(Group))), 
                             stringsAsFactors = F)
  colnames(OTUsTaxAbun) <- "Size"
  
  krakLab <- read.delim("kraken.labels", header = F, stringsAsFactors = F)
  formatKrakenRows <- function(x){
  t <- paste(unlist(strsplit(x,";"))[3:9], 
             collapse = ";")
  t1 <- unlist(strsplit(t," "))[1]
               return(t1)
  }
  taxonomyCol <- data.frame(ldply(lapply(krakLab$V2,formatKrakenRows))$V1,
                            stringsAsFactors = F)
  rownames(taxonomyCol) <- krakLab$V1
  colnames(taxonomyCol) <- "Taxonomy"
  
     formatTaxFields <- function(x){
       taxOut <- vector()
       taxOut[1] <- "Bacteria"
      taxClassif <- unlist(strsplit(x,";"))
      for (k in 2:7){
        isDomThere <- taxClassif[k]
        if (is.na(isDomThere)){
          k1=k-1
          alreadyUnclassified <- grep(taxClassif[k1],"unclassified")
          if (length(alreadyUnclassified) > 0 ) {
            taxOut[k] <- taxOut[k1]
          } else{
            taxOut[k] <- paste(taxClassif[k1],"unclassified",sep = "_")
          }
        }else {
          taxOut[k]  <- taxClassif[k]
        }
      }
      return(paste(taxOut, collapse = ";"))
     }
      
  taxonomyCol$Taxonomy <-  sapply(taxonomyCol$Taxonomy,formatTaxFields, USE.NAMES = F)
  finalOTUTaxTable <- merge(OTUsTaxAbun,taxonomyCol,by="row.names")
  colnames(finalOTUTaxTable)[1] <- "OTU"
  colsToKeep <- finalOTUTaxTable$OTU
  finalOTUCountTable <- OTUsCountTable[,c("Group",colsToKeep)]
  ## write final tables ready for phyloseq
  write.table(finalOTUTaxTable, "tempTaxaFile.txt",sep = "\t",
              col.names = T, row.names = F, quote = F)
  write.table(finalOTUCountTable, "tempCountFile.txt",sep = "\t",
              col.names = T, row.names = F, quote = F)
  ## rename output files
  #1) 
  oldContigFile <- "megahitResults/final.contigs.fa"
  newContigFile <- basename(output$getColumn("contigFile"))
  ezSystem(paste("mv",oldContigFile,newContigFile))
  #1)binned
  oldContigFile <- "binnedContigs.fasta"
  newContigFile <- basename(output$getColumn("binnedContigsFile"))
  ezSystem(paste("mv",oldContigFile,newContigFile))
  #2) 
  oldProdigalFile <- "prodigalAnnotation.gff"
  newProdigalFile <- basename(output$getColumn("prodigalPredictionFile"))
  ezSystem(paste("mv",oldProdigalFile,newProdigalFile))
  #3) 
  oldIPSAnnFile <- "interProScanOut.gff3"
  newIPSAnnFile <- basename(output$getColumn("interproscanFile"))
  ezSystem(paste("mv",oldIPSAnnFile,newIPSAnnFile))
  #4) 
  oldTaxaFile <- "tempTaxaFile.txt"
  newTaxaFile <- basename(output$getColumn("OTUsToTaxonomyFile"))
  ezSystem(paste("mv",oldTaxaFile,newTaxaFile))
  #5) 
  oldCountFile <- "tempCountFile.txt"
  newCountFile <- basename(output$getColumn("OTUsCountTable"))
  ezSystem(paste("mv",oldCountFile,newCountFile))
}

##' @template app-template
##' @templateVar method ezMethodMegahit()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMegahit <-
  setRefClass("EzAppMegahit",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMegahit
                  name <<- "EzAppMegahit"
                  appDefaults <<- rbind(megahitKmerList = ezFrame(Type="character",
                                                                DefaultValue="69,79,89",
                                                                Description="Comma-separated list of k-mer for the assembly."),
                                        ### TO DO: allow diamond option
                                        diamondEvalue = ezFrame(Type="numeric",
                                                                  DefaultValue="0.05",
                                                                  Description="Blast e-value cut-off."),
                                        diamondMaxSeqs = ezFrame(Type="integer",
                                                                DefaultValue="30",
                                                                Description="Blast maximum number of sequences to report.")
                                        
                                        
                  )
                }
              )
  )

