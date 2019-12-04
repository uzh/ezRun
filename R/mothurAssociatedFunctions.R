###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Mothur Summary Table
##' @description Create summary table from mothur summary files.
##' @param  summary, mothur summary files.
##' @return Returns a data.frame.

createSummaryTable <- function(summary){
part2 <- vector()
part1 <- apply(subset(summary,select=start:polymer),2,function(x)quantile(x, probs = c(0, 0.025,0.25, 0.5, 0.75, 0.975,1)))
k=1
part2[k] = 1
for (i in c(2.5,25,50,75,97.5,100)) {
  k=k+1
  part2[k] <- nrow(summary)*i/100
}
part2 <- data.frame(numSeqs=part2)
rawDataSummaryTable <- round(cbind(part1,data.frame(part2)), digits = 0)
rownames(rawDataSummaryTable) <- c("Mininmun","2.5%-tile","25%-tile","Median","75%-tile","97.5%-tile","Maximum")
return(rawDataSummaryTable)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Phyloseq OTU object
##' @description Create Phyloseq OTU object mothur OTU files.
##' @param  otuFileName, mothur shared clustered OTU  files.
##' @return Returns a Phyloseq OTU object.

phyloSeqOTU <- function(otuFileName){
otuFile <- read.table(otuFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
rownames(otuFile) <- otuFile$Group
colToDrop <- c("label","Group","numOtus")
otuFile1 <- as.matrix(otuFile1[,!names(otuFile1)%in%colToDrop])
otuObject <- otu_table(otuFile, taxa_are_rows = FALSE)
return(otuObject)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Phyloseq Taxa object
##' @description Create Phyloseq taxa object from mothur taxonomy files.
##' @param  taxaFileName, mothur taxonomy file.
##' @return Returns a Phyloseq Taxa object.

phyloSeqTaxa <- function(taxaFileName,technology){
taxaFile <- read.table(taxaFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
tempList <- lapply(taxaFile$Taxonomy,function(y) unlist(strsplit(y,";")))
taxaMatrix <- as.matrix(ldply(tempList))
rownames(taxaMatrix) <- paste(taxaFile$OTU, technology, sep = "_")
colnames(taxaMatrix) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")[1:(ncol(taxaMatrix))]
taxaObject <- tax_table(taxaMatrix)
return(taxaObject)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Phyloseq Sample object
##' @description Create Phyloseq sample object from dataset.
##' @param  taxaFileName, mothur taxonomy file.
##' @return Returns a Phyloseq Taxa object.
phyloSeqSample <- function(sampleFileName){
sampleFile <- ezRead.table(sampleFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
colToKeep <- grep("Factor",colnames(sampleFile))
colnames(sampleFile) <- sub('\\s\\[Factor\\]',"",colnames(sampleFile))
sampleObject <- sample_data(sampleFile)
return(sampleObject)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Creates Mothur input files from  Sushi dataset.
##' @description Converts Sushi dataset into Mothur input.
##' @param  taxaFileName, mothur taxonomy file.
##' @return Returns the .groups and .fasta files
datasetToMothur <- function(sushiInputDataset, param){
for (i in (1:nrow(sushiInputDataset))) {
filePathInDatset <- paste0(param$dataRoot,"/",sushiInputDataset$`Read1 [File]`[i])
techID <- sushiInputDataset$`Technology [Factor]`[i]
groupID <- rownames(sushiInputDataset)[i]
fastqFile <- readFastq(filePathInDatset)
x=data.frame(fastqFile@id)
readID <- data.frame(apply(x,1,function(y) unlist(strsplit(y," "))[[1]]))
groupFile <- data.frame(apply(readID,1,function(y) gsub(":","_",y)))
groupFile$group <- groupID
if (techID == "Illumina"){
write.table(groupFile, 'Illumina.groups', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
writeFasta(fastqFile,'Illumina.fasta', mode = 'a')
}else{
  write.table(groupFile, 'PacBio.groups', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE)
  writeFasta(fastqFile,'PacBio.fasta', mode = 'a')
}
}
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Error rate
##' @description Summarizes error rate from error count file
##' @param  errorCountFile, mothur taxonomy file.
##' @return Returns a number.
  errorRateSummaryPlot <- function(errorCountFileName){
  errorTable <- read.table(errorCountFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  
  errorTable$wrongBases = errorTable$Sequences*errorTable$Errors
  errorTable$totBases = errorTable$Sequences*1450
  errorTable$Errors <- as.factor(errorTable$Errors)
  overallErrorRate <- round(sum(errorTable$wrongBases)/sum(errorTable$totBases)*100,2)
  pct <- cumsum(round(errorTable$Sequences/sum(errorTable$Sequences)*100,2))
  lbls <- paste(errorTable$Errors, " (",pct, "%)", sep = "") # add percents to labels 
  col=rainbow(length(lbls))
  titleText <- "Error distribution in the sequences"
  subtitleText <- paste0("Overall error rate = ", overallErrorRate, "%")
  bp <- ggplot(errorTable, aes(x="", y=Sequences, fill=Errors)) + geom_bar(width = 1, stat = "identity") + 
    scale_fill_manual(breaks=c(0,5,10,50), labels=lbls[c(1,4,9,49)], values=col)
  pieVersion <- bp + coord_polar("y", start=0)
  finalVersion <- pieVersion +  labs(title=titleText, subtitle=subtitleText) + 
    theme(plot.title=element_text(size=15, face="bold",hjust=0.5),  
          plot.subtitle=element_text(size=10, face="bold",hjust=0.5))
  return(finalVersion)
}

  ###################################################################
  # Functional Genomics Center Zurich
  # This code is distributed under the terms of the GNU General
  # Public License Version 3, June 2007.
  # The terms are available here: http://www.gnu.org/licenses/gpl.html
  # www.fgcz.ch
  
  
  ##' @title Clustering steps
  ##' @description Converegence iteration 
  ##' @param  convStepFile, mothur steps file.
  ##' @return Returns a table
  convStepTable <- function(convStepFile){
    stepTable <- ezRead.table(convStepFile)
    colToKeep <- c("iter","num_otus","sensitivity","specificity","fdr","accuracy")
    stepTable <- stepTable[,colToKeep]
    return(stepTable)
  }
  

 
  ###################################################################
  # Functional Genomics Center Zurich
  # This code is distributed under the terms of the GNU General
  # Public License Version 3, June 2007.
  # The terms are available here: http://www.gnu.org/licenses/gpl.html
  # www.fgcz.ch
  
  
  ##' @title OTUs saturation table
  ##' @description HOw many OTUs do we really have?
  ##' @param  sharedFile, mothur shared abundance  file.
  ##' @return Returns a table
  otuSaturationTable <- function(sharedFile){
    sharedAbund <- read.table(sharedFile, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
    sharedAbund <- t(sharedAbund)
    totOtus <- sharedAbund[rownames(sharedAbund) == "numOtus",]
    rowToKeep <- grepl("^Otu.*$",rownames(sharedAbund))
    sharedAbundDF <- data.frame(data.matrix(data.frame(sharedAbund[rowToKeep,], stringsAsFactors = FALSE)))
    colnames(sharedAbundDF) <- sharedAbund[rownames(sharedAbund) == "Group",]
    cumSumTransform <- data.frame(apply(sharedAbundDF,2,cumsum))
    tempList <- apply(cumSumTransform,2,function(y) 
      ldply(lapply(seq(from = 10, to = 100, by = 10),function(x) y[x]/y[nrow(cumSumTransform)]*100)))
    finalSaturationTable <- data.frame(matrix(unlist(tempList), nrow=10, byrow=F),stringsAsFactors=FALSE)
    colnames(finalSaturationTable) <- names(tempList)
    rownames(finalSaturationTable) <- seq(from = 10, to = 100, by = 10)
    return(finalSaturationTable)
  }
  
  ###################################################################
  # Functional Genomics Center Zurich
  # This code is distributed under the terms of the GNU General
  # Public License Version 3, June 2007.
  # The terms are available here: http://www.gnu.org/licenses/gpl.html
  # www.fgcz.ch
  
  
  ##' @title Mothur fasta summary
  ##' @description Summarizes read count across samples for a fasta file
  ##' @param  sharedFile, mothur fasta  file
  ##' @return Returns a data frame

  countAndAssignSeqsFromFasta <- function(fastaFile,filterStep,groupFile){
    listOfReads <- readDNAStringSet(fastaFile)
    groupDesc <- ezRead.table(groupFile, sep = " ", header = F)
    actualReadNames <- sapply(names(listOfReads), function(x) 
      unlist(strsplit(x," "))[1])
    if (fastaFile == "Mothur.fasta"){
      finaldDF <- cbind(data.frame(table(groupDesc)), fStep = filterStep)
    } else {
   rownames(groupDesc) <- gsub(":","_",rownames(groupDesc))
   names(groupDesc) <- "sample"
   relIndex <- which(rownames(groupDesc)%in%actualReadNames)
   finaldDF <- cbind(data.frame(table(groupDesc[relIndex,])), fStep = filterStep)
    }
    names(finaldDF)[1] <- "sample"
   return(finaldDF)
  }
    
  
  ###################################################################
  # Functional Genomics Center Zurich
  # This code is distributed under the terms of the GNU General
  # Public License Version 3, June 2007.
  # The terms are available here: http://www.gnu.org/licenses/gpl.html
  # www.fgcz.ch
  
  
  ##' @title Chimera identification 
  ##' @description Summarizes chimera rates from chimera file for sample
  ##' @param  chimerFile mothur chimera file.
  ##' @return Returns a data frame.
  
  chimeraSummaryTable <- function(chimFile,groupFile){
    groupDesc <- ezRead.table(groupFile, sep = " ", header = F)
    rownames(groupDesc) <- gsub(":","_",rownames(groupDesc))
    names(groupDesc) <- "sample"
    chimeraFile <- read.delim(chimFile, header = F, stringsAsFactors = F)
    chimeraFile <-  chimeraFile[!duplicated(chimeraFile$V2),]
    actualReadNames <- chimeraFile$V2
    listOfSamples <- unique(groupDesc$sample)
    getChimPerc <- function(sample){
    relIndex1 <- which(rownames(groupDesc)%in%actualReadNames)
    temp <- groupDesc[relIndex1,,drop=FALSE]
    usedReads <-  rownames(temp[temp$sample == sample,,drop=FALSE])
    chimeraFileSam <- chimeraFile[chimeraFile$V2%in%usedReads,]
    BL <-  table(chimeraFileSam$V18)["?"]
    chim <- table(chimeraFileSam$V18)["Y"]
    noChim <- table(chimeraFileSam$V18)["N"]
    chimeraDF <- data.frame(rbind(chim,noChim,BL),stringsAsFactors = FALSE)
    colnames(chimeraDF) = "Freq"
    chimeraDF[is.na(chimeraDF$Freq),]= 0
    chimeraDF$Type  = as.factor(c("Chimeric","Not chimeric","Borderline"))
    chimeraDF$pct  = round(chimeraDF$Freq/sum(chimeraDF$Freq)*100,2)
    chimeraDF$sample = sample
    return(chimeraDF)
    }
    listOfChimFiles <- lapply(listOfSamples,getChimPerc)
    finalDF <- do.call("rbind",listOfChimFiles)
    return(finalDF)
  }