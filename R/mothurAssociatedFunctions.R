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
errorRateSummary <- function(errorCountFileName){
  errorTable <- read.table(errorCountFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  errorTable$wrongBases = errorTable$Sequences*errorTable$Errors
  errorTable$totBases = errorTable$Sequences*1450
  errorRate <- sum(errorTable$wrongBases)/sum(errorTable$totBases)*100
  return(errorRate)
}
