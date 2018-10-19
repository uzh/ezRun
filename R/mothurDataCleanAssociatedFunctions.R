###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Mothur Summary Table
##' @description Create summary table from mothur summary files.
##' @param  summary, ezTable from mothur summary files.
##' @return Returns a data.frame.

createSummaryTable <- function(summary){
numUniqReads <- nrow(summary)
rawDataSummaryTableTitle <- paste("Numer of unique sequences =",numUniqReads , sep = " ")
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
return(list(rawDataSummaryTable=rawDataSummaryTable, rawDataSummaryTableTitle=rawDataSummaryTableTitle))
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
mixedDatasetToMothur <- function(sushiInputDataset, param){
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
write.table(groupFile, 'Illumina.groups', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE, sep = "\t")
writeFasta(fastqFile,'Illumina.fasta', mode = 'a')
}else{
  write.table(groupFile, 'PacBio.groups', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE, sep = "\t")
  writeFasta(fastqFile,'PacBio.fasta', mode = 'a')
}
}
}

IlluminaDatasetToMothur <- function(sushiInputDataset, param){
  for (i in (1:nrow(sushiInputDataset))) {
    filePathInDatset <- paste0(param$dataRoot,"/",sushiInputDataset$`Read1 [File]`[i])
    groupID <- rownames(sushiInputDataset)[i]
    fastqFile <- readFastq(filePathInDatset)
    x=data.frame(fastqFile@id)
    readID <- data.frame(apply(x,1,function(y) unlist(strsplit(y," "))[[1]]))
    groupFile <- data.frame(apply(readID,1,function(y) gsub(":","_",y)))
    groupFile$group <- groupID
      write.table(groupFile, 'Illumina.groups', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE, sep = "\t")
      writeFasta(fastqFile,'Illumina.fasta', mode = 'a')
  }
}

prepareFilesLocallyForMothur <- function(sushiInputDataset, param){
  for (i in (1:nrow(sushiInputDataset))) {
    nameInDataset <- sushiInputDataset$Name
    file1PathInDatset <- paste0(param$dataRoot,"/",sushiInputDataset$`Read1 [File]`[i])
    file2PathInDatset <- paste0(param$dataRoot,"/",sushiInputDataset$`Read2 [File]`[i])
    initialTable <- cbind(nameInDataset,file1PathInDatset,file2PathInDatset)
    write.table(initialTable, 'Illumina.files', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE, sep = "\t")
  }
}
###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title OTUs saturation table conversion for VAMPs
##' @description writeOTUgzFileForVamps
##' @param  sharedFile,taxafile mothur shared abundance and taxa files.
##' @return Writes a gzipped file
writeOTUgzFileForVamps <- function(sharedFile, taxaFile){
  sharedAbund <- read.table(sharedFile, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  taxaAssign <- read.table(taxaFile, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  sharedAbund <- t(sharedAbund)
  rowToKeep <- grepl("^Otu.*$",rownames(sharedAbund))
  sharedAbundDF <- data.frame(data.matrix(data.frame(sharedAbund[rowToKeep,], stringsAsFactors = FALSE)))
  sharedAbundDF <- data.frame(cbind(taxaAssign$OTU,sharedAbundDF,taxaAssign$Taxonomy))
  colnames(sharedAbundDF) <- c("Cluster_ID",sharedAbund[rownames(sharedAbund) == "Group",],"Taxonomy")
  gz1 <- gzfile("OTUfileForVamps.txt.gz", "w")
  write.table(sharedAbundDF, gz1, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
  close(gz1)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Chimera identification 
##' @description Summarizes chimera rates from chimera file 
##' @param  chimerFile ezTable from  mothur chimera file.
##' @return Returns a pie chart plot.
chimeraSummaryPlot <- function(chimeraFile){
  BL <- sum(chimeraFile[chimeraFile$V18 == "?",]$V13)
  chim <- sum(chimeraFile[chimeraFile$V18 == "Y",]$V13)
  noChim <- sum(chimeraFile[chimeraFile$V18 == "N",]$V13)
  chimeraDF <- data.frame(rbind(chim,noChim,BL),stringsAsFactors = FALSE)
  colnames(chimeraDF) = "Freq"
  chimeraDF$Type  = as.factor(c("Chimeric","Not chimeric","Borderline"))
  pct <- round(chimeraDF$Freq/sum(chimeraDF$Freq)*100,2)
  titleText <- "Chimeric sequences in the sample"
  bp <- ggplot(chimeraDF, aes(x="", y=Freq, fill=Type)) + 
    geom_bar(position = position_stack(),width = 1, stat = "identity") + 
    geom_text(aes(label = pct), position = position_stack(vjust = 0.5),  size = 5)
      pieVersion <- bp + coord_polar("y", start=0)
  finalVersionChimeraPlot <- pieVersion +  labs(title=titleText, y="") + 
    theme(plot.title=element_text(size=15, face="bold",hjust=0.5))
  return(finalVersionChimeraPlot)
}
 