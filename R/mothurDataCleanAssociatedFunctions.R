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
#rawDataSummaryTableTitle <- paste("Number of unique sequences =",numUniqReads , sep = " ")
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
    nameInDataset <- rownames(sushiInputDataset)[i]
    file1PathInDatset <- paste0(param$dataRoot,"/",sushiInputDataset$`Read1 [File]`[i])
    file2PathInDatset <- paste0(param$dataRoot,"/",sushiInputDataset$`Read2 [File]`[i])
#   initialTable <- cbind(nameInDataset,file1PathInDatset,file2PathInDatset)
#    write.table(initialTable, 'Illumina.files', row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE, sep = "\t")
    k=0
    for (file in c(file1PathInDatset,file2PathInDatset)){
      k=k+1
      cpCmd <- paste0("gunzip -c ", file, "  > ", nameInDataset,".R",k,".fastq")
    ezSystem(cpCmd)
  }
  }
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
##' 

chimeraSummaryPlot <- function(x){
    nameRawFile <- basename(x)
    plotLables <- gsub(".chimPlot.txt","",nameRawFile)
    chimeraFile <- read.delim(nameRawFile, header = F)
  BL <-  table(chimeraFile$V18)["?"]
  chim <- table(chimeraFile$V18)["Y"]
  noChim <- table(chimeraFile$V18)["N"]
  chimeraDF <- data.frame(rbind(chim,noChim,BL),stringsAsFactors = FALSE)
  colnames(chimeraDF) = "Freq"
  chimeraDF[is.na(chimeraDF$Freq),]= 0
  chimeraDF$Type  = as.factor(c("Chimeric","Not chimeric","Borderline"))
  chimeraDF$pct  = round(chimeraDF$Freq/sum(chimeraDF$Freq)*100,2)
  bp <- ggplot(chimeraDF, aes(x=Type,Freq, fill=Type)) 
  facetSampleBar <- bp  + geom_bar(stat = "identity",  position = 'dodge') 
  finalVersionChimeraPlot  <- facetSampleBar +   
    theme(axis.title.x=element_blank(), axis.text.x=element_blank()) 
  finalVersionChimeraPlot <- finalVersionChimeraPlot +
    geom_text(aes(y = Freq + 500, label = paste0(pct, '%')),
              position = position_dodge(width = .9),size = 3)
  return(finalVersionChimeraPlot)
}



###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Merges summary tables for Kable_extra output 
##' @description Merges summary tables for Kable_extra output 
##' @param  x list of files for which to create the combined table
##' @return Returns a ktable and the above header
createSummaryTableForKableExtra <- function(x) {
  rawSummaryTable <- list()
  multiTableHeader <- vector()
  for (file in x){
    nameRawFile <- basename(file)
    tableTitle <- unlist(strsplit(nameRawFile,"\\."))[1]
    rawFileDF <- ezRead.table(nameRawFile)
    rawSummaryTable[[file]] <- as.matrix(createSummaryTable(rawFileDF))
    multiTableHeader[tableTitle] <- "7"
  }
  ktables <- do.call(rbind, rawSummaryTable)
  return(list(mergedTable = ktables, aboveHeader = multiTableHeader))
}
 

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Creates step convergence table for Kable_extra output 
##' @description Merges summary tables for Kable_extra output 
##' @param  x list of files for which to create the combined table
##' @return Returns a ktable and the above header
createStepConvTableForKableExtra <- function(x) {
  rawSummaryTable <- list()
  multiTableHeader <- vector()
  for (file in x){
    nameRawFile <- basename(file)
    rawSummaryTable[[file]] <- convStepTable(nameRawFile)
  }
  finNumRows <- max(sapply(rawSummaryTable, nrow))
  for (file in x){
    nameRawFile <- basename(file)
    tableTitle <- unlist(strsplit(nameRawFile,"\\."))[1]
    multiTableHeader[tableTitle] <- finNumRows
  }
  fillEmptyDF <- function(inDF,totRows){
    ncols <- ncol(inDF)
    presRows <- nrow(inDF)
    df <- data.frame(matrix(ncol = ncols, nrow = totRows-presRows), stringsAsFactors = FALSE)
    names(df) <- names(inDF)
    boundDF <- rbind(inDF,df)
    return(as.matrix(boundDF))
  }
  rawSummaryTableFilled <- lapply(rawSummaryTable,function(x)fillEmptyDF(x,totRows = finNumRows))
  ktables <- do.call(rbind, rawSummaryTableFilled)
  return(list(mergedTable = ktables, aboveHeader = multiTableHeader))
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Creates saturation  convergence table for Kable_extra output 
##' @description Merges summary tables for Kable_extra output 
##' @param  x list of mothur shared files for which to create the combined table
##' @return Returns a ktable and the above header
createSaturationTableForKableExtra <- function(x) {
  rawSummaryTable <- list()
  multiTableHeader <- vector()
  for (file in x){
    nameRawFile <- basename(file)
    tableTitle <- unlist(strsplit(nameRawFile,"\\."))[1]
    rawSummaryTable[[file]] <- as.matrix(otuSaturationTable(nameRawFile))
    multiTableHeader[tableTitle] <- "10"
  }
  ktables <- do.call(rbind, rawSummaryTable)
  return(list(mergedTable = ktables, aboveHeader = multiTableHeader))
}
