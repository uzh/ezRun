
###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Error rate
##' @description Summarizes error rate from error count file
##' @param  errorCountFile, mothur taxonomy file.
##' @return Returns a pie chart plot.
  errorRateSummaryPlot <- function(errorCountFileName,avReadLen){
  errorTable <- read.table(errorCountFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  errorTable$wrongBases = errorTable$Sequences*errorTable$Errors
  errorTable$totBases = errorTable$Sequences*avReadLen
  errorTable$Errors <- as.factor(errorTable$Errors)
  overallErrorRate <- round(sum(errorTable$wrongBases)/sum(errorTable$totBases)*100,2)
  pct <- cumsum(round(errorTable$Sequences/sum(errorTable$Sequences)*100,2))
  lbls <- paste(errorTable$Errors, " (",pct, "%)", sep = "") # add percents to labels 
  col=rainbow(length(lbls))
  titleText <- "Error distribution in the sequences"
  subtitleText <- paste0("Overall error rate = ", overallErrorRate, "%")
  bp <- ggplot(errorTable, aes(x="", y=Sequences, fill=Errors)) + geom_bar(width = 1, stat = "identity") + 
    scale_fill_manual(breaks=c(0,5,10,50), labels=lbls[c(1,6,11,51)], values=col)
  pieVersion <- bp + coord_polar("y", start=0)
  finalVersion <- pieVersion +  labs(title=titleText, subtitle=subtitleText, y="") + 
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
    stepTable <- read.table(convStepFile,stringsAsFactors = FALSE, sep = "\t", header = TRUE)
   # stepTable<- ezRead.table(convStepFile)
    colnames(stepTable) <- c(colnames(stepTable)[2:length(colnames(stepTable))],"dum")
    stepTable$iteration <- rownames(stepTable)
    colToKeep <- c("iteration","num_otus","sensitivity","specificity","fdr","accuracy")
    stepTable <- stepTable[,colToKeep]
    stepTable$iteration <- as.integer(stepTable$iteration)
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
      ldply(lapply(seq(from = 20, to = 200, by = 20),function(x) y[x]/y[nrow(cumSumTransform)]*100)))
    finalSaturationTable <- data.frame(matrix(unlist(tempList), nrow=10, byrow=F),stringsAsFactors=FALSE)
    finalSaturationTable <- data.frame(cbind(seq(from = 20, to = 200, by = 20),finalSaturationTable))
    colnames(finalSaturationTable) <- c("numOTUs","PercOTUsCovered")
    finalSaturationTableDF <- data.frame(PercOTUsCovered = finalSaturationTable$PercOTUsCovered, stringsAsFactors = F)
    rownames(finalSaturationTableDF) <- finalSaturationTable$numOTUs
    return(finalSaturationTableDF)
  }