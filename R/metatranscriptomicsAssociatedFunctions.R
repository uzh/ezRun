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

convertDiamondAnnotationToAbund <- function(annFile,feature){
  annDF <- read.delim(annFile, header = F, stringsAsFactors = F)
  foundIds <- annDF$V2
  if (feature == "function"){
  funcCatTemp <- sapply(foundIds, function(x) unlist(strsplit(x,"\\.1_"))[2], USE.NAMES = F)
  funcCatTemp <- sapply(funcCatTemp, function(x) unlist(strsplit(x,"_\\["))[1], USE.NAMES = F)
  funcCatTemp <- gsub("MULTISPECIES:_","",funcCatTemp)
  funcCat <- gsub(",_partial","",funcCatTemp)
  funcAbundDF <- data.frame(table(funcCat))
  rownames(funcAbundDF) <- funcAbundDF$funcCat
  funcAbundDF <- subset(funcAbundDF,select=-c(funcCat))
  names(funcAbundDF) <- gsub(".RefSeq.annotated.txt","",basename(annFile))
  return(data.frame(funcAbundDF, stringsAsFactors = F))
  } else if (feature == "organism"){
  orgTemp <- sapply(foundIds, function(x) unlist(strsplit(x,"_\\["))[2], USE.NAMES = F)
  orgTemp2 <- gsub("\\]","",orgTemp)
  org <- gsub("\\[","",orgTemp2)
  orgAbundDF <- data.frame(table(org))
  rownames(orgAbundDF) <- orgAbundDF$org
  orgAbundDF <- subset(orgAbundDF,select=-c(org))
  names(orgAbundDF) <- gsub(".RefSeq.annotated.txt","",basename(annFile))
  return(data.frame(orgAbundDF, stringsAsFactors = F))
  }
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title listOfAbundMerge
##' @description Create summary table from list of individual abundance.
##' @param  a list of data frame, one per sample.
##' @return Returns a data.frame.

listOfAbundMerge <- function(listOfAbund,names){
mergedDF <- Reduce(function(x, y) {
  z <- merge(x, y, all=TRUE,by=0)
  rownames(z) <- z$Row.names
  z <- subset(z, select=-c(Row.names))
  return(z)},
  listOfAbund)
names(mergedDF) <- names
mergedDF[is.na(mergedDF)] = 0
mergedDF <- mergedDF[apply(mergedDF,1,sd)>0,]
return(mergedDF)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title diffAbundReport
##' @description Differential abundance analysis between groups.
##' @param  a count and a design matrix
##' @return Returns a list of tables and plots.
metagMetatrDeseq2_tableAndPlots <- function(countMatrix,designMatrix,group,sampleGroup,refGroup){
  design <- model.matrix(~0+designMatrix[[group]])
  specifContrast <- paste0("C=",sampleGroup,"-",refGroup)
  myContrast <- makeContrasts(specifContrast,levels=design)
  filteredMatrix <- filterCountsForDiffExpr(countMatrix,designMatrix,group1,group2)
  y <- DGEList(counts=filteredMatrix, group=designMatrix$group)
  y <- calcNormFactors(y)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, contrast = myContrast[,x])
  DF <- qlf@.Data[[17]]
  DF <- DF[order(DF$PValue),]
  DF <- round(DF,5)
  DF$ID <- rownames(DF)
  addTaxa <- DF
  addTaxaOut <- addTaxaOut[order(addTaxaOut$padj),]
  colsToKeep <- c("logFC","PValue","ID")
  addTaxa <- addTaxaOut[,colsToKeep]
  names(addTaxa)[2] <- padj
  ## select top 20 to report in the table

  tableToReport <- addTaxa[1:20,colsToReport]
  
  ## mark significance
  addTaxa$Significance <- "Significant"
  addTaxa[addTaxa$padj > 0.05,]$Significance <- "nonSignificant"
  addTaxa$Significance <- as.factor(addTaxa$Significance)
  addTaxa$Significance <- factor(addTaxa$Significance, levels = rev(levels(addTaxa$Significance)))

  ### volcano plot
  title <- "Volcano plot (p-value threshold  = 0.05)"
  volcanoPlot <- ggplot(addTaxa, aes(y=-log10(pvalue), x=log2FoldChange)) +
    geom_point(aes(color=Significance),size=3) 
  volcanoPlot <- volcanoPlot + labs(title=title) + 
    theme(plot.title=element_text(size=10,hjust=0.5))
  volcanoPlot <- volcanoPlot + geom_hline(yintercept=1.3, color="blue") + 
    labs(color=Significance)
  ### Add barplots stacking org on x=samples for the top N diff. expr. func

  }
  return(list(vPlot=volcanoPlot,tableToReport=tableToReport,
              fullTable=addTaxa))
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title filterCountsForDiffExpr
##' @description Prepares matrix for differential analysis with specific contrast.
##' @param  a count matrix, a design matrix, the groups in the contrast 
##' @return Returns data.frame
filterCountsForDiffExpr <- function(DF, desMat, group1, group2){
  samplesGroup1 <- desMat[desMat$group  == group1 ,]$sample
  samplesGroup2 <- desMat[desMat$group  == group2,]$sample
  countGroup1 <- DF[,colnames(DF)%in%samplesGroup1]
  countGroup2 <- DF[,colnames(DF)%in%samplesGroup2]
  isPresGroup1 <- apply(countGroup1,1,function(x) sum(x >5)>1)
  isPresGroup2 <- apply(countGroup2,1,function(x) sum(x >5)>1)
  subsetDF <- DF[isPresGroup1 |isPresGroup2,]
  return(subsetDF)
}
