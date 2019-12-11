###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Mothur Summary Table
##' @description Create summary table from mothur summary files.
##' @param  summary, mothur summary files.
##' @return Returns a list of data.frames.

convertDiamondAnnotationToAbund <- function(annFile){
  annDF <- read.delim(annFile, header = F, stringsAsFactors = F)
  foundIds <- annDF$V2
  listToTurnIntoMapDF <- lapply(foundIds,function(x){
    sepPoint1 <- stri_locate_last(x,regex ="_\\[")[1]-1
    sepPoint2 <- stri_locate_last(x,regex ="_\\[")[2] 
    funcPoint <- unlist(strsplit(substr(x,1,sepPoint1),"\\.1_"))[2]
    orgPoint <- gsub("\\[","",substr(x,sepPoint2,nchar(x)))
    orgPoint <- gsub("\\]","",orgPoint)
    c(funcPoint,orgPoint)
  })
  fullListDF <- data.frame(ldply(listToTurnIntoMapDF))
  names(fullListDF) <- c("function","organism")
  fullListDF[["function"]] <- gsub("MULTISPECIES:_","",fullListDF[["function"]])
    ### get the function and organism abundance
  funcAbundDF <- data.frame(table(fullListDF[["function"]]), 
                            stringsAsFactors = F,
                            row.names = names(table(fullListDF[["function"]])))
  funcAbundDF <- subset(funcAbundDF,select=-c(Var1))
  orgAbundDF <- data.frame(table(fullListDF[["organism"]]), 
                            stringsAsFactors = F,
                            row.names = names(table(fullListDF[["organism"]])))
  orgAbundDF <- subset(orgAbundDF,select=-c(Var1))  
  names(funcAbundDF) <- gsub(".RefSeq.annotated.txt","",basename(annFile))
  names(orgAbundDF) <- gsub(".RefSeq.annotated.txt","",basename(annFile))
  
  return(list(orgAbundDF = orgAbundDF,
              funcAbundDF = funcAbundDF,
              funcOrgMap = fullListDF))
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
metagMetatrDifferentialAnalysis_tableAndPlots <- function(countMatrix,designMatrix,
                                            group,sampleGroup,refGroup,
                                            mode="expression",DB=NULL,N=30){
  colnames(designMatrix) <- gsub(" \\[Factor\\]","",colnames(designMatrix))
  design <- model.matrix(~0+designMatrix[[group]])
  colnames(design) <- gsub("designMatrix\\[\\[group\\]\\]","",colnames(design))
  filteredMatrix <- filterCountsForDiffExpr(countMatrix,designMatrix,
                                            group,sampleGroup,refGroup)
  if (mode == "expression"){
    ddsTxi <- DESeqDataSetFromMatrix(filteredMatrix,
                                     colData = data.frame(designMatrix, row.names = rownames(designMatrix)),
                                     design =  design)
    resObj <- DESeq(ddsTxi)
    resDF <- data.frame(results(resObj))
    addTaxa <- resDF[order(resDF$padj),]
    addTaxa <- addTaxa[!is.na(addTaxa$padj),]
  colsToKeep <- c("log2FoldChange","padj")
  addTaxa <- addTaxa[,colsToKeep]
  ## select top 20 to report in the table

  tableToReport <- addTaxa[1:20,]
  
  ## mark significance
  addTaxa$Significance <- "Significant"
  addTaxa[addTaxa$padj > 0.05,]$Significance <- "nonSignificant"
  addTaxa$Significance <- as.factor(addTaxa$Significance)
  addTaxa$Significance <- factor(addTaxa$Significance, 
                                 levels = rev(levels(addTaxa$Significance)))

  ### volcano plot
  title <- "Volcano plot (padj threshold  = 0.05)"
  volcanoPlot <- ggplot(addTaxa, aes(y=-log10(padj), x=log2FoldChange)) +
    geom_point(aes(color=Significance),size=3) 
  volcanoPlot <- volcanoPlot + labs(title=title) + 
    theme(plot.title=element_text(size=10,hjust=0.5))
  volcanoPlot <- volcanoPlot + geom_hline(yintercept=1.3, color="blue") + 
    labs(color="Significance")
  return(list(vPlot=volcanoPlot,tableToReport=tableToReport))
  } else {
    ### p-value
    toGetValue <- function(x,y) {
        result  = round(t.test(x,y)$p.value, 3)
          return(result)
        }
    filteredMatrixLog <- log10(filteredMatrix+0.001)
    samplesGroup1 <- rownames(designMatrix[designMatrix[[group]]  == sampleGroup,])
    samplesGroup2 <- rownames(designMatrix[designMatrix[[group]]  == refGroup,])
    countGroup1 <- data.matrix(filteredMatrixLog[,samplesGroup1])
    countGroup2 <- data.matrix(filteredMatrixLog[,samplesGroup2])
    filteredMatrix$pvalue <- mapply(toGetValue, alply(countGroup1, 1),alply(countGroup2,1))
    filteredMatrixSort <- filteredMatrix[order(filteredMatrix$pvalue),]
    filteredMatrixSort <- head(filteredMatrixSort,N)
    return(filteredMatrixSort)
    ### Add barplots stacking org on x=samples for the top N diff. expr. func
  }
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
filterCountsForDiffExpr <- function(DF, desMat, groupName, group1, group2){
  samplesGroup1 <- rownames(desMat[desMat[[groupName]]  == group1,])
  samplesGroup2 <- rownames(desMat[desMat[[groupName]]  == group2,])
  countGroup1 <- DF[,samplesGroup1]
  countGroup2 <- DF[,samplesGroup2]
  isPresGroup1 <- apply(countGroup1,1,function(x) sum(x >5)>1)
  isPresGroup2 <- apply(countGroup2,1,function(x) sum(x >5)>1)
  subsetDF <- DF[isPresGroup1 |isPresGroup2,]
  subsetDF <- subsetDF[apply(subsetDF,1,sd)>0,]
  return(subsetDF)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Organism/function heatmap
##' @description Prepares organsim/function matrix for heatmap
##' @param  a list of annotated org/func abundance table
##' @return Returns a data.frame
orgFuncHeatmapPrep <- function(listOfAnnotatedAbundTable,diffTableFunc){
selectMap <- listOfAnnotatedAbundTable[["funcOrgMap"]]
topDiff <- grep("hypothetical_protein",rownames(diffTableFunc$tableToReport)[1:20],
                invert = T, value = T)
ff <- selectMap[selectMap[["function"]]%in%topDiff,]
ff1 <- lapply(unique(ff[["organism"]]), function(x){
                     pp <- data.frame(table(ff[ff$organism == x,"function"]))
                     rownames(pp) <- pp$Var1
                     pp <- subset(pp, select=-c(Var1))
                     names(pp) = x
                     return(pp)
                                })
finalDFforTopFuncAbundPlot <- listOfAbundMerge(ff1,unique(ff[["organism"]]))
return(finalDFforTopFuncAbundPlot) 
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Organism/function heatmap-to-boxplot preparation
##' @description Gathers organsim/function from heatmap to boxplot DF
##' @param  a heatmap-like dataframe
##' @return Returns a data.frame
orgFuncHeatmaptoBoxPlotPrep <- function(x,designMatrix,group,N){
  colnames(designMatrix) <- gsub(" \\[Factor\\]","",colnames(designMatrix))
  designMatrixSub <- data.frame(designMatrix[,group],
                                row.names =rownames(designMatrix))
  names(designMatrixSub) <- group
  orgDFforHeatmap2 <- cbind(x,sdV=apply(x,1,sd))
  orgDFforHeatmap2 <- head(orgDFforHeatmap2[order(orgDFforHeatmap2$sdV,decreasing = T),],N)
  orgDFforHeatmap2 <- subset(orgDFforHeatmap2,select=-c(sdV))
  orgDFforHeatmap2$meanV <- apply(orgDFforHeatmap2,1,mean)
  orgDFforHeatmap3 <-  orgDFforHeatmap2[order(orgDFforHeatmap2$meanV,decreasing = T),]
  orgDFforHeatmap3  <- subset(orgDFforHeatmap3 ,select=-c(meanV))
  orgDFforHeatmap3$ID <- rownames(orgDFforHeatmap3)
  boxPlotReadyDF <- gather(orgDFforHeatmap3,"sample","abund",-ID)
  dd <- merge(boxPlotReadyDF,designMatrixSub,by.x="sample",by.y=0,all=TRUE)
  return(dd) 
}
