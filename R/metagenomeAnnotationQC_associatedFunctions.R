
###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Prepares all the prodigal-associated files
##' @description Prepares all the prodigal-associated plots
##' @param  a pridgal.gff prediction file
##' @return Returns the full DF and the subest to partial =00

### get input files
prodigalFileReport <- function(x,meth){
prodigalGffImport <- import.gff(x)
prodigalSummaryDF <- data.frame(mcols(prodigalGffImport), stringsAsFactors = F)
prodigalSummaryDF$gc_cont <- as.numeric(prodigalSummaryDF$gc_cont)
prodigalSummaryDF$conf <- as.numeric(prodigalSummaryDF$conf)
prodigalSummaryDF$method <- meth
subsetDataToPartial00DF <- prodigalSummaryDF[prodigalSummaryDF$partial =="00" 
                                           & (prodigalSummaryDF$start_type == "ATG"|
                                                prodigalSummaryDF$start_type == "GTG"),]
subsetDataToPartial00DF$method <- meth
return(list(fullSumm=prodigalSummaryDF,subsetDataToPartial00=subsetDataToPartial00DF))
}




###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Prepares all the interproscan-associated files; the sceond extracts topN from the list
##' @description Prepares all the interproscan-associated plots
##' @param  a interproscan.gff prediction file
##' @return Returns the full DF
## extract N entries with top frequency 
extractTopN <- function(DF,column,N){
  col <- vector()
  tabNoNa <- DF[DF[[column]] != "NA",]
  tab <- table(tabNoNa[[column]])
  tab_s <- sort(tab)                                           
  col <- data.frame(tail(names(tab_s), N), stringsAsFactors = F)
  colnames(col) <- column
  topN <- data.frame(cbind(col, abundance = tail(as.data.frame(tab_s)$Freq, N)),
                     stringsAsFactors = F)
  topN <- topN[order(topN$abundance), ]
  topN[[column]] <- gsub("\"","",topN[[column]])
  return(topN)
}

interproscanFileReport <- function(x,meth,N){
IPSGffImport <- import.gff(x)
description <- mcols(IPSGffImport)$signature_desc
description[sapply(description,function(x) length(x)==0)] <- "NA"
description <- sapply(description,function(x)x[1])
ontology <- mcols(IPSGffImport)$Ontology_term
ontology[sapply(ontology,function(x) length(x)==0)] <- "NA" 
ontology <- sapply(ontology,function(x)x[1])
IPSGffSummaryDF <- data.frame(score = as.numeric(mcols(IPSGffImport)$score),
                              description = description, 
                              GOterm = ontology,
                              type = mcols(IPSGffImport)$type,
                              stringsAsFactors = F)
IPSGffSummaryDF <- IPSGffSummaryDF[IPSGffSummaryDF$type == "protein_match",
                                   c("score","description","GOterm")]
IPSGffSummaryDF_topN_GO <- extractTopN(IPSGffSummaryDF,"GOterm",N)
IPSGffSummaryDF_topN_desc <- extractTopN(IPSGffSummaryDF,"description",N)
full_GO <- extractTopN(IPSGffSummaryDF,"GOterm",nrow(IPSGffSummaryDF))
full_descrip <- extractTopN(IPSGffSummaryDF,"description",nrow(IPSGffSummaryDF))
IPSGffSummaryDF$method <- meth
IPSGffSummaryDF_topN_GO$method <- meth
IPSGffSummaryDF_topN_desc$method <- meth
  return(list(summDF=IPSGffSummaryDF,
              topN_GO=IPSGffSummaryDF_topN_GO,
              topN_desc=IPSGffSummaryDF_topN_desc,
              full_GO=full_GO,
              full_descrip=full_descrip))
}


###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Prepares all the prodigal-associated plots
##' @description Prepares all the prodigal-associated plots
##' @param  a pridgal.gff prediction file
##' @return Returns ggplots

summaryScorePlot <- function(x){
p <-  ggplot(x,aes(x=score,fill=method)) + geom_histogram(binwidth=10) +  
  facet_grid(rows = vars(start_type), cols = vars(partial)) + 
  labs(title="Summary of the gene prediction scores")
return(p)
}

summaryConfPlot <- function(x){
  p <-  ggplot(x,aes(x=conf,fill=method)) + geom_histogram(binwidth=5, position = "dodge") +  
    facet_grid(rows = vars(start_type), cols = vars(partial),labeller = label_both)+ 
    labs(title="Summary of the gene prediction confidence")
  return(p)
  }

summaryGcContPlot <- function(x){
  p <- ggplot(x,aes(x=gc_cont,fill=method)) + geom_histogram(binwidth=0.001) +  
    facet_grid(rows = vars(start_type), cols = vars(partial),labeller = label_both) + labs(title="Summary of the GC-content")
  return(p)
}

### summary rbs_spacer hist
summaryRBSSpacePlot <-  function(x){
  p<- ggplot(x,aes(x=rbs_spacer,fill=method)) + geom_bar(position = "dodge") +  
    facet_grid(rows = vars(start_type), cols = vars(partial), labeller = label_both) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="RBS-spacer distribution")
  return(p)
}
### summary rbs_motif hist
summaryRBSMotifPlot <- function(x){
  p<- ggplot(x,aes(x=rbs_motif,fill=method)) + 
    geom_bar(position = "dodge") +  facet_grid(rows = vars(start_type), cols = vars(partial), labeller = label_both) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="RBS-motif distribution")
  return(p)
} 


###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Prepares all the IPS-associated plots
##' @description Prepares all the prodigal-associated plots
##' @param  a pridgal.gff prediction file
##' @return Returns ggplots

summaryMatchScorePlot <- function(x){
  p<- ggplot(x,aes(x=-log(score), fill=sample)) + geom_histogram(binwidth = 1) +
    labs(title="Summary of protein match score") + facet_wrap(vars(sample)) +
    theme(legend.position = "none")
  return(p)
}

### topNcateg plots: GO
summaryGOPlot <- function(x,numberOfTopNCategories){
  DFforSummProtGO <- x
  nEntries <- nrow(DFforSummProtGO)
  if (nEntries < 25){
    yLabelSize = 5
  }else if (nEntries > 25 & nEntries <31) {
    yLabelSize = 4
  }else if ( nEntries >31){
    yLabelSize = 3
  }
GOdesc <- sapply(DFforSummProtGO$GOterm,
                 function(x)Term(GOTERM)[names(Term(GOTERM))%in%x], USE.NAMES = F)
DFforSummProtGO$GOterm <- paste(DFforSummProtGO$GOterm,GOdesc)
## expand palette colour to numberOfTopNCategories
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
expandedPalette <- getPalette(numberOfTopNCategories)
basicPlot <- ggplot(DFforSummProtGO,aes(x=reorder(GOterm, -abundance),y=abundance, fill = method)) 
withThemes <- basicPlot + geom_bar(stat = "Identity", position = "dodge") +
  theme(axis.title = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y= element_text(size =yLabelSize),
        legend.position = "right",legend.text = element_text(size =10)) 
p <- withThemes  +scale_color_manual(expandedPalette) +
  labs(title="Most represented GO terms") +
  guides(fill=guide_legend(ncol=1, byrow=F,title.position = "top"))+
  aes(stringr::str_wrap(GOterm,30)) + ylab(NULL) + coord_flip()
  return(p)
}
### topNcateg plots: protein family
summaryFamilyPlot <- function(x,numberOfTopNCategories){
  DFforSummProtFamilies <- x
  nEntries <- nrow(DFforSummProtFamilies)
  if (nEntries < 25){
    yLabelSize = 5
  }else if (nEntries > 25 & nEntries <31) {
    yLabelSize = 4
  }else if ( nEntries >31){
    yLabelSize = 3
  }
  ## expand palette colour to numberOfTopNCategories
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  expandedPalette <- getPalette(numberOfTopNCategories)
  basicPlot <- ggplot(DFforSummProtFamilies,aes(x=reorder(description, -abundance),y=abundance, fill = method)) 
  withThemes <- basicPlot + geom_bar(stat = "Identity", position = "dodge") +
    theme(axis.title = element_blank(),axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y= element_text(size =yLabelSize),
          legend.position = "right",legend.text = element_text(size =10)) 
  p <- withThemes  +scale_color_manual(expandedPalette) +
    labs(title="Most represented  protein families") +
    guides(fill=guide_legend(ncol=1, byrow=F,title.position = "top"))+
    aes(stringr::str_wrap(description,30)) + ylab(NULL) + coord_flip()
  
return(p)
}



###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title summaryBins
##' @description Create summary files for bins
##' @param  a kraken.labels files and fasta bin files 
##' @return Returns a DF summarizing bins 

### get input files
summaryMetagenomeBins <- function(krakenFile,binFiles){
  taxRanks <- c("Life","Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  krakenLabels <- read.delim(krakenFile, stringsAsFactors = F, header = F)
  names(krakenLabels) <- c("contigID","orgName")
  binIDmapList <- lapply(binFiles,function(x) {
    y <- readDNAStringSet(x)
    data.frame(contigID=names(y), 
               binID=unlist(strsplit(basename(x),"\\."))[2],
               stringsAsFactors = F)
               })
  binIDmap <- do.call("rbind",binIDmapList)
  taxonSummary <- lapply(krakenLabels$orgName,
                                 function(x){
                                   y <- unlist(strsplit(x,";"))
                                   y <- c(y[1:2],tail(y,7))
                                   taxLevel <- min(length(y),length(taxRanks))
                                   orgTemp <- gsub("\\[","",y[length(y)])
                                   org <- gsub("\\]","",orgTemp)
                                   return(list(rank=taxRanks[taxLevel], organism=org))}
                                 )
  taxonSummaryDF <- ldply(taxonSummary,unlist)
  krakenLabels$coverage <- sapply(krakenLabels$contigID, function(x){
    round(as.numeric(unlist(strsplit(x,"cov_"))[2]),2)
  })
  krakenLabels$ID <- sapply(krakenLabels$contigID, function(x){
    unlist(strsplit(x,"_length"))[1]
  })
  krakenLabels$length<- sapply(krakenLabels$contigID, function(x){
    as.numeric(unlist(strsplit(x,"_"))[4])
  })
  fullDF <- subset(cbind(krakenLabels,taxonSummaryDF),select=-c(orgName))
  mergingBinTax <-  merge(fullDF,binIDmap,by="contigID", all=FALSE)
  finalDF <- subset(mergingBinTax, select=-c(contigID))
  return(finalDF)
}


###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title mergeSummaryBinFile 
##' @description Merged summary bin files for plots
##' @param  a list of summary bin files generated from the summaryMetagenomeBins fun
##' @return Returns a data frame
  
mergeSummaryBinFiles <- function(binFile){
  binDF <- read.delim(binFile, stringsAsFactors = F)
  binDF$sample <- as.factor(gsub(".binSummaryFile.txt","",binFile))
  binDF$binID <- as.factor(binDF$binID)
  return(binDF)
}


###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title createAbundTable 
##' @description It creates abundance table from the organims annotations to bins
##' @param  the data frame generated by mergeSummaryBinFiles
##' @return Returns a data frame
createAbundTable <- function(mergedSummaryBinDF){
  DFforHeatmap <- list()
  allSamples <- levels(mergedSummaryBinDF$sample)
  for (sample in allSamples){
    sampleSpecDF <- mergedSummaryBinDF[mergedSummaryBinDF$sample == sample,]
    DFforHeatmapTemp <- data.frame(table(sampleSpecDF$organism))
    rownames(DFforHeatmapTemp) <- DFforHeatmapTemp$Var1
    DFforHeatmapTemp <- subset(DFforHeatmapTemp,select=-c(Var1)) 
    DFforHeatmap[[sample]] <- data.frame(DFforHeatmapTemp)
  }
  finaLMergedAbundTable <- listOfAbundMerge(DFforHeatmap,names(DFforHeatmap))
  return(finaLMergedAbundTable)
}
###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title makeHeatmapFromSummbinFile 
##' @description taxa-by-sample heatmap
##' @param  the data frame generated by mergeSummaryBinFiles
##' @return Returns a heatmap

summaryHeatmap <- function(mergedBinAbundTable,isGroupThere,
                           dataset=NULL,numberOfTopNCategories,plotTitle="plotTitle"){
  plot_heatmap_Pheatmap <- function() {
    mergedBinAbundTable$sdV <-   apply(mergedBinAbundTable,1,sd)
    mergedBinAbundTableOrdered <- mergedBinAbundTable[order(mergedBinAbundTable$sdV,
                                                            decreasing = T),]
    mergedBinAbundTableNoSdV <- subset(mergedBinAbundTableOrdered, select=-c(sdV))
    fontsize = 8
    if (isGroupThere) {
      mergedBinAbundTableTopNToPlot <- head(mergedBinAbundTableNoSdV,
                                            numberOfTopNCategories)
      cellwidth = max(1,25-(nrow(dataset)))
    } else {
      mergedBinAbundTableTopN <- head(mergedBinAbundTableNoSdV,
                                            numberOfTopNCategories)
      mergedBinAbundTableTopNTran <- data.frame(t(mergedBinAbundTableTopN))
      rowsToKeep <- which(apply(mergedBinAbundTableTopNTran,1,sd) >0)
      mergedBinAbundTableTopNTran <- mergedBinAbundTableTopNTran[rowsToKeep,]
    mergedBinAbundTableTopNToPlot <- mergedBinAbundTableTopNTran
    cellwidth = max(1,25-(numberOfTopNCategories))
    }

    if (isGroupThere){
      gr <- list()
      nCols <- list()
      pal <- list()
         colsToKeep <- grep("Factor",colnames(dataset), value = T)
      for (col in colsToKeep){
       factVar <- as.factor(dataset[[col]])
      gr[[col]] <- gsub(" \\[Factor\\]","",col)
      nCols[[col]] <- nlevels(factVar)
      pal[[col]] <- colorRampPalette(brewer.pal(11, "Blues"))(nCols[[col]])
      names(pal[[col]]) <- levels(factVar)
      }
      mat_colors <- pal
      names(mat_colors) <- gr
      mat_col <- data.frame(dataset[,colsToKeep], row.names = rownames(dataset))
      names(mat_col) <- gr
      ## heatmap
      pheatmap(mergedBinAbundTableTopNToPlot,show_rownames = TRUE,
               show_colnames     = TRUE,
               annotation_col    = mat_col,
               annotation_colors = mat_colors,
               cluster_rows = FALSE, 
               cluster_cols = TRUE, 
               scale="row", 
               method = "average",
               cellwidth=cellwidth,
               fontsize = fontsize,
               cellheight = fontsize)
    } else {
     pheatmap(mergedBinAbundTableTopNToPlot,show_rownames = FALSE,
               show_colnames     = TRUE,
               cluster_rows = TRUE,  
               scale="row",
               method = "average",
               cellwidth=cellwidth,
               fontsize=fontsize,main=plotTitle)
    }
  }
}


