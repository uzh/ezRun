
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
otuFile1 <- as.matrix(otuFile[,!names(otuFile)%in%colToDrop])
otuObject <- otu_table(otuFile1, taxa_are_rows = FALSE)
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

phyloSeqTaxa <- function(taxaFileName){
taxaFile <- read.table(taxaFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
tempList <- lapply(taxaFile$Taxonomy,function(y) unlist(strsplit(y,";")))
taxaMatrix <- as.matrix(ldply(tempList))
rownames(taxaMatrix) <- taxaFile$OTU
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
##' @param  path to the tsv desing file 
##' @return Returns a Phyloseq Taxa object.
phyloSeqSample <- function(sampleFileName){
#sampleFile <- ezRead.table(sampleFileName)
  sampleFile <- sampleFileName
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


##' @title Phyloseq preprocess
##' @description Preprocesses a phyloseq object.
##' @param  phyloseqObj, a phyloseq object.
##' @return Returns a  filtered Phyloseq  object.
phyloSeqPreprocess <- function(phyloseqObj){
  ## Standardize abundances to the median sequencing depth
  total = median(sample_sums(phyloseqObj))
  standf = function(x, t=total) round(t * (x / sum(x)))
  gps = transform_sample_counts(phyloseqObj, standf)
  ## transformed to relative abundance and filtered by abundance
  GPr  = transform_sample_counts(gps, function(x) x / sum(x) )
  GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
  filteredPhyloseqObj = filter_taxa(GPfr, function(x) sum(x > 3*1e-5) > (0.2*length(x)), TRUE)
  filteredPhyloseqObj = subset_taxa(filteredPhyloseqObj, Domain=="Bacteria")
  filteredPhyloseqObj = subset_taxa(gps, Domain=="Bacteria")
  return(filteredPhyloseqObj)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Differential abundance analysis between groups
##' @description Comaprison of metagenomics communities from files stored as a phyloseq object.
##' @param  phyloseqObj, a phyloseq object.
##' @return Returns a list of tables and plots.
phyloSeqToDeseq2_tableAndPlots <- function(phyloseqObj){
  ## Convert to Deseq obj and analyze
  ## to do: add selector for group and test
  phyloseqObjNoMock <- prune_samples(sample_data(phyloseqObj)$Group != "Mock")
  phyloseq_to_deseq2 = phyloseq_to_deseq2(phyloseqObj, ~ Group)
  DEseqPhyRes <- DESeq(phyloseq_to_deseq2, test="Wald", fitType="parametric")
  res = results(DEseqPhyRes, cooksCutoff = FALSE)
  addTaxa <- cbind(data.frame(res),t(otu_table(phyloseqObj)), tax_table(phyloseqObj))
  addTaxaOut <- cbind(data.frame(res),t(otu_table(phyloseqObj)), tax_table(phyloseqObj))
  ## sort and prepare fpr plot
  x = tapply(addTaxa$log2FoldChange, addTaxa$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  addTaxa$Phylum = factor(as.character(addTaxa$Phylum), levels=names(x))
  # Genus order
  x = tapply(addTaxa$log2FoldChange, addTaxa$Genus, function(x) max(x))
  x = sort(x, TRUE)
  addTaxa$Genus = factor(as.character(addTaxa$Genus), levels=names(x))
  addTaxa <- na.omit(addTaxa)
  addTaxa$Significance <- "Significant"
  addTaxa[addTaxa$padj > 0.05,]$Significance <- "nonSignificant"
  addTaxa$Significance <- as.factor(addTaxa$Significance)
  addTaxa$Significance <- factor(addTaxa$Significance, levels = rev(levels(addTaxa$Significance)))
  ### log2fold plot
  title <- "Abundance changes between the groups"
  plotLogFoldVsTaxon <- ggplot(addTaxa, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
    geom_hline(aes(yintercept=1),color="red") + geom_text(aes(1,1,label = 1, vjust = -1), color = "red", size =3) + 
    geom_hline(aes(yintercept=-1),color="red") + geom_text(aes(1,-1,label = -1, vjust = 1), color = "red", size =3)
  plotLogFoldVsTaxon <- plotLogFoldVsTaxon + labs(title=title) + theme(plot.title=element_text(size=15, face="bold",hjust=0.5))
  ### volcano plot
  title <- "Volcano plot (padj  = 0.05)"
  volcanoPlot <- ggplot(addTaxa, aes(y=-log10(pvalue), x=log2FoldChange)) +
    geom_point(aes(shape=Significance, color=Phylum),size=3) 
  volcanoPlot <- volcanoPlot + labs(title=title) + theme(plot.title=element_text(size=15, face="bold",hjust=0.5))
  ### Diff.expr. pie chart
  OTUsToPlot <- na.omit(addTaxa[addTaxa$padj < 0.05,])
  tableTaxa <- data.frame(table(OTUsToPlot[,"Genus"]))
  colnames(tableTaxa)[1] <- "Genus"
  colRain=rainbow(nrow(tableTaxa))
  titleText = "Differentially abundant genera"
  bp <- ggplot(tableTaxa, aes(x="", y=Freq, fill=Genus)) + geom_bar(width = 10, stat = "identity") + 
    scale_fill_manual(values=colRain)
  pieVersion <- bp + coord_polar("y", start=0)
  finalVersionPie <- pieVersion +  labs(title=titleText, y="") + 
    theme(plot.title=element_text(size=15, face="bold",hjust=0.5))
  
  return(list(logPlot=plotLogFoldVsTaxon,vPlot=volcanoPlot,pieChart=finalVersionPie,table=addTaxaOut))
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Richness plot for a certain taxa 
##' @description Create pie chart richness plot 
##' @param  taxaFileName,rank mothur taxonomy file.
##' @return Returns a pie chart plot.

phyloSeqDivPlotAndPercUnclassified <- function(taxaFileName, rank){
  taxaFile <- read.table(taxaFileName, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  tempList <- lapply(taxaFile$Taxonomy,function(y) unlist(strsplit(y,";")))
  taxaMatrix <- as.matrix(ldply(tempList))
  rownames(taxaMatrix) <- taxaFile$OTU
  colnames(taxaMatrix) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")[1:(ncol(taxaMatrix))]
  taxaDF <- data.frame(taxaMatrix)
  tableTaxa <- data.frame(table(taxaDF[,rank]))
  colnames(tableTaxa)[1] <- rank
  labs <- as.character(tableTaxa[order(tableTaxa[,rank]),][,rank][1:3])
  colRain=rainbow(nrow(tableTaxa))
  percUnclassified <- length(grep("unclassified",taxaDF[,rank]))/length(taxaDF[,rank])*100
  percUnclassified <- paste(round(percUnclassified,2),"%", sep = " ")
  titleText = paste("Community diversity at the", rank, "level", sep = " ")
  subtitleText = paste("Percentage of unclassified", rank, "=", percUnclassified, sep = " ")
  bp <- ggplot(tableTaxa, aes(x="", y=Freq, fill=get(rank))) + geom_bar(width = 10, stat = "identity") + 
    scale_fill_manual(values=colRain)
  pieVersion <- bp + coord_polar("y", start=0)
  finalVersion <- pieVersion +  labs(title=titleText, subtitle=subtitleText, y="") + 
    theme(legend.position="none",plot.title=element_text(size=15, face="bold",hjust=0.5),  
          plot.subtitle=element_text(size=10, face="bold",hjust=0.5))
  return(finalVersion)
}