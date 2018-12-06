
###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Phyloseq OTU object
##' @description Create Phyloseq OTU object mothur OTU files.
##' @param  a data frame in the format of mothur shared clustered OTU  files.
##' @return Returns a Phyloseq OTU object.

phyloSeqOTU <- function(otuDF){
rownames(otuDF) <- otuDF$Group
colToDrop <- c("Group")
otuFile1 <- as.matrix(otuDF[,!names(otuDF)%in%colToDrop])
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
##' @param  taxaDB, a DF in the format of mothur taxonomy file.
##' @return Returns a Phyloseq Taxa object.

phyloSeqTaxa <- function(taxaDB){
tempList <- lapply(taxaDB$Taxonomy,function(y) unlist(strsplit(y,";")))
taxaMatrix <- as.matrix(ldply(tempList))
rownames(taxaMatrix) <- taxaDB$OTU
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
  sampleFile <- sampleFileName
#colToKeep <- grep("Factor",colnames(sampleFile))
#colnames(sampleFile) <- sub('\\s\\[Factor\\]',"",colnames(sampleFile))
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
  GPfr <- phyloseqObj
  ## Standardize abundances to the median sequencing depth
 # total = median(sample_sums(phyloseqObj))
 # standf = function(x, t=total) round(t * (x / sum(x)))
#  gps = transform_sample_counts(phyloseqObj, standf)
#  ## transformed to relative abundance and filtered by abundance
#  GPr  = transform_sample_counts(gps, function(x) x / sum(x) )
#  GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
  filteredPhyloseqObj = filter_taxa(GPfr, function(x) sum(x > 3*1e-5) > (0.2*length(x)), TRUE)
  filteredPhyloseqObj = subset_taxa(filteredPhyloseqObj, Domain=="Bacteria")
 # filteredPhyloseqObj = subset_taxa(GPfr, Domain=="Bacteria")
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
  phyloseqObjNoMock <- prune_samples(sample_data(phyloseqObj)$Group != "Mock", phyloseqObj)
  phyloseq_to_deseq2 = phyloseq_to_deseq2(phyloseqObjNoMock, ~ Group)
  DEseqPhyRes <- DESeq(phyloseq_to_deseq2, test="Wald", fitType="parametric")
  res = results(DEseqPhyRes, cooksCutoff = FALSE)
  addTaxa <- cbind(data.frame(res),t(otu_table(phyloseqObjNoMock)), tax_table(phyloseqObjNoMock))
  addTaxaOut <- cbind(data.frame(res),t(otu_table(phyloseqObjNoMock)), tax_table(phyloseqObjNoMock))
  addTaxaOut <- addTaxaOut[order(addTaxaOut$padj),]
  addTaxaOut <- head(addTaxaOut,20)
  colsToKeep <- grep("baseMean|lfcSE", colnames(addTaxaOut), invert = T)
  addTaxaOut <- addTaxaOut[,colsToKeep]
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
  plotLogFoldVsTaxon <- ggplot(addTaxa, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
    geom_point(size=3) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
    geom_hline(aes(yintercept=1),color="red") + geom_text(aes(1,1,label = 1, vjust = -1), color = "red", size =3) + 
    geom_hline(aes(yintercept=-1),color="red") + geom_text(aes(1,-1,label = -1, vjust = 1), color = "red", size =3)
  plotLogFoldVsTaxon <- plotLogFoldVsTaxon + labs(title=title) + 
    theme(plot.title=element_text(size=10, face="bold",hjust=0.5))
  ### volcano plot
  title <- "Volcano plot (padj  = 0.05)"
  volcanoPlot <- ggplot(addTaxa, aes(y=-log10(pvalue), x=log2FoldChange)) +
    geom_point(aes(shape=Significance, color=Phylum),size=3) 
  volcanoPlot <- volcanoPlot + labs(title=title) + 
    theme(plot.title=element_text(size=10, face="bold",hjust=0.5))
  ### Diff.expr. pie chart
  OTUsToPlot <- na.omit(addTaxaOut[addTaxaOut$padj < 0.05,])
  tableTaxa <- data.frame(table(droplevels(OTUsToPlot[,"Genus"])))
  colnames(tableTaxa)[1] <- "Genus"
  pct <- round(tableTaxa$Freq/sum(tableTaxa$Freq)*100,2)
  pct = paste0(pct,"%")
  titleText = "Differentially abundant genera"
  bp <- ggplot(tableTaxa, aes(x="", y=Freq, fill=Genus)) + 
    geom_bar(position = position_stack(),width = 1, stat = "identity") 
  pieVersion <- bp + coord_polar("y", start=0)
  finalVersionPie <- pieVersion +  labs(title=titleText, y="") + 
    theme(plot.title=element_text(size=10, face="bold",hjust=0.5))
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

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title  Rank summary plots  
##' @description Summarized the top present genera  
##' @param  phySeq,rank a phyloseq object and the rank to summarize
##' @return Returns a stacked bar  plot.

phyloSeqCommunityComp <- function(physeq){
  sampleToKeep <- rownames(sample_data(physeq))[grep("mock",rownames(sample_data(physeq)))]
  mockObj <- prune_samples(sampleToKeep, physeq)
  otus <- data.frame(t(otu_table(mockObj)))
  joinedData <- data.frame(cbind(Freq=otus[,sampleToKeep],Genus=data.frame(tax_table(physeq), stringsAsFactors = FALSE)$Genus),
                           stringsAsFactors = FALSE)
  joinedData$Freq <- as.numeric(joinedData$Freq)
  fractions <- lapply(unique(joinedData$Genus), 
                             function(x) cbind(x,sum(joinedData[joinedData$Genus == x,]$Freq)/sum(joinedData$Freq)*100))
  fractions <- data.frame(matrix(unlist(fractions), byrow = TRUE, ncol = 2), stringsAsFactors = FALSE)
  colnames(fractions) <- c("Genus","Freq")
  fractions$Freq <- round(as.numeric(fractions$Freq),2)
  fractions <- fractions[fractions$Freq >2 ,]
  titleText <- "Community composition (Genera, min. 2 %)"
  bp <- ggplot(fractions, aes(x="", y=Freq, fill=Genus)) +  
    geom_bar(position = position_stack(),width = 1, stat = "identity") + 
    geom_text(aes(label = paste0(Freq,"%")), position = position_stack(vjust = 0.5),  size = 3)
  finalVersion <- bp +  labs(title=titleText, y="") + 
    theme(plot.title=element_text(size=15, face="bold",hjust=0.5))
  return(finalVersion)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title  Alternative pca plot for phyloseq abundnce-taxonimy matrix  
##' @description Alternative pca plot for phyloseq abundnce-taxonimy matrix  
##' @param  phySeq,rank a phyloseq otu object and the number or top ranked categpries
##' @return Returns an MDS   plot.
##' 
### PCA plot function
pcaForPhylotseqPlot <- function(input,groups){
  ## calculate MDS values
  input <- t(input)
  mds = plotMDS(log2(input+10), plot=FALSE, ndim=2)
  ggg <- data.frame(comp1 = mds$cmdscale.out[,1] ,comp2 = mds$cmdscale.out[,2], group = groups)
  
  ## calculate var explained
  s <- svd(input-rowMeans(input))
  cc <- s$d^2/sum(s$d^2)
  PC1varExpl <- round(cc[1],4)*100
  PC2varExpl <- round(cc[2],4)*100
  xAxisLabel <- paste0("PC1 (", PC1varExpl,"% explained var.)")
  yAxisLabel <- paste0("PC2 (", PC2varExpl ,"% explained var.)")
  
  ## plot
  g <- ggplot(ggg, aes(comp1,comp2, group = group)) + geom_point(aes(colour = group),size =3)
  g <- g + xlab(xAxisLabel) + ylab(yAxisLabel)
  plot(g)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title  Alternative heatmap plot for phyloseq abundnce-taxonimy matrix  
##' @description Alternative heatmap plot for phyloseq abundnce-taxonimy matrix  
##' @param   a phyloseq object and the rank to summarize
##' @return Returns a stacked bar  plot.
##' 
### Heatmap function
heatmapForPhyloseqPlot <- function(phyloseqOtuObj){
  plot_heatmap <- function() {
    ## clust funct
  distCor <- function(x) {as.dist(1-cor(x))}
  zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {
    if (scale=="row") z <- t(scale(t(x)))
    if (scale=="col") z <- scale(x)
    z <- pmin(pmax(z, zlim[1]), zlim[2])
    hcl_row <- hclust(distCor(t(z)), method=method)
    hcl_col <- hclust(distCor(z), method=method)
    return(list(data=z, hcl_r=hcl_row,hcl_c=hcl_col, 
                Rowv=as.dendrogram(hcl_row), Colv=as.dendrogram(hcl_col)))
  }
  z <- zClust(t(phyloseqOtuObj))
  cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
  ## heatmap
    heatmap.2(z$data,dendrogram=c("col"),Rowv=FALSE,Colv=z$Colv,col=rev(cols), 
              trace='none',density.info=c("none"),keysize = 0.8, 
              labRow=NA,cexCol = 1)
  }
}


###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Phyloseq OTU object
##' @description Create Phyloseq OTU object mothur OTU files.
##' @param  a data frame in the format of mothur shared clustered OTU  files.
##' @return Returns a Phyloseq OTU object.

phyloSeqOTUFromFile <- function(otuFile){
  otuDF <- read.delim(otuFile, header = T,stringsAsFactors = F, check.names = F)
  rownames(otuDF) <- otuDF$Group
  colToDrop <- c("Group","label","numOTUs")
  otuFile1 <- as.matrix(otuDF[,!names(otuDF)%in%colToDrop])
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
##' @param  taxaDB, a DF in the format of mothur taxonomy file.
##' @return Returns a Phyloseq Taxa object.

phyloSeqTaxaFromFile  <- function(taxaFile){
  taxaDB <- read.delim(taxaFile, header = T,stringsAsFactors = F)
  tempList <- lapply(taxaDB$Taxonomy,function(y) unlist(strsplit(y,";")))
  taxaMatrix <- as.matrix(ldply(tempList))
  rownames(taxaMatrix) <- taxaDB$OTU
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


##' @title  Alternative heatmap with pheatmap plot for phyloseq abundnce-taxonimy matrix  
##' @description Alternative heatmap plot for phyloseq abundnce-taxonimy matrix  
##' @param   a phyloseq object and the rank to summarize
##' @return Returns a stacked bar  plot.
##' 
### Heatmap function
heatmapForPhylotseqPlotPheatmap <- function(phyloseqOtuObj, matrix){
    plot_heatmap_Pheatmap <- function() {
  ## clust funct
  distCor <- function(x) {as.dist(1-cor(x))}
  zClust <- function(x, scale="row", zlim=c(-3,3), method="average") {
    if (scale=="row") z <- t(scale(t(x)))
    if (scale=="col") z <- scale(x)
    z <- pmin(pmax(z, zlim[1]), zlim[2])
    hcl_row <- hclust(distCor(t(z)), method=method)
    hcl_col <- hclust(distCor(z), method=method)
    return(list(data=z, hcl_r=hcl_row,hcl_c=hcl_col, 
                Rowv=as.dendrogram(hcl_row), Colv=as.dendrogram(hcl_col)))
  }
  z <- zClust(t(phyloseqOtuObj))
  mat_col <- matrix
  ncols <- min(3,nlevels(as.factor(mat_col$Group)))
  mat_colors <- list(group = c("red","blue"))
  names(mat_colors$group) <- unique(mat_col$Group)
  ## heatmap
  pheatmap(z$data,show_rownames = FALSE,
           show_colnames     = TRUE,
           annotation_col    = mat_col,
           annotation_colors = mat_colors,
           cluster_rows = TRUE,  scale="column", method = "average")
   }
}
