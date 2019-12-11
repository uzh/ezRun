
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
  colToDrop <- c("Group","label","numOtus")
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
  colnames(taxaMatrix) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[1:(ncol(taxaMatrix))]
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
phyloSeqPreprocess <- function(phyloseqObj,rawCount,sampleFraction){
  ### First remove taxa not seen at least rowCount times in at least sampleFraction of the samples
filteredTaxa <- filter_taxa(phyloseqObj, function(x) sum(x > rawCount) > (sampleFraction*length(x)), TRUE)
  ### then remove samples which have zero observations
samplesToKeep <- which(apply(otu_table(filteredTaxa),1,sum)>0)
filteredTaxaAndSamples <- prune_samples(names(samplesToKeep),filteredTaxa)
  return(filteredTaxaAndSamples)
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
phyloSeqToDeseq2_tableAndPlots <- function(phyloseqObj,rank,group,sampleGroup,refGroup){
  ## Convert to Deseq obj and analyze
  ## to do: add selector for group and test
  phyloseqObjNoMock <- prune_samples(sample_data(phyloseqObj)[[group]] != "Mock", phyloseqObj)
  phyloseqToDeseq2Obj = phyloseq_to_deseq2(phyloseqObjNoMock,  as.formula(paste0("~",group)))
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(phyloseqToDeseq2Obj), 1, gm_mean)
  phyloseqToDeseq2Obj = estimateSizeFactors(phyloseqToDeseq2Obj, geoMeans = geoMeans)
  DEseqPhyRes <- DESeq(phyloseqToDeseq2Obj, test="Wald", fitType="parametric")
  res = results(DEseqPhyRes, cooksCutoff = FALSE,contrast = c(group,sampleGroup,refGroup))
  otuObj <- data.frame(t(otu_table(phyloseqObjNoMock)@.Data), check.names = F, stringsAsFactors = F)
  taxObj <-  data.frame(tax_table(phyloseqObjNoMock)@.Data, check.names = F,stringsAsFactors = F)
  taxObj[is.na(taxObj)] = "NA"
  addTaxaOut <- cbind(data.frame(res),otuObj,taxObj)
  addTaxaOut <- addTaxaOut[!is.na(addTaxaOut$Kingdom) & !is.na(addTaxaOut$padj),]
  addTaxaOut$id <- paste0("otu",seq(1,nrow(addTaxaOut)))
  addTaxaOut <- addTaxaOut[order(addTaxaOut$padj),]
  colsToKeep <- grep("baseMean|lfcSE", colnames(addTaxaOut), invert = T, value = T)
  addTaxa <- addTaxaOut[,colsToKeep]
  ## select fields to report in the table
  colsToRemove <- rownames(sample_data(phyloseqObjNoMock))
  colsToReport <- !colsToKeep%in%colsToRemove
  tableToReport <- addTaxa[1:20,colsToReport]
  
  ##
  addTaxa$Significance <- "Significant"
  addTaxa[addTaxa$padj > 0.05,]$Significance <- "nonSignificant"
  addTaxa$Significance <- as.factor(addTaxa$Significance)
  addTaxa$Significance <- factor(addTaxa$Significance, levels = rev(levels(addTaxa$Significance)))
  ## sort and prepare fpr plot
  x = tapply(addTaxa$log2FoldChange, addTaxa[[rank]], function(x) max(x))
  x = sort(x, TRUE)
  ### log2fold plot
  title <- "Abundance changes between the groups"
  plotLogFoldVsTaxon <- ggplot(addTaxa, aes(x=addTaxa[[rank]], y=log2FoldChange, color=addTaxa[[rank]])) + 
    geom_violin() + geom_point(size=3) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), axis.title.x = element_blank()) +
    geom_hline(aes(yintercept=1),color="blue")  + 
    geom_hline(aes(yintercept=-1),color="blue") 
  plotLogFoldVsTaxon <- plotLogFoldVsTaxon + labs(title=title) + 
    theme(plot.title=element_text(size=10,hjust=0.5)) + labs(color=rank)
  plotLogFoldVsTaxon <- plotLogFoldVsTaxon +guides(color = guide_legend(nrow = 10))
  ### volcano plot
  title <- "Volcano plot (p-value threshold  = 0.05)"
  volcanoPlot <- ggplot(addTaxa, aes(y=-log10(pvalue), x=log2FoldChange)) +
    geom_point(aes(shape=Significance, color=addTaxa[[rank]]),size=3) 
  volcanoPlot <- volcanoPlot + labs(title=title) + 
    theme(plot.title=element_text(size=10,hjust=0.5))
  volcanoPlot <- volcanoPlot + geom_hline(yintercept=1.3, color="blue") + labs(color=rank)
  ### Diff.expr. pie chart
  OTUsToPlot <- addTaxa[addTaxa$Significance == "Significant",]
  isAllNa <- all(names(table(OTUsToPlot[[rank]])) == "NA")
  if (isAllNa){
    isAllNaMsg <- paste("No differentially abundant  OTUs are annotated to the rank",rank, ". No pie chart to plot.")
    finalVersionPie <- NULL
  } else {
    isAllNaMsg <- NULL
  OTUsToPlot[[rank]]  <- as.factor(OTUsToPlot[[rank]])
  tableTaxa <- data.frame(table(droplevels(OTUsToPlot[,rank])))
  colnames(tableTaxa)[1] <- rank
  pct <- round(tableTaxa$Freq/sum(tableTaxa$Freq)*100,2)
  pct = paste0(pct,"%")
  titleText = paste("Distribution of differentially abundant taxa at rank", rank,".")
  bp <- ggplot(tableTaxa, aes(x="", y=Freq, fill=tableTaxa[[rank]])) + 
    geom_bar(position = position_stack(),width = 1, stat = "identity") 
  pieVersion <- bp + coord_polar("y", start=0)
  finalVersionPie <- pieVersion +  labs(title=titleText, y="") + 
    theme(plot.title=element_text(size=10,hjust=0.5)) + labs(fill=rank)
  }
  return(list(logPlot=plotLogFoldVsTaxon,vPlot=volcanoPlot,pieChart=finalVersionPie,tableToReport=tableToReport,
              fullTable=addTaxa,isAllNaMsg=isAllNaMsg,isAllNa=isAllNa))
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
    theme(legend.position="none",plot.title=element_text(size=15,hjust=0.5),  
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
pcaForPhyloseqPlot <- function(phySeqObject,type,group){
  ## calculate MDS values
  if (type =="samples") {
    input <- t(phySeqObject@otu_table@.Data)
    groups <- phySeqObject@sam_data[[group]]
  }else if (type =="taxa"){
    input <- phySeqObject@otu_table@.Data
    groups <- as.vector(phySeqObject@tax_table@.Data)
  } else {
    stop("type must be either samples or taxa")
  }

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
  g <- ggplot(ggg, aes(comp1,comp2, group = ggg[[group]])) + geom_point(aes(colour = ggg[[group]]),size =3)
  g <- g + xlab(xAxisLabel) + ylab(yAxisLabel)
  plot(g)
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
heatmapForPhylotseqPlotPheatmap <- function(phyloseqOtuObj,areThereMultVar,isGroupThere,rank){
  input <- data.frame(t(phyloseqOtuObj@otu_table@.Data), check.names = F)
  taxDF <- data.frame(phyloseqOtuObj@tax_table@.Data)
  input$rank <- taxDF[[rank]]
  colToAggregate <- grep("rank",colnames(input), value = T,invert = T)
  dd <- aggregate(x = input[colToAggregate],by = list(input$rank),sum)
  rownames(dd) <- dd$Group.1
  dd <-  subset(dd, select = -c(Group.1))
  dd <- dd[apply(dd,1,sd)>0,]
  input <- dd
  plot_heatmap_Pheatmap <- function() {
    if (isGroupThere){
  gr1 <- colnames(sample_data(phyloseqOtuObj))[1]
  gr2 <- colnames(sample_data(phyloseqOtuObj))[2]
  fact1 <- as.factor(sample_data(phyloseqOtuObj)@.Data[[1]])
  fact2 <- as.factor(sample_data(phyloseqOtuObj)@.Data[[2]])
  nColsGr1 <- nlevels(fact1)
  nColsGr2 <- nlevels(as.factor(sample_data(phyloseqOtuObj)@.Data[[2]]))
  pal1 <- colorRampPalette(brewer.pal(11, "Blues"))(nColsGr1)
  names(pal1) <- levels(fact1)
  pal2 <- colorRampPalette(brewer.pal(11, "RdYlGn"))(nColsGr2)
  names(pal2) <- levels(fact2)
  mat_colors <- list(pal1,pal2)
  names(mat_colors) <- c(gr1,gr2)
  mat_col_temp <- data.frame(sample_data(phyloseqOtuObj))
  mat_col <- data.frame(lapply(mat_col_temp,as.factor), row.names = rownames(mat_col_temp))
  if (!areThereMultVar){
    mat_colors <- mat_colors[gr1]
    mat_col <- data.frame(sample_data(phyloseqOtuObj)[,gr1])
  }
  ## heatmap
  pheatmap(input,show_rownames = TRUE,
           show_colnames     = TRUE,
           annotation_col    = mat_col,
           annotation_colors = mat_colors,
           cluster_rows = TRUE,  scale="row", method = "average")
    } else {

      pheatmap(input,show_rownames = TRUE,
               show_colnames     = TRUE,
               cluster_rows = TRUE,  scale="row", method = "average")
   }
  }
}  

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title OTUs saturation plot
##' @description HOw many OTUs do we really have?
##' @param  x, mothur shared abundance  file or already read-in table (dep on sec. param).
##' @return Returns a grid of plots
rarefactionPlot <- function(adundDF, type){
  if (type == 1){
    yLabel <- "Community saturation"
  }else{
    yLabel <- "Community rarefaction"
  }
  bb <- iNEXT(adundDF, q=0, datatype="abundance")
  fortifiedObj <- fortify(bb, type=type) 
  fortifiedObjPoint <- fortifiedObj[which(fortifiedObj$method=="observed"),]
  fortifiedObjLine <- fortifiedObj[which(fortifiedObj$method!="observed"),]
  fortifiedObjLine$method <- factor(fortifiedObjLine$method, 
                           c("interpolated", "extrapolated"),
                           c("int", "ext"))
  saturationPlot <- ggplot(fortifiedObj, aes(x=x, y=y, colour=site)) + 
    geom_point(aes(shape=site), size=4, data=fortifiedObjPoint) + 
    scale_shape_manual(values=rep(seq(1,23),3)) +
    geom_line(aes(linetype=method), lwd=1, data=fortifiedObjLine) +  
    guides(shape = guide_legend(nrow = 6), linetype= guide_legend(nrow = 2))

  saturationPlot <- saturationPlot + labs(x="Number of OTUs", 
                                          y=yLabel, 
                                          shape="Samples", colour="Samples",
                                          linetype="Method")
  return(saturationPlot)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title modified phyloseq barplot
##' @description Removed stacks border
##' @param  physeqFullObject (phyloseq object)
##' @return A ggplot
plotBarMod <- function(xx, x, fill = NULL, title = NULL, facet_grid = NULL,group) 
{
 mdf = psmelt(xx)
 if (x=="S"){
   xAxisVar ="Sample"
 mdf$relFract <- 0
 for (sample in unique(mdf$Sample)){
   tot <- sum(mdf[mdf$Sample == sample,]$Abundance)
 mdf[mdf$Sample == sample,]$relFract <- mdf[mdf$Sample == sample,]$Abundance/tot
 }
 }else if (x=="G") { 
   xAxisVar = group
   mdf$relFract <- 0
 for (g in levels(mdf[[group]])){
   tot <- sum(mdf[mdf[[group]] == g,]$Abundance)
   mdf[mdf[[group]] == g,]$relFract <-  mdf[mdf[[group]] == g,]$Abundance/tot
 }
 }
 p = ggplot(mdf, aes(x=mdf[[xAxisVar]],y=relFract, fill = mdf[[fill]]))
 p = p + geom_bar(stat = "identity", position = "stack") + xlab(xAxisVar) +
   ylab("relative fraction") + labs(fill = fill) + guides(fill = guide_legend(nrow = 12))
   p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
 if (!is.null(facet_grid)) {
     p <- p + facet_grid(facet_grid)
   }
   if (!is.null(title)) {
      p <- p + ggtitle(title)
}
return(p)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Rank-specific abundance plot
##' @description Abundance distribution for a specific rank
##' @param  physeqFullObject (phyloseq object),x (rank)
##' @return Returns a ggplot
abundPlot <- function(rank,physeqFullObject,xAesLogic,numTopRanks,group) {
  naRmoved <- subsetTaxMod(physeqFullObject, rank)
  if (naRmoved$toStop == TRUE) {
    return(list(abPlot=NULL,stop=TRUE))
  }else{
  naRmovedTrimmed <- subsetRankTopN(naRmoved$pObj, rank,numTopRanks)
  p <- plotBarMod(naRmovedTrimmed,x=xAesLogic, fill=rank,group=group)  
  p <- p+  theme(legend.key.size = unit(0.3, "cm"),legend.key.width = unit(0.3,"cm"))
  return(list(abPlot=p,stop=FALSE))
}
}
###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Modifed subset taxa
##' @description Removes NA for a specific rank for better plotting
##' @param  physeqFullObject (phyloseq object),x (rank)
##' @return Returns a filtered phyloseqobj
subsetTaxMod <- function (physeq, x) 
{
  if (is.null(tax_table(physeq))) {
    cat("Nothing subset. No taxonomyTable in physeq.\n")
    return(physeq)
  } else {
    oldMA <- as(tax_table(physeq), "matrix")
    oldDF <- data.frame(oldMA)
    newDF <- data.frame(oldDF[!is.na(oldDF[[x]]),])
    if(nrow(newDF) == 0){
      return(list(pObj=physeq,toStop=TRUE))
    } else{
    colnames(newDF) <- attr(physeq@tax_table@.Data, "dimnames")[[2]]
    newMA <- as(newDF, "matrix")
    if (inherits(physeq, "taxonomyTable")) {
      return(tax_table(newMA))
    } else {
      tax_table(physeq) <- tax_table(newMA)
      return(list(pObj=physeq,toStop=FALSE))
    }
  }
  }
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Rank-specific tax-based ordination plot
##' @description ordination plot for a specific rank
##' @param  fullObject (phyloseq object),x (rank)
##' @return Returns a ggplot
ordPlot <- function(rank,fullObject,type,areThereMultVar,numTopRanks,isGroupThere) {
  if (type=="taxa"){
    if (all(is.na(tax_table(fullObject)[,rank]))){
    naRmovedTrimmedOrd <- ordinate(fullObject, "NMDS", "bray")
    p1 = plot_ordination(fullObject, naRmovedTrimmedOrd, type = "taxa")    
    }else{
    naRmovedTrimmed <- subsetRankTopN(fullObject, rank,numTopRanks)
    naRmovedTrimmedOrd <- ordinate(naRmovedTrimmed, "NMDS", "bray")
    p1 = plot_ordination(naRmovedTrimmed, naRmovedTrimmedOrd, type = "taxa", color=rank)
    }
  }else if (type=="samples") {
    GP.ord <- ordinate(fullObject, "NMDS", "bray")
    if (isGroupThere) {
      gr1 <- colnames(sample_data(fullObject))[1]
     if (areThereMultVar){
     gr2 <- colnames(sample_data(fullObject))[2]
     p1 = plot_ordination(fullObject, GP.ord, type="samples", color=gr1,shape=gr2)
     }else{
     p1 = plot_ordination(fullObject, GP.ord, type="samples", color=gr1) + geom_point(size=6)
     }
    }  else {
      p1 = plot_ordination(fullObject, GP.ord, type="samples") + geom_point(size=6)
    }
  } else{
    stop("type must be either samples or taxa")
  }
  return(p1)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Select only top10 rank
##' @description subset taxa to retain only top10 renks. Useful for plotting
##' @param  physeqFullObject (phyloseq object),x (rank)
##' @return A subsetted phyloseq object
subsetRankTopN <- function(physeqFullObject,rank,N){
phylAsum = tapply(taxa_sums(physeqFullObject), tax_table(physeqFullObject)[, rank], sum, na.rm=TRUE)
topN = names(sort(phylAsum, TRUE))[1:N]
physeqFullObjectTrimmed = prune_taxa((tax_table(physeqFullObject)[, rank] %in% topN), physeqFullObject)
return(physeqFullObjectTrimmed)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Modified phyloseq richness plot 
##' @description Modified richness plot to include box plots statistical comparison when groups are present
##' @param  physeqFullObject (phyloseq object)
##' @return A ggplot
groupModRichPlot <- function(physeq, x, color = NULL, shape = NULL, 
          title = NULL, scales = "free_y", nrow = 1, shsi = NULL, 
          measures = c("Shannon"), sortby = NULL) 
{
  erDF = estimate_richness(physeq, split = TRUE, measures = measures)
  measures = colnames(erDF)
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  measures = measures[!measures %in% ses]
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    DF <- data.frame(erDF, sample_data(physeq))
  }else {
    DF <- data.frame(erDF)
  }
  if (!"samples" %in% colnames(DF)) {
    DF$samples <- sample_names(physeq)
  }
  if (!is.null(x)) {
    if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
      x <- "samples"
    }
  }else {
    x <- "samples"
  }
  mdf = reshape2::melt(DF, measure.vars = measures)
  p <- ggplot(mdf,aes(mdf[[x]],value))+ geom_boxplot() +
    stat_compare_means(method = "wilcox.test",hjust = -0.5, vjust = -0.5) +
    xlab(x)
  return(p)
}


###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Formats ggplot title
##' @description Centered and right size
##' @param  p (ggplot), text (title)
##' @return A plot.grid
add_centered_title <- function(p, text){
  grid.arrange(p, ncol = 1, top = text)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Summary plot of chimeras
##' @description Generate chimera plot from chimera file
##' @param  a data.frame
##' @return A ggplot
chimeraSummaryPlot <- function(chimeraDF){
bp <- ggplot(chimeraDF, aes(x=Type,Freq, fill=Type)) 
facetSampleBar <- bp  + geom_bar(stat = "identity",  position = 'dodge') 
finalVersionChimeraPlot  <- facetSampleBar +   
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) 
finalVersionChimeraPlot <- finalVersionChimeraPlot +
  geom_text(aes(y = Freq + 500, label = paste0(pct, '%')),
            position = position_dodge(width = .9),size = 3)
finalVersionChimeraPlot <- finalVersionChimeraPlot + facet_wrap(vars(sample))
finalVersionChimeraPlot <- finalVersionChimeraPlot + ylab("Count")
return(finalVersionChimeraPlot)
}

###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

##' @title Summary of community composition
##' @description It calculates the percentages of the organisms in the community
##' @param  a phyloseq object
##' @return A data.frame

communityPercSummTable <- function(phyloseqOtuObj,rank) {
  input <- data.frame(t(phyloseqOtuObj@otu_table@.Data))
  taxDF <- data.frame(phyloseqOtuObj@tax_table@.Data)
  input$rank <- taxDF[[rank]]
  colToAggregate <- grep("rank",colnames(input), value = T,invert = T)
  aggrDF <- aggregate(x = input[colToAggregate],by = list(input$rank),sum)
  rownames(aggrDF) <- aggrDF$Group.1
  aggrDF <-  subset(aggrDF, select = -c(Group.1))
  aggrDF <- aggrDF[apply(aggrDF,1,sd)>0,]
  percTable <- apply(aggrDF,2,function(y)sapply(y,function(x) round(x/sum(y)*100,2)))
  percTable <- data.frame(percTable)
  percTable$meanPc <- apply(percTable,1,mean)
  percTable <- percTable[order(percTable$meanPc, decreasing = T),]
  percTableToShow <- head(subset(percTable,  select = -c(meanPc)),20)
  return(percTableToShow)
}