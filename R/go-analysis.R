###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


doGo = function(param, seqAnno){
  flag = param$runGO && (hasGeneMapping(param, seqAnno) || param$featureLevel == "gene") && hasGoAnnotation(seqAnno)
  message("doGo: ", flag)
  return(flag)
}

hasGoAnnotation = function(seqAnno){
  goColumns <- c("GO BP", "GO MF", "GO CC")
  hasGO <- all(goColumns %in% colnames(seqAnno))
  validGO <- all(colSums(seqAnno[, goColumns] == "") != nrow(seqAnno))
  validGO2 <- !any(is.na(seqAnno[, goColumns]))
  return(hasGO && validGO && validGO2)
}

##' @title Separates GO ID's by ontology
##' @description Separates GO ID's by ontology
##' @param goIdStrings the GO ID's to separate.
##' @template roxygen-template
##' @return Returns the separated GO ID's
##' @examples separateGoIdsByOnto(c("GO:0008150", "GO:0005575"))
separateGoIdsByOnto = function(goIdStrings){
  require("GO.db", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  bpGos = keys(GOBPPARENTS)
  mfGos = keys(GOMFPARENTS)
  ccGos = keys(GOCCPARENTS)
  result = data.frame(row.names=1:length(goIdStrings))
  goList = strsplit(trimWhiteSpace(gsub("; ", ";", goIdStrings)), ";") ## support both separators "; " and ";" 
  result[["GO BP"]] = sapply(lapply(goList, function(x,y){intersect(x,y)}, y=bpGos), ezCollapse)
  result[["GO MF"]] = sapply(lapply(goList, function(x,y){intersect(x,y)}, y=mfGos), ezCollapse)
  result[["GO CC"]] = sapply(lapply(goList, function(x,y){intersect(x,y)}, y=ccGos), ezCollapse)
  return(result)
}

##' @title Gets the GO parents
##' @description Gets the GO parents.
##' @param id the ID to check
##' @param onto the ontology to use
##' @template roxygen-template
##' @return Returns the ancestor terms if they are specified.
##' @examples
##' getGOparents("GO:0034767")
##' addGoParents(c("GO:0034767", "GO:0034768"), "BP")
getGOparents = function(id, onto="BP"){
  require(GO.db)
  go2Ancestor = switch(onto, BP=as.list(GOBPANCESTOR), MF=as.list(GOMFANCESTOR), CC=as.list(GOCCANCESTOR))
  ancestorIds = setdiff(go2Ancestor[[id]], "all") ## remove the all category
  if (is.null(ancestorIds)){
    return()
  } else {
    ancestorTerms = Term(GOTERM[ancestorIds])
    return(ancestorTerms)
  }
}

##' @describeIn getGOparents Adds the GO parents.
addGoParents = function(gene2goList, onto){
  require(GO.db)
  goParents = switch(onto, BP=as.list(GOBPPARENTS),
                     CC=as.list(GOCCPARENTS),
                     MF=as.list(GOMFPARENTS))
  return(lapply(gene2goList, function(x){setdiff(union(x, unlist(goParents[x])), "all")}))
}

addGeneNamesEnrich <- function(resEnrich, se){
  gene_ids <- strsplit(resEnrich$Genes, "; ")
  ids2names <- setNames(rowData(se)$gene_name, rowData(se)$gene_id)
  gene_names <- relist(ids2names[unlist(gene_ids)], gene_ids)
  gene_names <- sapply(gene_names, paste, collapse="; ")
  resEnrich$GenesNames <- gene_names
  return(resEnrich)
}

compileEnrichmentInput <- function(param, se){
  require(SummarizedExperiment)
  seqAnno <- data.frame(rowData(se), row.names=rownames(se),
                        check.names = FALSE, stringsAsFactors=FALSE)
  logSignal <- log2(shiftZeros(assays(se)$xNorm, param$minSignal))
  groupMeans <- cbind(rowMeans(logSignal[ , param$grouping == param$sampleGroup, 
                                          drop=FALSE]),
                      rowMeans(logSignal[ , param$grouping == param$refGroup, 
                                          drop=FALSE])
  )
  colnames(groupMeans) = c(param$sampleGroup, param$refGroup)
  normalizedAvgSignal=rowMeans(groupMeans)
  
  if (param$featureLevel != "gene"){
    genes = getGeneMapping(param, seqAnno)
    seqAnno = aggregateGoAnnotation(seqAnno, genes)
    if (!is.null(normalizedAvgSignal)){ ## if its not an identity mapping
      normalizedAvgSignal = tapply(normalizedAvgSignal[names(genes)], genes, mean)
      normalizedAvgSignal = normalizedAvgSignal[rownames(seqAnno)]
    }
    if (is.null(genes)){
      stop("no probe 2 gene mapping found found for ")
    }
  } else {
    genes = rownames(seqAnno)
    names(genes) = genes
  }
  
  isSig = rowData(se)$pValue < param$pValThreshGO & rowData(se)$usedInTest
  isUp = rowData(se)$log2Ratio > param$log2RatioThreshGO & isSig
  isDown = rowData(se)$log2Ratio < -param$log2RatioThreshGO & isSig
  probes = rownames(groupMeans)
  presentGenes = na.omit(unique(genes[probes[rowData(se)$isPresentProbe]]))
  upGenes = na.omit(unique(genes[probes[isUp]]))
  downGenes = na.omit(unique(genes[probes[isDown]]))
  bothGenes = union(upGenes, downGenes)
  normalizedAvgSignal = normalizedAvgSignal[presentGenes]
  if (length(presentGenes) == 0 | length(bothGenes) == 0){
    ezWrite("presentGenes: ", length(presentGenes), " up: ", length(upGenes), 
            " down: ", length(downGenes), " both: ", length(bothGenes))
  }
  ans <- list(upGenes=upGenes, downGenes=downGenes, bothGenes=bothGenes,
              presentGenes=presentGenes, 
              normalizedAvgSignal=normalizedAvgSignal, seqAnno=seqAnno)
  return(ans)
}

twoGroupsGOSE = function(param, enrichInput, method="Wallenius"){

  upGenes <- enrichInput$upGenes
  downGenes <- enrichInput$downGenes
  bothGenes <- enrichInput$bothGenes
  presentGenes <- enrichInput$presentGenes
  normalizedAvgSignal <- enrichInput$normalizedAvgSignal
  seqAnno = enrichInput$seqAnno
  
  job = ezJobStart("twoGroupsGO")
  require("GOstats", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  require("annotate", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)

  ontologies = c("BP", "MF", "CC")
  #goResults = list()
  #for (onto in ontologies){
  goResults = ezMclapply(ontologies, function(onto){
    gene2goList = goStringsToList(seqAnno[[paste("GO", onto)]],
                                  listNames=rownames(seqAnno))[presentGenes]
    if (param$includeGoParentAnnotation){
      gene2goList = addGoParents(gene2goList, onto) 
    }
    enrichUp = ezGoseq(param, selectedGenes=upGenes, allGenes=presentGenes, 
                       gene2goList=gene2goList, method=method, 
                       normalizedAvgSignal=normalizedAvgSignal, onto=onto)
    enrichDown = ezGoseq(param, selectedGenes=downGenes, allGenes=presentGenes, 
                         gene2goList=gene2goList, method=method, 
                         normalizedAvgSignal=normalizedAvgSignal, onto=onto)
    enrichBoth = ezGoseq(param, selectedGenes=bothGenes, allGenes=presentGenes, 
                         gene2goList=gene2goList, method=method, 
                         normalizedAvgSignal=normalizedAvgSignal, onto=onto)
    result = list(enrichUp=enrichUp, enrichDown=enrichDown, enrichBoth=enrichBoth)
    result <- lapply(result, addGeneNamesEnrich, se)
    return(result)
  }, mc.cores=1)
  names(goResults) = ontologies
  #     goResults[[onto]] = list(enrichUp=enrichUp, enrichDown=enrichDown, enrichBoth=enrichBoth,
  #                              countsUp=countsUp, countsDown=countsDown, countsBoth=countsBoth)
  ezWriteElapsed(job)
  return(goResults)
}


## ezGoseq considers only the genes that have annotations; genes without annotation are removed from the selectedGenes and allGenes
##' @describeIn twoGroupsGO Performs the GO analysis and returns a list of results.
ezGoseq = function(param, selectedGenes, allGenes, gene2goList=NULL, 
                   method=c("Wallenius", "Sampling", "Hypergeometric"),
                   onto=NULL, normalizedAvgSignal=NULL){
  method = match.arg(method)
  require("GO.db", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  #if (length(selectedGenes) <= 1){
  #  return(NA)
  #}
  
  stopifnot(names(gene2goList) %in% allGenes)
  stopifnot(!is.null(onto))
  ### consider only genes with annotation in the currently selected ontology!!!!
  allGos = switch(onto,
                  BP=keys(GOBPPARENTS),
                  MF=keys(GOMFPARENTS),
                  CC=keys(GOCCPARENTS),
                  NA)
  if (!all(unlist(gene2goList) %in% allGos)){
    gene2goList = lapply(gene2goList, function(x){intersect(x, allGos)})
  }
  
  gene2goList = gene2goList[lengths(gene2goList) > 0]
  
  ## GO analysis
  allGenes = intersect(allGenes, names(gene2goList))
  selectedGenes = intersect(selectedGenes, names(gene2goList))
  
  go2GenesList = inverseMapping(gene2goList)
  go2GenesList = go2GenesList[lengths(go2GenesList) >= param$minCountFisher]    
  goSizes = lengths(go2GenesList)
  go2SelectedGenes = lapply(go2GenesList, function(x, y){intersect(x, y)}, selectedGenes)  
  goCounts = lengths(go2SelectedGenes)
  
  gene.vector = as.integer(allGenes %in% selectedGenes)
  names(gene.vector) = allGenes
  pwf.counts = data.frame(DEgenes=gene.vector, bias.data=1, pwf=1, row.names=names(gene.vector))
  if (!is.null(normalizedAvgSignal)){
    stopifnot(names(gene.vector) %in% names(normalizedAvgSignal))
    tryCatch({pwf.counts = goseq::nullp(gene.vector, bias.data=2^normalizedAvgSignal[names(gene.vector)], plot.fit=FALSE)}, error=function(e){message("nullp failed")})
    #tryCatch({pwf.counts = goseq::nullp(gene.vector, bias.data=2^normalizedAvgSignal[names(gene.vector)], plot.fit=FALSE)})
    
  }
  #   else {
  #     stopifnot(method == "Hypergeometric") ## if there is no bias --> method must be Hypergeometric
  #     pwf.counts = data.frame(DEgenes=gene.vector, bias.data=1, pwf=1, row.names=names(gene.vector))
  #   }
  go.counts = goseq::goseq(pwf.counts, gene2cat=gene2goList, method=method)
  pvalues = go.counts[match(names(go2GenesList), go.counts$category), "over_represented_pvalue"]
  
  ## compile the result
  result = data.frame(row.names=names(go2GenesList), stringsAsFactors=FALSE, 
                      check.names=FALSE)
  result$Pvalue = pvalues
  result$fdr = p.adjust(pvalues, method="fdr")
  result$Count = goCounts
  result$Size = goSizes
  result$Term = Term(GOTERM[rownames(result)])
  result$Genes = sapply(go2SelectedGenes, function(x){paste(sort(x), collapse="; ")})
  return(result)
}


### -----------------------------------------------------------------
### ezGSEA
###
ezGSEA <- function(param, se){
  require(clusterProfiler)
  require(GO.db)
  godata <- prepareGOData(param, se)
  presentGenes <- godata$presentGenes
  seqAnno <- data.frame(rowData(se), row.names=rownames(se),
                        check.names = FALSE, stringsAsFactors=FALSE)
  geneid2name <- setNames(seqAnno$gene_name, seqAnno$gene_id)
  geneList <- setNames(rowData(se)$log2Ratio, rowData(se)$gene_id)
  geneList <- sort(geneList, decreasing = TRUE)
  
  ontologies = c("BP", "MF", "CC")
  
  goResults = ezMclapply(ontologies, function(onto){
    message("GSEA: ", onto)
    gene2goList = goStringsToList(seqAnno[[paste("GO", onto)]], 
                                  listNames=rownames(seqAnno))[presentGenes]
    if (param$includeGoParentAnnotation){
      gene2goList = addGoParents(gene2goList, onto) 
    }
    ### consider only genes with annotation in the currently selected ontology!!!!
    allGos = switch(onto,
                    BP=keys(GOBPPARENTS),
                    MF=keys(GOMFPARENTS),
                    CC=keys(GOCCPARENTS),
                    NA)
    if (!all(unlist(gene2goList) %in% allGos)){
      gene2goList = lapply(gene2goList, function(x){intersect(x, allGos)})
    }
    gene2goList = gene2goList[lengths(gene2goList) > 0]
    goIDs <- unlist(gene2goList, use.names=FALSE)
    go2geneDF <- data.frame(ont=goIDs,
                            gene=rep(names(gene2goList), lengths(gene2goList)),
                            stringsAsFactors = FALSE)
    resGSEA <- GSEA(gene=geneList, TERM2GENE=go2geneDF,
                    by="fgsea")
    tempTable <- resGSEA@result
    if(nrow(tempTable) != 0L){
      tempTable$Description <- substr(Term(GOTERM[tempTable$ID]), 1, 30)
      tempTable$geneName <- sapply(relist(geneid2name[unlist(strsplit(tempTable$core_enrichment, "/"))], strsplit(tempTable$core_enrichment, "/")), paste, collapse="/")
      resGSEA@result <- tempTable
    }
    return(resGSEA)
  }, mc.cores=1)
  names(goResults) = ontologies
  return(goResults)
}

##' @title Groups GO terms and information
##' @description Groups GO terms and information.
##' @param selectedGenes a character vector containing the selected genes.
##' @param go2GeneList a character vector.
##' @param onto a character representing the ontology to use.
##' @param levels an integer vector selecting the levels to get GO information by.
##' @param goSlim a character vector.
##' @template roxygen-template
##' @return Returns the results of grouping the GO information.
ezGroupGO = function(selectedGenes, go2GeneList, onto="CC", levels = 2:4, goSlim=NULL) {
  goByLevel = ezGetGoByLevels(onto, levels, goSlim=goSlim)
  levelCounts = list()
  for (l in levels){
    goIds = goByLevel[[paste("level", l)]]
    go2selGenes = lapply(go2GeneList[goIds], function(genesInCategory) intersect(selectedGenes, genesInCategory))
    result = data.frame("GO ID" = goIds, Level=rep(paste("level", l), length(goIds)), stringsAsFactors=FALSE, check.names=FALSE)
    result[["Term"]] = Term(GOTERM[goIds])
    result[["Count"]] = sapply(go2selGenes, length)
    result[["Size"]] = sapply(go2GeneList[goIds], length)
    result[["Genes"]] = sapply(go2selGenes, paste, collapse="; ")
    levelCounts[[paste("level", l)]] =result
  }
  levelCounts[["unknown"]] = data.frame("GO ID"="unknown", Level="unknown", Term="unknown",
                                        Count=length(setdiff(selectedGenes, unlist(go2GeneList))),
                                        Size=NA, Genes="", stringsAsFactors=FALSE, check.names=FALSE
  )
  result = do.call("rbind", levelCounts)
  return(result)
}

## copied from clusterProfiler
##' @describeIn ezGroupGO Gets the gene ontology by level.
ezGetGoByLevels = function(onto, levels, goSlim=NULL) {
  require("GO.db", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  switch(onto, MF = {
    topNode <- "GO:0003674"
    children <- GOMFCHILDREN
  }, BP = {
    topNode <- "GO:0008150"
    children <- GOBPCHILDREN
  }, CC = {
    topNode <- "GO:0005575"
    children <- GOCCCHILDREN
  })
  if (!is.null(goSlim)){
    goAncestor = switch(onto, BP=as.list(GOBPANCESTOR),
                        CC=as.list(GOCCANCESTOR),
                        MF=as.list(GOMFANCESTOR))
    goSlim = goSlim[goSlim %in% names(goAncestor)]
    goSlim = setdiff(c(goSlim, unlist(goAncestor[goSlim])), "all")
  }
  
  nodes = topNode
  nodesByLevel = list("level 1"=nodes)
  for (i in seq_len(max(levels) - 1)) {
    nodes <- keys(children[nodes])
    nodes <- na.omit(unique(unlist(nodes)))
    goIds = as.character(nodes)
    if (!is.null(goSlim)){
      goIds = intersect(goIds, goSlim)
    }
    nodesByLevel[[paste("level", i+1)]] = goIds
  }
  return(nodesByLevel[paste("level", levels)])
}

##' @title Parses GO strings to a list
##' @description Parses GO strings to a list.
##' @param goStrings a character vector containing the GO strings
##' @param listNames an optional character vector naming the list elements.
##' @template roxygen-template
##' @return Returns the GO strings as a list.
##' @examples
##' goStrings = c("GO 1;GO 2;GO 3", "GO 4")
##' listNames = letters[1:2]
##' goStringsToList(goStrings, listNames)
goStringsToList = function(goStrings, listNames=NULL){
  x = strsplit(goStrings, "; ")
  names(x) = listNames
  #  x = lapply(x, function(x){ifelse(x == "", character(0), x)})
  ## "" after strsplit is character(0)
  return(x)
}

##' @describeIn goClusterTable Gets the child terms from a GO term and returns them together.
getChildTerms = function(x, subset, goRelatives, indent="", childEnvir){
  
  result = character()
  if (is.element(x, subset)){
    term = getGOTerm(x)[[1]]
    displayLabel = paste0(indent, term)
    result[displayLabel] = x
    subset = setdiff(subset, x) ## this is the modification to make the table non-redundant!!
  }
  kids = intersect(childEnvir[[x]], goRelatives)
  indent = paste(indent, ".")
  for (kid in kids){
    kidResult = getChildTerms(kid, subset, goRelatives, indent=indent, childEnvir)
    subset = setdiff(subset, kidResult) ## again the non-redundancy modification
    result = c(result, kidResult)
  }
  return(result)
}

##' @describeIn clusterHeatmap Initializes and returns a list of results used by \code{clusterHeatmap()} and \code{goClusterResults()}.
clusterResults = function(x, nClusters=5, clusterColors=rainbow(nClusters), 
                          d=NULL, method="ward.D2"){
  if (is.null(d)){
    xx = x
    xx[is.na(x)] = min(x, na.rm=TRUE)
    d = dist(xx)
  }
  hcl = hclust(d, method=method)
  clusterNumbers = cutree(hcl, k=nClusters)
  result = list()
  result$nClusters = nClusters
  result$clusterNumbers = clusterNumbers
  result$clusterColors = clusterColors
  result$hcl = hcl
  return(result)
}

##' @title Plots the cluster heatmap
##' @description Plots the cluster heatmap.
##' @param x the data to plot.
##' @param param a list of parameters to extract the logical \code{showGeneClusterLabels} from.
##' @param result a list of cluster result properties obtained with \code{clusterResults()}.
##' @param file a character representing the name of the .png file to derive tables from.
##' @param method a character representing the method to pass to \code{hclust()}.
##' @param doClusterColumns a logical indicating whether to do cluster columns.
##' @param columnDist the distance matrix to use for the cluster columns.
##' @param colColors used for the \code{ColSideColors} argument passed to \code{heatmap.2()}.
##' @param lim two integers used for \code{breaks} argument passed to \code{heatmap.2()}.
##' @param cexRow an integer passed to \code{heatmap.2()}.
##' @param cexCol an integer passed to \code{heatmap.2()}.
##' @param labRow a character vector, possibly modified, then passed to \code{heatmap.2()}.
##' @param margins an integer, possibly modified, then passed to \code{heatmap.2()}.
##' @template colors-template
##' @param maxGenesWithLabel an integer specifying the maximum amount of genes with labels.
##' @template addargs-template
##' @templateVar fun heatmap.2()
##' @template roxygen-template
##' @seealso \code{\link[stats]{hclust}}
##' @seealso \code{\link[gplots]{heatmap.2}}
clusterHeatmap = function(x, param, result, file="cluster-heatmap.png", 
                          method="ward.D2",
                          doClusterColumns=FALSE, columnDist=NULL, 
                          colColors=NULL, lim=c(-4, 4),
                          lwid=c(1, 4), lhei=c(1,5),
                          cexRow=1.0, cexCol=1.5, labRow=rownames(x), 
                          margins=c(14,9),
                          colors=getBlueRedScale(), 
                          maxGenesWithLabel=50, ...){
  require("gplots", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  probeDendro = as.dendrogram(result$hcl)
  probeDendro = reorder(probeDendro, rowMeans(x, na.rm=TRUE))
  if (doClusterColumns){
    if (is.null(columnDist)){
      columnDist = dist(t(x))
    }
    result$hcl = hclust(columnDist, method=method)
    colDendro = as.dendrogram(result$hcl)
    colDendro = reorder(colDendro, colMeans(x, na.rm=TRUE))
    showDendro = "both"
  } else {
    colDendro = FALSE
    showDendro = "row"
  }
  if (param$showGeneClusterLabels & nrow(x) < maxGenesWithLabel){
    if (is.null(labRow)){
      labRow = rownames(x)
    }
  } else {
    labRow = ""
    margins[2] = 0.5 
  }
  if (!is.null(colColors)){
    if (result$nClusters > 1){
      heatmap.2(x,
                ColSideColors=colColors, RowSideColors=result$clusterColors[result$clusterNumbers],
                scale="none", dendrogram=showDendro, labRow=labRow,
                breaks=seq(from=lim[1], to=lim[2], length.out=257),
                Colv=colDendro, Rowv=probeDendro, col=colors,
                key=TRUE, density.info="none", trace="none", margins=margins, cexCol=cexCol, cexRow=cexRow, lwid=lwid, lhei=lhei, ...)
    } else {
      heatmap.2(x,
                ColSideColors=colColors,
                scale="none", dendrogram=showDendro, labRow=labRow,
                breaks=seq(from=lim[1], to=lim[2], length.out=257),
                Colv=colDendro, Rowv=probeDendro, col=colors,
                key=TRUE, density.info="none", trace="none", margins=margins, cexCol=cexCol, cexRow=cexRow, lwid=lwid, lhei=lhei, ...)
    }
  } else {
    if (result$nClusters > 1){
      heatmap.2(x,
                RowSideColors=result$clusterColors[result$clusterNumbers],
                scale="none", dendrogram=showDendro, labRow=labRow,
                breaks=seq(from=lim[1], to=lim[2], length.out=257),
                Colv=colDendro, Rowv=probeDendro, col=colors,
                key=TRUE, density.info="none", trace="none", margins=margins, cexCol=cexCol, cexRow=cexRow, lwid=lwid, lhei=lhei,...)
    } else {
      heatmap.2(x,
                scale="none", dendrogram=showDendro, labRow=labRow,
                breaks=seq(from=lim[1], to=lim[2], length.out=257),
                Colv=colDendro, Rowv=probeDendro, col=colors,
                key=TRUE, density.info="none", trace="none", margins=margins, cexCol=cexCol, cexRow=cexRow, lwid=lwid, lhei=lhei, ...)
    }
  }
  if (!is.null(file)){
    ## export also the clusters
    xTmp = as.data.frame(result$clusterNumbers)
    colnames(xTmp) = "Cluster"
    xTmp[ ,1] = result$clusterColors[xTmp[,1]]
    xTmp = xTmp[order(xTmp$Cluster), ,drop=FALSE]
    ezWrite.table(xTmp, file=sub(".png$", "-clusterMembers.txt", file))
  }
}

clusterPheatmap <- function(x, design, param, 
                            clusterColors=c("red", "yellow", "orange", 
                                            "green", "blue", "cyan"), 
                            method="ward.D2", doClusterColumns=FALSE,
                            colors=getBlueRedScale(),
                            colColors=NULL, lim=c(-4, 4),
                            maxGenesWithLabel=50,
                            sampleColors=NULL){
  require(pheatmap)
  nClusters <- length(clusterColors)
  
  if(param$showGeneClusterLabels & nrow(x) < maxGenesWithLabel){
    isShowRowNames <- TRUE
  }else{
    isShowRowNames <- FALSE
  }
  
  callback = function(hc, mat){
    dend = reorder(as.dendrogram(hc), wts = rowMeans(mat, na.rm=TRUE))
    as.hclust(dend)
  }
  clusterInfo <- pheatmap(x, color=colors, clustering_method=method,
                          breaks=seq(from=lim[1], to=lim[2], length.out=257),
                          scale="none",
                          clustering_callback = callback,
                          silent=TRUE)
  
  clusters <- as.factor(cutree(clusterInfo$tree_row, nClusters))
  annotation_row = data.frame(Clusters=clusters)
  if(doClusterColumns){
    colDendro <- clusterInfo$tree_col
  }else{
    colDendro <- FALSE
  }
  
  # define annotation_colors list to avoid default pheatmap colors.
  ann_colors <- list(setNames(clusterColors, levels(clusters)),
                     setNames(unique(sampleColors), unique(design[[1]])))
  names(ann_colors) = c(colnames(annotation_row),colnames(design[1]))
  
  p <- pheatmap(x, color=colors, clustering_method=method,
           breaks=seq(from=lim[1], to=lim[2], length.out=257),
           scale="none", cluster_rows=clusterInfo$tree_row,
           cluster_cols=colDendro,
           show_rownames=isShowRowNames,
           annotation_col = design, annotation_row=annotation_row,
           annotation_colors = ann_colors)
  
  ans <- list(nClusters=nClusters, clusterNumbers=clusters,
              clusterColors=clusterColors, hcl=clusterInfo$tree_row,
              pheatmap=p)
  invisible(ans)
}

##' @describeIn clusterHeatmap Applies a GO analysis to the cluster results if GO should be done.
goClusterResults = function(x, param, result, ontologies=c("BP", "MF", "CC"), seqAnno=NULL,
                            universeGeneIds=NULL, universeProbeIds=NULL, keggOrganism=NA){
  require("GOstats", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  require("annotate", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)

  if (is.null(universeGeneIds)){
    universeGeneIds =na.omit(unique(unlist( genes[universeProbeIds]) ))
  }

  genesByCluster = tapply(genes[rownames(x)], result$clusterNumbers, function(x){na.omit(unique(x))}, simplify=FALSE)
  goClusterResults = ezMclapply(ontologies, function(onto){
    gene2go = goStringsToList(seqAnno[[paste("GO", onto)]], listNames=rownames(seqAnno))[universeGeneIds]
    if (param$includeGoParentAnnotation){
      gene2go = addGoParents(gene2go, onto)
    }
    result = lapply(genesByCluster, function(selectedGenes){
      ezGoseq(param, selectedGenes, universeGeneIds, gene2goList=gene2go,
              method="Hypergeometric",
              onto=onto, normalizedAvgSignal=NULL)})
  }, mc.cores=1)
  names(goClusterResults) = ontologies
  result$GO = goClusterResults
  return(invisible(result))
}



## still used?
myFisherTest = function(goGenes, selGenes, allGenes, alternative="greater"){ 
  fisher.test(factor(allGenes %in% selGenes, levels=c("FALSE", "TRUE")), 
              factor(allGenes %in% goGenes, levels=c("FALSE", "TRUE")), alternative=alternative)$p.value
}

## NOTE: the GSEABase has also GO to slim mapping but does not include the GO parents
## at the geneontology there seems to be a map2slim perl script but
## a) the perl script installation seems cumbersome; lots of dependencies
## b) it requires a gene-assocation-file as input; a more complex than necessary file
mapGoToSlim = function(goList, ontology, slimGo){
  require("GO.db", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  #   require(GSEABase)
  #   slim = getOBOCollection(oboFile)
  #   slimGo = ids(slim)
  goAncestor = switch(ontology, BP=as.list(GOBPANCESTOR),
                      CC=as.list(GOCCANCESTOR),
                      MF=as.list(GOMFANCESTOR))
  goOffspring = switch(ontology, MF=as.list(GOMFOFFSPRING),
                       BP=as.list(GOBPOFFSPRING), CC=as.list(GOCCOFFSPRING))
  slimGo = slimGo[slimGo %in% names(goAncestor)]
  slimGo = setdiff(union(slimGo, unlist(goAncestor[slimGo])), "all") ## this automatically filters on the ontology  
  stopifnot(slimGo %in% names(goOffspring))
  slim2full = sapply(slimGo, function(sg){c(sg, goOffspring[[sg]])}, simplify=FALSE)
  full2slim = inverseMapping(slim2full)
  ## remove all parents
  full2slim = lapply(full2slim, function(slimIds){setdiff(slimIds, unlist(goAncestor[slimIds]))})
  slimList = lapply(goList, function(goIds){unique(unlist(full2slim[goIds]))})
  return(slimList)
}



