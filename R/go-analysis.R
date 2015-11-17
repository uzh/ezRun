###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Do GO?
##' @description Decides whether to do a gene ontologies analysis.
##' @param param a list of parameters to extract the logical \code{runGO} and the character \code{featureLevel} from.
##' @param seqAnno the sequence annotation.
##' @template roxygen-template
##' @return Returns a logical.
##' @examples
##' param = ezParam()
##' param$ezRef@@refFeatureFile = "./inst/extdata/genes.gtf"
##' param$ezRef@@refAnnotationFile = "./inst/extdata/genes_annotation_example.txt"
##' param$ezRef@@refFastaFile = "/srv/GT/reference/Saccharomyces_cerevisiae/Ensembl/EF4/Sequence/WholeGenomeFasta/genome.fa"
##' seqAnno = writeAnnotationFromGtf(param)
##' doGo(param, seqAnno)
doGo = function(param, seqAnno){
  flag = param$runGO && (hasGeneMapping(param, seqAnno) || param$featureLevel == "gene") && hasGoAnnotation(seqAnno)
  message("doGo: ", flag)
  return(flag)
}

##' @describeIn doGo Checks if the Annotation contains gene ontologies.
hasGoAnnotation = function(seqAnnoDF){
  any(c("GO BP", "GO MF", "GO CC") %in% colnames(seqAnnoDF))
}

##' @title Separates GO ID's by ontology
##' @description Separates GO ID's by ontology
##' @param goIdStrings the GO ID's to separate.
##' @template roxygen-template
##' @return Returns the separated GO ID's
## TODOEXAMPLE: example
separateGoIdsByOnto = function(goIdStrings){
  library(GO.db, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  bpGos = keys(GOBPPARENTS)
  mfGos = keys(GOMFPARENTS)
  ccGos = keys(GOCCPARENTS)
  result = data.frame(row.names=1:length(goIdStrings))
  goList = strsplit(trimWhiteSpace(gsub("; ", ";", goIdStrings)), ";") ## support both separators "; " and ";" 
  result[["GO BP"]] = sapply(lapply(goList, function(x,y){intersect(x,y)}, y=bpGos), collapse)
  result[["GO MF"]] = sapply(lapply(goList, function(x,y){intersect(x,y)}, y=mfGos), collapse)
  result[["GO CC"]] = sapply(lapply(goList, function(x,y){intersect(x,y)}, y=ccGos), collapse)
  return(result)
}

##' @title Gets the GO parents
##' @description Gets the GO parents.
##' @param id the ID to check
##' @param onto the ontology to use
##' @template roxygen-template
##' @return Returns the ancestor terms if they are specified.
getGOparents = function(id, onto="BP"){
  go2Ancestor = switch(onto, BP=as.list(GOBPANCESTOR), MF=as.list(GOMFANCESTOR), CC=as.list(GOCCANCESTOR))
  ancestorIds = setdiff(go2Ancestor[[id]], "all") ## remove the all category
  if (is.null(ancestorIds)){
    return()
  } else {
    ancestorTerms = Term(ancestorIds[GOTERM])
    return(ancestorTerms)
  }
}

## seems very similar to getGOparents
addGoParents = function(gene2goList, onto){
  goParents = switch(onto, BP=as.list(GOBPPARENTS),
                     CC=as.list(GOCCPARENTS),
                     MF=as.list(GOMFPARENTS))
  return(lapply(gene2goList, function(x){setdiff(union(x, unlist(goParents[x])), "all")}))
}

##' @title Performs the GO analysis for two groups
##' @description Performs the GO analysis for two groups.
##' @param param a list of parameters:
##' \itemize{
##'   \item{featureLevel}{ which feature level to use.}
##'   \item{pValThreshGO}{ a numeric specifying the threshold for the GO p-Value.}
##'   \item{log2RatioThreshGO}{ a numeric specifying the threshold for the GO log2 ratios.}
##'   \item{includeGoParentAnnotation}{ a logical indicating whether to include the annotation of the GO parents.}
##' }
##' @param testResult a list containing the results of an earlier test.
##' @param seqAnno the sequence annotation.
##' @param normalizedAvgSignal an optional normalized average signal.
##' @param method which method to pass to \code{ezGoseq()}.
##' @template roxygen-template
##' @return Returns the results of the GO analysis.
##' @seealso \code{\link{ezGoseq}}
twoGroupsGO = function(param, testResult, seqAnno, normalizedAvgSignal=NULL, method="Wallenius"){
  
  job = ezJobStart("twoGroupsGO")
  library(GOstats, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  library(annotate, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  
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
  
  isSig = testResult$pValue < param$pValThreshGO & testResult$usedInTest
  isUp = testResult$log2Ratio > param$log2RatioThreshGO & isSig
  isDown = testResult$log2Ratio < -param$log2RatioThreshGO & isSig
  probes = rownames(testResult$groupMeans)
  presentGenes = na.omit(unique(genes[probes[testResult$isPresentProbe]]))
  upGenes = na.omit(unique(genes[probes[isUp]]))
  downGenes = na.omit(unique(genes[probes[isDown]]))
  bothGenes = union(upGenes, downGenes)
  normalizedAvgSignal = normalizedAvgSignal[presentGenes]
  if (length(presentGenes) == 0 | length(bothGenes) == 0){
    ezWrite("presentGenes: ", length(presentGenes), " up: ", length(upGenes), " down: ", length(downGenes), " both: ", length(bothGenes))
  }
  ontologies = c("BP", "MF", "CC")
  #goResults = list()
  #for (onto in ontologies){
  goResults = ezMclapply(ontologies, function(onto){
    gene2goList = goStringsToList(seqAnno[[paste("GO", onto)]], listNames=rownames(seqAnno))[presentGenes]
    if (param$includeGoParentAnnotation){
      gene2goList = addGoParents(gene2goList, onto) 
    }
    enrichUp = ezGoseq(param, selectedGenes=upGenes, allGenes=presentGenes, gene2goList=gene2goList, method=method, normalizedAvgSignal=normalizedAvgSignal, onto=onto)
    enrichDown = ezGoseq(param, selectedGenes=downGenes, allGenes=presentGenes, gene2goList=gene2goList, method=method, normalizedAvgSignal=normalizedAvgSignal, onto=onto)
    enrichBoth = ezGoseq(param, selectedGenes=bothGenes, allGenes=presentGenes, gene2goList=gene2goList, method=method, normalizedAvgSignal=normalizedAvgSignal, onto=onto)
    result = list(enrichUp=enrichUp, enrichDown=enrichDown, enrichBoth=enrichBoth)
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
  library(GO.db, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  library(goseq, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  if (length(selectedGenes) <= 1){
    return(NA)
  }
  
  stopifnot(names(gene2goList) %in% allGenes)
  stopifnot(!is.null(onto))
  ### consider only genes with annotation in the currently selected ontology!!!!
  allGos = switch(onto,
                  BP=keys(GOBPPARENTS),
                  MF=keys(GOMFPARENTS),
                  CC=keys(GOCCPARENTS),
                  NA)
  if (!all(unlist(gene2goList) %in% allGos)){
    gene2goList = lapply(gene2goList, function(x){ intersect(x, allGos)})
  }
  
  gene2goList = gene2goList[sapply(gene2goList, length) > 0]
  
  ## GO analysis
  allGenes = intersect(allGenes, names(gene2goList))
  selectedGenes = intersect(selectedGenes, names(gene2goList))
  
  go2GenesList= inverseMapping(gene2goList)
  go2GenesList = go2GenesList[sapply(go2GenesList, length) >= param$minCountFisher]    
  goSizes = sapply(go2GenesList, length)
  go2SelectedGenes = lapply(go2GenesList, function(x, y){ intersect(x, y)}, selectedGenes)  
  goCounts = sapply(go2SelectedGenes, length)
  
  gene.vector = as.integer(allGenes %in% selectedGenes)
  names(gene.vector) = allGenes
  pwf.counts = data.frame(DEgenes=gene.vector, bias.data=1, pwf=1, row.names=names(gene.vector))
  if (!is.null(normalizedAvgSignal)){
    stopifnot(names(gene.vector) %in% names(normalizedAvgSignal))
    tryCatch({pwf.counts = nullp(gene.vector, bias.data=2^normalizedAvgSignal[names(gene.vector)], plot.fit=FALSE)}, error=function(e){message("nullp failed")})
    #tryCatch({pwf.counts = nullp(gene.vector, bias.data=2^normalizedAvgSignal[names(gene.vector)], plot.fit=FALSE)})
    
  }
  #   else {
  #     stopifnot(method == "Hypergeometric") ## if there is no bias --> method must be Hypergeometric
  #     pwf.counts = data.frame(DEgenes=gene.vector, bias.data=1, pwf=1, row.names=names(gene.vector))
  #   }
  go.counts = goseq(pwf.counts, gene2cat=gene2goList, method=method)
  pvalues = go.counts[match(names(go2GenesList), go.counts$category), "over_represented_pvalue"]
  
  ## compile the result
  result = data.frame(row.names=names(go2GenesList),stringsAsFactors=FALSE, check.names=FALSE)
  result$Pvalue = pvalues
  result$fdr = p.adjust(pvalues, method="fdr")
  result$Count = goCounts
  result$Size = goSizes
  result$Term = Term(GOTERM[rownames(result)])
  result$Genes = sapply(go2SelectedGenes, function(x){paste(sort(x), collapse="; ")})
  return(result)
}

##' @title 1
##' @description 1
##' @param selectedGenes a character vector containing the selected genes.
##' @param go2GeneList 1
##' @param onto a character representing the ontology to use.
##' @param levels 1
##' @param goSlim 1
##' @template roxygen-template
##' @return Returns the results of the GO analysis.
## TODOP: finish documenting
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
ezGetGoByLevels = function(ont, levels, goSlim=NULL) {
  library(GO.db, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  switch(ont, MF = {
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
    goAncestor = switch(ont, BP=as.list(GOBPANCESTOR),
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
##' listNames = letters[1:4]
##' goStringsToList(goStrings, listNames)
goStringsToList = function(goStrings, listNames=NULL){
  x = strsplit(goStrings, "; ")
  names(x) = listNames
  x = lapply(x, function(x){ifelse(x == "", character(0), x)})
  return(x)
}

##' @describeIn goClusterTable Gets the child terms from a GO term and returns them together.
getChildTerms = function(x, subset, goRelatives, indent="", childEnvir){
  
  result = character()
  if (is.element(x, subset)){
    term = getGOTerm(x)[[1]]
    displayLabel = paste0(indent, substr(term, 1, 30))
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

##' @describeIn clusterHeatmap Performs some initializations and returns a list of results used by \code{clusterHeatmap()} and \code{goClusterResults()}.
clusterResults = function(x, nClusters=5, clusterColors=rainbow(nClusters), d=NULL, method="ward.D2"){
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
##' @param x a data vector to plot.
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
##' @templateVar fun heatmap.2
##' @template roxygen-template
##' @seealso \code{\link[stats]{hclust}}
##' @seealso \code{\link[gplot]{heatmap.2}}
clusterHeatmap = function(x, param, result, file="cluster-heatmap.png", method="ward.D2",
                          doClusterColumns=FALSE, columnDist=NULL, colColors=NULL, lim=c(-4, 4),
                          cexRow=1.0, cexCol=1.5, labRow=rownames(x), margins=c(14,9),
                          colors=getBlueRedScale(256), maxGenesWithLabel=50, ...){
  require(gplots, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
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
                key=TRUE, density.info="none", trace="none", margins=margins, cexCol=cexCol, cexRow=cexRow, ...)
    } else {
      heatmap.2(x,
                ColSideColors=colColors,
                scale="none", dendrogram=showDendro, labRow=labRow,
                breaks=seq(from=lim[1], to=lim[2], length.out=257),
                Colv=colDendro, Rowv=probeDendro, col=colors,
                key=TRUE, density.info="none", trace="none", margins=margins, cexCol=cexCol, cexRow=cexRow, ...)
    }
  } else {
    if (result$nClusters > 1){
      heatmap.2(x,
                RowSideColors=result$clusterColors[result$clusterNumbers],
                scale="none", dendrogram=showDendro, labRow=labRow,
                breaks=seq(from=lim[1], to=lim[2], length.out=257),
                Colv=colDendro, Rowv=probeDendro, col=colors,
                key=TRUE, density.info="none", trace="none", margins=margins, cexCol=cexCol, cexRow=cexRow, ...)
    } else {
      heatmap.2(x,
                scale="none", dendrogram=showDendro, labRow=labRow,
                breaks=seq(from=lim[1], to=lim[2], length.out=257),
                Colv=colDendro, Rowv=probeDendro, col=colors,
                key=TRUE, density.info="none", trace="none", margins=margins, cexCol=cexCol, cexRow=cexRow, ...)
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

##' @describeIn clusterHeatmap Applies a GO analysis to the cluster results if GO should be done.
goClusterResults = function(x, param, result, ontologies=c("BP", "MF", "CC"), seqAnno=NULL,
                            universeGeneIds=NULL, universeProbeIds=NULL, keggOrganism=NA){
  require(GOstats, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  require(annotate, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  if (param$featureLevel != "gene"){
    genes = getGeneMapping(param, seqAnno)
    if (is.null(genes)){
      stop("no probe 2 gene mapping found found")
    }
    seqAnno = aggregateGoAnnotation(seqAnno, genes)
  } else {
    genes = rownames(seqAnno)
    names(genes) = genes
  }
  
  if (is.null(universeGeneIds)){
    universeGeneIds =na.omit(unique(unlist( genes[universeProbeIds]) ))
  }
  keggClusterResults = NULL
  if (!is.na(keggOrganism)){
    kegg2GenesList = getPathway2GenesList(param, keggOrganism)
    if (all(unlist(kegg2GenesList)) %in% universeGeneIds){
      keggClusterResults = list()
    }
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
collectChildren = function(go2geneList, onto){
  goOffsprings = switch(onto, BP=as.list(GOBPOFFSPRING),
                        CC=as.list(GOCCOFFSPRING),
                        MF=as.list(GOMFOFFSPRING))
  goNames = union(names(go2geneList), names(goOffsprings))
  names(goNames) = goNames
  doCollect = function(goName, goOffsprings=NULL, go2genes=NULL){
    return(unique(unlist(go2genes[c(goOffsprings[[goName]], goName)])))
  }
  go2geneList = lapply(goNames, doCollect, goOffsprings=goOffsprings, go2genes=go2geneList)
  return(go2geneList)
}


## still used?
myFisherTest = function(goGenes, selGenes, allGenes, alternative="greater"){ 
  fisher.test(factor(allGenes %in% selGenes, levels=c("FALSE", "TRUE")), 
              factor(allGenes %in% goGenes, levels=c("FALSE", "TRUE")), alternative=alternative)$p.value
}


writeGOClusterResult = function(html, param, clusterResult){
  
  ontologies = names(clusterResult$GO)
  tables = ezMatrix("", rows=paste("Cluster", 1:clusterResult$nClusters), cols=ontologies)
  for (onto in ontologies){
    for (i in 1:clusterResult$nClusters){
      x = clusterResult$GO[[onto]][[i]]
      tables[i, onto] = goResultToHtmlTable(x, param$pValThreshFisher, param$minCountFisher, onto=onto);
    }
  }
  writeTableToHtml(tables, con=html, border=2,
                   bgcolors=matrix(gsub("FF$", "", clusterResult$clusterColors), nrow=clusterResult$nClusters, ncol=1))
}


# does not work anymore
writeGOTables = function(html , param, goResult){
  
  ezWrite("<h3>GO Enrichment Analysis</h3>", con=html)
  ezWrite("<p>Red GO categories are overrepresented among the significantly upregulated genes<br>", con=html)
  ezWrite("<p>Blue GO categories are overrepresented among the significantly downregulated genes<br>", con=html)
  ezWrite("<p>Black GO categories are overrepresented among all signifcantly regulated genes</p>", con=html)
  ezWrite("<p>Maximum number of terms displayed: ", param$maxNumberGroupsDisplayed, "</p>", con=html)
  tables = ezMatrix("", rows="Cats", cols=names(goResult))
  txtFiles = character()
  for (onto in names(goResult)){
    x = goResult[[onto]]
    tables[1,onto] = goResultToHtmlTable2(x, param$pValThreshFisher, 
                                          param$minCountFisher, onto=onto, maxNumberOfTerms=param$maxNumberGroupsDisplayed)  ### function deprecated
    for (sub in names(x)){ #c("enrichUp", "enrichDown", "enrichBoth")){
      xSub = x[[sub]]
      if (is.data.frame(xSub)){
        name = paste0(onto, "-", param$comparison, "-", sub)
        if (!is.null(xSub$Pvalue)){
          xSub = xSub[order(xSub$Pvalue), ]
          xSub = cbind("GO ID"=rownames(xSub), xSub)
        }
        txtFile = ezValidFilename(paste0(name, ".txt"), replace="-")
        ezWrite.table(xSub, file=txtFile, row.names=FALSE)
        if (param$doZip){
          txtFiles[name] = zipFile(txtFile)
        } else {
          txtFiles[name] = txtFile
        }
      }
    }
  }
  writeTableToHtml(tables, con=html, border=2, valign="top")
  if (param$doZip){
    writeTxtLinksToHtml(txtFiles, con=html, mime="application/zip")
  } else {
    writeTxtLinksToHtml(txtFiles, con=html, mime="application/txt")  
  }
}


## NOTE: the GSEABase has also GO to slim mapping but does not include the GO parents
## at the geneontology there seems to be a map2slim perl script but
## a) the perl script installation seems cumbersome; lots of dependencies
## b) it requires a gene-assocation-file as input; a more complex than necessary file
mapGoToSlim = function(goList, ontology, slimGo){
  library(GO.db, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
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



#
#writeChildren = function(x, subset, goRelatives, indent=""){
#
#  if (length(intersect(x, subset)>0)){
#    ezWrite(indent, x, " ", getGOTerm(x)[[1]])
#  }
#  kids = get(x, envir=eval(paste0("GO", onto, "CHILDREN")))
#  indent = paste(indent, " ")
#  for (kid in kids){
#    if (is.element(kid, goRelatives)){
#      writeChildren(kid, subset, goRelatives, indent=indent)
#    }
#  }
#}
#


### deprecated should no longer be used
# makeChipAnno = function(param, probeAnno, chip){
#   
#   if (missing(chip)){
#     stop("wrong call to makeChipAnno")
#   }
#   if (hasGeneMapping(param, probeAnno)){
#     message("make chip anno for ", chip, " -- have gene mapping")
#   } else {
#     message("make chip anno for ", chip, " -- gene mapping is missing in: ", paste(colnames(probeAnno), collapse=", "))
#     return(NULL)
#   }
#   
#   require(GOstats, quietly=TRUE)
#   require(GO.db, quietly=TRUE)
#   require(annotate, quietly=TRUE)
#   
#   #job = ezJobStart("entrez")
#   entrezid = new.env(hash=TRUE, parent=emptyenv(), size=nrow(probeAnno))
#   #  apply(probeAnno, 1, function(x){  assign(x["Probe"], unname(x["Gene Name"]), envir=entrezid)})
#   genes = getGeneMapping(param, probeAnno)
#   genes[genes == ""] = NA
#   probe2Gene = as.list(genes)
#   l2e(probe2Gene, entrezid)
#   #  apply(probeAnno, 1, function(x){  assign(x["Probe"], unname(x[ "Entrez Gene ID"]), envir=entrezid)})
#   #ezWriteElapsed(job)
#   
#   #job = ezJobStart("go")
#   if (all(c("GO BP", "GO CC", "GO MF") %in% colnames(probeAnno))){
#     go = new.env(hash=TRUE, parent=emptyenv(), size=nrow(probeAnno))
#     goIds = paste(probeAnno[ , "GO BP"], probeAnno[ ,"GO CC"], probeAnno[ ,"GO MF"], sep="; ")
#     goIds = sub("; ; ", "; ", goIds)
#     goIds = sub("^; ", "", goIds)
#     goIds = sub("; $", "", goIds)
#     probe2GoList = strsplit(goIds, "; ", fixed=TRUE)
#     
#     names(probe2GoList) = rownames(probeAnno)
#     ##goList = lapply(goList, function(x){ if (length(x)== 0){NA} else {paste0("GO:", x)}})
#     l2e(probe2GoList, go)
#     #ezWriteElapsed(job)
#     
#     #job = ezJobStart("go prep")
#     
#     #  goTermVec = unlist(probe2GoList)
#     #  probeVec = rep(rownames(probeAnno), sapply(probe2GoList, length))
#     #  probesByGO = split(probeVec, goTermVec)
#     probesByGo = inverseMapping(probe2GoList)
#     #ezWriteElapsed(job)
#     
#     #job = ezJobStart("go MF")
#     
#     allGo = names(probesByGo)
#     collectProbes = function(goName, goOffsprings, probesByGo){
#       return(unique(unlist(probesByGo[c(goName, goOffsprings)])))
#     }
#     
#     
#     goBpOff = as.list(GOBPOFFSPRING)
#     bpList = mapply(collectProbes, names(goBpOff), goBpOff, MoreArgs=list(probesByGo=probesByGo))
#     
#     goMfOff = as.list(GOMFOFFSPRING)
#     mfList = mapply(collectProbes, names(goMfOff), goMfOff, MoreArgs=list(probesByGo=probesByGo))
#     
#     goCcOff = as.list(GOCCOFFSPRING)
#     ccList = mapply(collectProbes, names(goCcOff), goCcOff, MoreArgs=list(probesByGo=probesByGo))
#     
#     goList = c(bpList, mfList, ccList)
#     stopifnot(!is.na(unlist(goList)))
#     
#     go2probe = new.env(hash=TRUE, parent=emptyenv(), size=length(goList))
#     l2e(goList, go2probe)
#     
#     l = list()
#     l[[paste0(chip, "GO")]] = go
#     l[[paste0(chip, "GO2ALLPROBES")]] = go2probe
#   } else {
#     l = list()
#   }
#   l[[paste0(chip, "ENTREZID")]] = entrezid
#   pkgName = paste0("package:", chip, ".db")
#   ezWrite("attach package: ", pkgName)
#   attach(l, name=pkgName)
#   return(chip)
# }
# 

# hasGoAnnotationInEnv = function(chip){
#   exists(paste0(chip, "GO2ALLPROBES"))
# }
