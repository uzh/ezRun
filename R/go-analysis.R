###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


# getGoCounts = function(genes, seqAnno)


doGo = function(param, seqAnno){
  flag = param$runGO && (hasGeneMapping(param, seqAnno) || param$featureLevel == "gene") && hasGoAnnotation(param, seqAnno)
  message("doGo: ", flag)
  return(flag)
}


hasGoAnnotation = function(param, seqAnnoDF){
  any(c("GO BP", "GO MF", "GO CC") %in% colnames(seqAnnoDF))
}


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



getGOparents = function(id, onto="BP"){
  go2Ancestor = switch (onto, BP=as.list(GOBPANCESTOR), MF=as.list(GOMFANCESTOR), CC=as.list(GOCCANCESTOR))
  ancestorIds = setdiff(go2Ancestor[[id]], "all") ## remove the all category
  if (is.null(ancestorIds)){
    return()
  } else {
    ancestorTerms = Term(ancestorIds[GOTERM])
    return(ancestorTerms)
  }
}

addGoParents = function(gene2goList, onto){
  goParents = switch(onto, BP=as.list(GOBPPARENTS),
                     CC=as.list(GOCCPARENTS),
                     MF=as.list(GOMFPARENTS))
  return(lapply(gene2goList, function(x){setdiff(union(x, unlist(goParents[x])), "all")}))
}


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


ezGroupGO = function (selectedGenes, go2GeneList, onto="CC", levels = 2:4, goSlim=NULL) {
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
ezGetGoByLevels = function (ont, levels, goSlim=NULL) {
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





goStringsToList = function(goStrings, listNames=NULL){
  x = strsplit(goStrings, "; ")
  names(x) = listNames
  x = lapply(x, function(x){ifelse(x == "", character(0), x)})
  return(x)
}


## ezGoseq considers only the genes that have annotations; genes without annotation are removed from the selectedGenes and allGenes
ezGoseq = function(param, selectedGenes, allGenes, gene2goList=NULL, 
                    method=c("Wallenius", "Sampling", "Hypergeometric"),
                    onto=NULL, normalizedAvgSignal=NULL, verbose=FALSE){
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



myFisherTest = function(goGenes, selGenes, allGenes, alternative="greater"){ 
  fisher.test(factor(allGenes %in% selGenes, levels=c("FALSE", "TRUE")), 
              factor(allGenes %in% goGenes, levels=c("FALSE", "TRUE")), alternative=alternative)$p.value
}


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
                                          param$minCountFisher, onto=onto, maxNumberOfTerms=param$maxNumberGroupsDisplayed)
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


goResultToHtmlTable2 = function(goResults, pThreshGo, minCount, onto=NA, maxNumberOfTerms=40) {
  
  upTerms = .getGoTermsAsTd(goResults[["enrichUp"]], pThreshGo, minCount, onto, color="red", maxNumberOfTerms=maxNumberOfTerms)
  bothTerms = .getGoTermsAsTd(goResults[["enrichBoth"]], pThreshGo, minCount, onto, color="black", maxNumberOfTerms=maxNumberOfTerms)
  downTerms = .getGoTermsAsTd(goResults[["enrichDown"]], pThreshGo, minCount, onto, color="blue", maxNumberOfTerms=maxNumberOfTerms)
  terms = c(upTerms, bothTerms, downTerms)
  lines = paste(paste("<tr>", terms, "</tr>"), collapse="\n")
  table = paste("<table border='1'><tr><th>Term</th><th>p</th><th>N</th></tr>", lines, "</table>")
  table
}

goResultToHtmlTable = function(x, pThreshGo, minCount, onto=NA, maxNumberOfTerms=40){
  
  terms = .getGoTermsAsTd(x, pThreshGo, minCount, onto, maxNumberOfTerms=maxNumberOfTerms)
  if (length(terms) >= 1){
    rows = paste(paste("<tr>", terms, "</tr>"), collapse="\n")
    table = paste("<table border='1'><tr><th>Term</th><th>p</th><th>N</th></tr>", rows, "</table>")
  } else {
    table = ""
  }
  
  return(table)
}



.getGoTermsAsTd = function(x, pThreshGo, minCount, onto, color="black", valign="top", maxNumberOfTerms=40){
  
  require(GO.db)
  
  rows = character(0) #"<td> </td><td> </td><td> </td>"
  if (!is.data.frame(x)){
    return(rows)
  }
  x = x[x$Count >= minCount & x$Pvalue < pThreshGo, ]
  x = x[order(x$Pvalue), ]
  if (nrow(x) > maxNumberOfTerms){
    x = x[1:maxNumberOfTerms, ]
  }
  if (nrow(x) == 0){
    #return("<tr><td><i>no significants</i></td><td></td><td></td></tr>")
    return(character(0))
  }
  
  if (onto == "CC"){
    ANCESTOR = GOCCANCESTOR
    OFFSPRING = GOCCOFFSPRING
    CHILDREN = GOCCCHILDREN
  }
  if (onto == "BP"){
    ANCESTOR = GOBPANCESTOR
    OFFSPRING = GOBPOFFSPRING
    CHILDREN = GOBPCHILDREN
  }
  if (onto == "MF"){
    ANCESTOR = GOMFANCESTOR
    OFFSPRING = GOMFOFFSPRING
    CHILDREN = GOMFCHILDREN
  }
  
  goIds = rownames(x)
  goAncestorList = as.list(ANCESTOR[goIds])
  
  goRoots = character()
  for (goId in goIds){
    if (length(intersect(goIds, goAncestorList[[goId]])) == 0){
      goRoots[goId] = goId
    }
  }
  goOffsprings = unique(as.list(OFFSPRING[goIds]))[[1]]
  goAncestors = unique(unlist(goAncestorList))
  goRelatives = union(intersect(goAncestors, goOffsprings), goIds)
  
  terms = character(length(goRoots))
  pValues = character(length(goRoots))
  counts = character(length(goRoots))
  for (i in 1:length(goRoots)){
    childTerms = getChildTerms(goRoots[i], goIds, goRelatives, indent="", CHILDREN)
    terms[i] = paste0("<span title='", getGOTerm(childTerms)[[1]], "'>", names(childTerms), "</span>",
                     collapse="<br>")
    pValues[i] = paste(signif(x[childTerms, "Pvalue"], 3), collapse="<br>")
    counts[i] = paste(x[childTerms, "Count"], x[childTerms, "Size"], sep="/", collapse="<br>")
  }
  rows = paste0("<td style='color:", color, "' valign='", valign, "'>", terms, "</td>",
               "<td style='color:", color, "' valign='", valign, "'>", pValues,"</td>",
               "<td style='color:", color, "' valign='", valign, "'>", counts, "</td>")
  rows
}



getChildTerms = function(x, subset, goRelatives, indent="", childEnvir){
  
  result = character()
  if (is.element(x, subset)){
    term = getGOTerm(x)[[1]]
    displayLabel = paste0(indent, substr(term, 1, 30))
    result[displayLabel] = x
    indent = paste0("<font color='#FFFFFF'>", indent, "</font>")
    subset = setdiff(subset, x) ## this is the modification to make the table non-redundant!!
  }
  kids = intersect(childEnvir[[x]], goRelatives)
  indent = paste(indent, ".")
  for (kid in kids){
    kidResult = getChildTerms(kid, subset, goRelatives, indent=indent, childEnvir)
    subset = setdiff(subset, kidResult) ## again the non-redundancy modification
    result = c(result, kidResult)
  }
  result
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

getTopGoHtmlTable = function(goTermStat, pValueThresh=NA, catSize=NA){
  
  use = goTermStat$Pvalue < pValueThresh & goTermStat$Annotated > catSize & goTermStat$Significant > goTermStat$Expected
  result = paste("<td>", goTermStat$Term[use], "</td><td>", signif(goTermStat$Pvalue[use], digits=3), "</td>")
  result = paste(paste("<tr>", result, "</tr>"), collapse="\n")
  result = paste("<table>", result, "</table")
  result
}


clusterHeatmap = function(param, x, file="cluster-heatmap.png", nClusters=5, lim=c(-4, 4),
                          colColors=NULL, d=NULL, columnDist=NULL, doClusterColumns=FALSE,
                          clusterColors=rainbow(nClusters), doGO=TRUE, ontologies=c("BP", "MF", "CC"), 
                          seqAnno=NULL,
                          universeGeneIds=NULL, universeProbeIds=NULL, method="ward.D2", 
                          cexRow=1.0, cexCol=1.5, labRow=rownames(x), width=max(800, 400 + 10 * ncol(x)), height=1000, margins=c(14,9),
                          colors=ezRedBlueScale(256), maxGenesWithLabel=50, keggOrganism=NA, ...){
  
  library(gplots, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  
  
  if (is.null(d)){
    xx = x
    xx[is.na(x)] = min(x, na.rm=TRUE)
    d = dist(xx)
  }
  hcl = hclust(d, method=method)
  probeDendro = as.dendrogram(hcl)
  probeDendro = reorder(probeDendro, rowMeans(x, na.rm=TRUE))
  clusterNumbers = cutree(hcl, k=nClusters)
  
  if (doClusterColumns){
    if (is.null(columnDist)){
      columnDist = dist(t(x))
    }
    hcl = hclust(columnDist, method=method)
    colDendro = as.dendrogram(hcl)
    colDendro = reorder(colDendro, colMeans(x, na.rm=TRUE))
    showDendro = "both"
  } else {
    colDendro =FALSE
    showDendro = "row"
  }
  
  
  result = list()
  result$nClusters = nClusters
  result$clusterNumbers = clusterNumbers
  result$clusterColors = clusterColors
  
  if (param$showGeneClusterLabels & nrow(x) < maxGenesWithLabel){
    if (is.null(labRow)){
      labRow = rownames(x)    
    }
  } else {
    labRow = ""
    margins[2] = 0.5 
  }
  
  
  if (!is.null(file)){
    png(file=file, width=width, height=height)
    on.exit(dev.off())
  }
  if (!is.null(colColors)){
    if (nClusters > 1){
      heatmap.2(x,
                ColSideColors=colColors, RowSideColors=clusterColors[clusterNumbers],
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
    if (nClusters > 1){		
      heatmap.2(x,
                RowSideColors=clusterColors[clusterNumbers],
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
  
  if (doGO && doGo(param, seqAnno)){
    library(GOstats, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
    library(annotate, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
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
      #ezWrite("universeProbeIds:\n", str(universeProbeIds))
      #ezWrite("p2g[universeProbeIds]:\n", str(unlist( probe2Gene[universeProbeIds])))
      universeGeneIds =na.omit(unique(unlist( genes[universeProbeIds]) ))
      #ezWrite("universeGeneids:\n", str(universeGeneIds))
    }
    keggClusterResults = NULL
    if (!is.na(keggOrganism)){
      kegg2GenesList = getPathway2GenesList(param, keggOrganism)
      if (all(unlist(kegg2GenesList)) %in% universeGeneIds){
        keggClusterResults = list()
      }
    }
    genesByCluster = tapply(genes[rownames(x)], clusterNumbers, function(x){na.omit(unique(x))}, simplify=FALSE)
    goClusterResults = ezMclapply(ontologies, function(onto){
      gene2go = goStringsToList(seqAnno[[paste("GO", onto)]], listNames=rownames(seqAnno))[universeGeneIds]
      if (param$includeGoParentAnnotation){
        gene2go = addGoParents(gene2go, onto)
      }
      result = lapply(genesByCluster, function(selectedGenes){
        ezGoseq(param, selectedGenes, universeGeneIds, gene2goList=gene2go,
                 method="Hypergeometric",
                 onto=onto, normalizedAvgSignal=NULL, verbose=FALSE)})
    }, mc.cores=1)
    names(goClusterResults) = ontologies
    result$GO = goClusterResults
  }    
  return(invisible(result))
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
