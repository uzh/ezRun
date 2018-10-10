goClusterTable = function(param, clusterResult, seqAnno){
  ontologies = names(clusterResult$GO)
  tables = ezMatrix("", rows=paste("Cluster", 1:clusterResult$nClusters), cols=ontologies)
  linkTable = ezMatrix("", rows = 1:clusterResult$nClusters, cols = ontologies)
  enrichrTable = ezMatrix("", rows = 1:clusterResult$nClusters, cols = "Enrichr")
  for (i in 1:clusterResult$nClusters){
    genesToUse = rownames(seqAnno) %in% names(clusterResult$clusterNumbers)[clusterResult$clusterNumbers==i]
    genesList = paste(seqAnno$gene_name[genesToUse], collapse="\\n")
    jsCall = paste0('enrich({list: "', genesList, '", popup: true});')
    enrichrTable[i, 1] = as.html(pot(paste0("<a href='javascript:void(0)' onClick='", jsCall, "'>Enrichr</a>")))
    for (onto in ontologies){
      x = clusterResult$GO[[onto]][[i]]
      goFrame = .getGoTermsAsTd(x, param$pValThreshFisher, param$minCountFisher, onto=onto)
      if (nrow(goFrame)==0) next
      linkTable[i, onto] = paste0("Cluster-", onto, "-", i, ".html")
      ezInteractiveTable(goFrame, tableLink=linkTable[i, onto], digits=3,
                         title=paste("GO categories of cluster", i, "and ontology", onto))
      linkTable[i, onto] = as.html(ezLink(linkTable[i, onto], target="_blank"))
      goFrame$Term = substr(goFrame$Term, 1, 30)
      if (nrow(goFrame) > 0){
        tables[i, onto] = as.html(ezFlexTable(goFrame, talign="right", header.columns = TRUE))
      }
    }
  }
  nameMap = c("BP"="Biological Proc. (BP)", "MF"="Molecular Func. (MF)", "CC"="Cellular Comp. (CC)")
  colnames(tables) = nameMap[colnames(tables)]
  ft = ezFlexTable(tables, border = 2, header.columns = TRUE, add.rownames=TRUE)
  bgColors = rep(gsub("FF$", "", clusterResult$clusterColors))
  ft = setFlexTableBackgroundColors(ft, j=1, colors=bgColors)
  return(list(ft=ft, linkTable=linkTable, enrichrTable=enrichrTable))
}


##' @title Adds the GO up-down results
##' @description Adds the GO up-down results to an html file.
##' @template doc-template
##' @templateVar object result
##' @param param a list of parameters to pass to \code{goUpDownTables()} and extract \code{doZip} from.
##' @param goResult the GO result to get the up-down results from. Can be obtained by \code{twoGroupsGO()}.
##' @seealso \code{\link{twoGroupsGO}}
##' @template roxygen-template
addGoUpDownResult = function(doc, param, goResult){
  udt = goUpDownTables(param, goResult)
  addParagraph(doc, paste("Maximum number of terms displayed:", param$maxNumberGroupsDisplayed))
  
  addFlexTable(doc, ezFlexTable(udt$linkTable, add.rownames=TRUE))
  addTitle(doc, "GO categories that are overrepresented among significantly upregulated genes.", 3)
  addFlexTable(doc, udt$flexTables[["enrichUp"]])
  addTitle(doc, "GO categories that are overrepresented among significantly downregulated genes.", 3)
  addFlexTable(doc, udt$flexTables[["enrichDown"]])
  addTitle(doc, "GO categories that are overrepresented among all significantly regulated genes.", 3)
  addFlexTable(doc, udt$flexTables[["enrichBoth"]])
  
  revigoLinks = ezMatrix("", rows=c('enrichBoth', 'enrichDown', 'enrichUp'), cols=c('BP', 'CC', 'MF'))
  for (col in colnames(revigoLinks)){
    for (row in rownames(revigoLinks)){
      goSubResult = goResult[[col]][[row]]
      if (all(is.na(goSubResult))) next
      goSubResult = goSubResult[which(goSubResult$Pvalue < param$pValThreshFisher),]
      if(nrow(goSubResult) > param$maxNumberGroupsDisplayed) {
        goSubResult = goSubResult[1:param$maxNumberGroupsDisplayed,]
      }
      revigoLinks[row, col] = paste0('http://revigo.irb.hr/?inputGoList=',
                                     paste(rownames(goSubResult), goSubResult[,'Pvalue'], collapse='%0D%0A'))
      revigoLinks[row, col] = as.html(pot("ReViGO", hyperlink = revigoLinks[row, col]))
    }
  }
  revigoTitle = "ReViGO"
  addTitle(doc, revigoTitle, 3, id=revigoTitle)
  addFlexTable(doc, ezFlexTable(revigoLinks, valign="middle", header.columns = TRUE, add.rownames = TRUE))
  addTxtLinksToReport(doc, udt$txtFiles, param$doZip)
  return(revigoTitle)
}

