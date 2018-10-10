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
