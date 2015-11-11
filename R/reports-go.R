###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


goClusterTable = function(param, clusterResult){
  ontologies = names(clusterResult$GO)
  tables = ezMatrix("", rows=paste("Cluster", 1:clusterResult$nClusters), cols=ontologies)
  for (onto in ontologies){
    for (i in 1:clusterResult$nClusters){
      x = clusterResult$GO[[onto]][[i]]
      tables[i, onto] = goResultToHtmlTable(x, param$pValThreshFisher, param$minCountFisher, onto=onto) ## TODOP: REFAC
    }
  }
  ft = ezFlexTable(tables, border = 2, header.columns = TRUE)
  bgColors = rep(gsub("FF$", "", clusterResult$clusterColors), each=ncol(tables))
  ft = setFlexTableBackgroundColors(ft, colors=bgColors)
  return(ft)
}


addGoUpDownResult = function(doc, param, goResult){
  udt = goUpDownTables(param, goResult)
  
  doc = addParagraph(doc, "Red GO categories are overrepresented among the significantly upregulated genes")
  doc = addParagraph(doc, "Blue GO categories are overrepresented among the significantly downregulated genes")
  doc = addParagraph(doc, "Black GO categories are overrepresented among all signifcantly regulated genes")
  doc = addParagraph(doc, paste("Maximum number of terms displayed:", param$maxNumberGroupsDisplayed))
  doc = addFlexTable(doc, udt$flexTable)
  
  if (param$doZip){
    addTxtLinksToReport(doc, udt$txtFiles, mime="application/zip")
  } else {
    addTxtLinksToReport(doc, udt$txtFiles, mime="application/txt")
  }
}


goUpDownTables = function(param, goResult){
  tables = ezMatrix("", rows="Cats", cols=names(goResult))
  txtFiles = character()
  for (onto in names(goResult)){
    x = goResult[[onto]]
    tables[1,onto] = goResultToHtmlTable2(x, param$pValThreshFisher, ## TODOP: REFAC
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
  ft = ezFlexTable(tables, border = 2, header.columns = TRUE)
  return(list(flexTable=ft, txtFiles=txtFiles))
}
