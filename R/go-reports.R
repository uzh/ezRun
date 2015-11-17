###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets the GO cluster table
##' @description Gets the GO cluster table 
##' @param param a list of parameters to extract \code{pValThreshFisher} and \code{minCountFisher} from
##' @param clusterResult a list containing the result of the analysis done by \code{goClusterResults()}.
##' @template roxygen-template
##' @seealso \code{\link{goClusterResults}}
##' @return Returns a flex table containing the GO information of the cluster result.
goClusterTable = function(param, clusterResult){
  ontologies = names(clusterResult$GO)
  tables = ezMatrix("", rows=paste("Cluster", 1:clusterResult$nClusters), cols=ontologies)
  for (onto in ontologies){
    for (i in 1:clusterResult$nClusters){
      x = clusterResult$GO[[onto]][[i]]
      goFrame = .getGoTermsAsTd(x, param$pValThreshFisher, param$minCountFisher, onto=onto)
      if (nrow(goFrame) > 0){
        tables[i, onto] = as.html(ezFlexTable(goFrame))
      }
    }
  }
  nameMap = c("BP"="Biological Proc. (BP)", "MF"="Molecular Func. (MF)", "CC"="Cellular Comp. (CC)")
  colnames(tables) = nameMap[colnames(tables)]
  ft = ezFlexTable(tables, border = 2, header.columns = TRUE, add.rownames=TRUE)
  bgColors = rep(gsub("FF$", "", clusterResult$clusterColors), each=ncol(tables)+1)
  ft = setFlexTableBackgroundColors(ft, colors=bgColors)
  return(ft)
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
  doc = addParagraph(doc, paste("Maximum number of terms displayed:", param$maxNumberGroupsDisplayed))
  
  doc = addTitle(doc, "GO categories that are overrepresented among significantly upregulated genes.", 3)
  doc = addFlexTable(doc, udt$flexTables[["enrichUp"]])
  doc = addTitle(doc, "GO categories that are overrepresented among significantly downregulated genes.", 3)
  doc = addFlexTable(doc, udt$flexTables[["enrichDown"]])
  doc = addTitle(doc, "GO categories that are overrepresented among all significantly regulated genes.", 3)
  doc = addFlexTable(doc, udt$flexTables[["enrichBoth"]])
  
  #     goFiles = list.files('.',pattern='enrich.*txt')
  #     revigoLinks = ezMatrix("", rows=c('Both', 'Down', 'Up'), cols=c('BP', 'CC', 'MF'))
  #     for(j in 1:length(goFiles)){
  #       goResult = ezRead.table(goFiles[j], sep='\t')[, c('GO ID','Pvalue')]
  #       goResult = goResult[which(goResult$Pvalue < param$pValThreshFisher),]
  #       if(nrow(goResult) > param$maxNumberGroupsDisplayed) {
  #         goResult = goResult[1:param$maxNumberGroupsDisplayed,]
  #       }
  #       revigoLinks[j] = paste0('http://revigo.irb.hr/?inputGoList=',
  #                               paste(goResult[,'GO ID'], goResult[,'Pvalue'], collapse='%0D%0A'))
  #       revigoLinks[j] = pot("ReViGO Link", hyperlink = revigoLinks[j])
  #     }
  #     titles[["ReViGO"]] = "ReViGO"
  #     addTitleWithAnchor(doc, titles[[length(titles)]], 3)
  #     doc = addFlexTable(doc, ezFlexTable(cbind(rownames(revigoLinks), revigoLinks), valign="middle", header.columns=TRUE))

    ## addTxtLinksToReport(doc, param, udt$txtFiles)
  if (param$doZip){
    addTxtLinksToReport(doc, udt$txtFiles, mime="application/zip")
  } else {
    addTxtLinksToReport(doc, udt$txtFiles, mime="application/txt")
  }
  
  
  doc
}

##' @describeIn addGoUpDownResult Gets the GO up-down tables.
goUpDownTables = function(param, goResult){
  goTable = ezMatrix("", rows="Cats", cols=names(goResult))
  resultList = list("enrichUp"=goTable, "enrichBoth"=goTable, "enrichDown"=goTable)
  txtFiles = character() ## TODO make a list of list; similar to resultList
  ## txtList = list("enrichUp"=list(), "enrichBoth"=list(), "enrichDown"=list())
  for (onto in names(goResult)){ ## BP, MF , CC
    x = goResult[[onto]]
    for (sub in names(x)){ #c("enrichUp", "enrichDown", "enrichBoth")){
      message("sub: ", sub)
      goFrame = .getGoTermsAsTd(x[[sub]], param$pValThreshFisher, param$minCountFisher, onto=onto,
                                maxNumberOfTerms=param$maxNumberGroupsDisplayed)
      if (nrow(goFrame) > 0){
        resultList[[sub]]["Cats", onto] = as.html(ezFlexTable(goFrame))
      } else {
        message("no rows")
      }
      xSub = x[[sub]]
      if (is.data.frame(xSub)){
        name = paste0(onto, "-", param$comparison, "-", sub)
        if (!is.null(xSub$Pvalue)){
          xSub = xSub[order(xSub$Pvalue), ]
          xSub = cbind("GO ID"=rownames(xSub), xSub)
        }
        txtFile = ezValidFilename(paste0(name, ".txt"), replace="-")
        # txtList[[sub]][[onto]] = ezValidFilename(paste0(name, ".txt"), replace="-")
        ezWrite.table(xSub, file=txtFile, row.names=FALSE)
        if (param$doZip){
          txtFiles[name] = zipFile(txtFile)
        } else {
          txtFiles[name] = txtFile
        }
      }
    }
  }
  flexTables = lapply(resultList, function(res){
    ezFlexTable(res, border = 2, header.columns = TRUE)
  })
  return(list(flexTables=flexTables, txtFiles=txtFiles))
}

##' @describeIn goClusterTable Gets the GO terms and pastes them into a table.
.getGoTermsAsTd = function(x, pThreshGo, minCount, onto=NA, maxNumberOfTerms=40){
  
  require(GO.db)
  
  if (!is.data.frame(x)){
    message("got not data frame")
    return(ezFrame("Term"=character(0), "p"=numeric(0), "N"=integer(0)))
  }
  x = x[x$Count >= minCount & x$Pvalue < pThreshGo, ]
  x = x[order(x$Pvalue), ]
  if (nrow(x) > maxNumberOfTerms){
    x = x[1:maxNumberOfTerms, ]
  }
  if (nrow(x) == 0){
    return(ezFrame("Term"=character(0), "p"=numeric(0), "N"=integer(0)))
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
  
  terms = character()
  pValues = character()
  counts = character()
  for (i in 1:length(goRoots)){
    childTerms = getChildTerms(goRoots[i], goIds, goRelatives, indent="", CHILDREN)
    for (term in childTerms){
      terms = append(terms, names(childTerms)[childTerms==term])
      pValues = append(pValues, signif(x[term, "Pvalue"], 3))
      counts = append(counts, paste(x[term, "Count"], x[term, "Size"], sep="/"))
    }
  }
  return(ezFrame("Term"=terms, "p"=pValues, "N"=counts))
}
