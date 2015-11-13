###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets the GO cluster table
##' @description Gets the GO cluster table 
##' @param param a list of parameters to extract \code{pValThreshFisher} and \code{minCountFisher} from
##' @param clusterResult a list containing the result of the analysis done by \code{clusterHeatmap()}.
##' @template roxygen-template
##' @seealso \code{\link{clusterHeatmap}}
##' @return Returns a flex table containing the GO information of the cluster result.
goClusterTable = function(param, clusterResult){
  ontologies = names(clusterResult$GO)
  tables = ezMatrix("", rows=paste("Cluster", 1:clusterResult$nClusters), cols=ontologies)
  for (onto in ontologies){
    for (i in 1:clusterResult$nClusters){
      x = clusterResult$GO[[onto]][[i]]
      tables[i, onto] = .getGoTermsAsTd(x, param$pValThreshFisher, param$minCountFisher, onto=onto)
    }
  }
  ft = ezFlexTable(tables, border = 2, header.columns = TRUE)
  bgColors = rep(gsub("FF$", "", clusterResult$clusterColors), each=ncol(tables))
  ft = setFlexTableBackgroundColors(ft, colors=bgColors)
  return(ft)
}

##' @title Adds the GO up-down results
##' @description Adds the GO up-down results to an html file.
##' @param doc an object of the class bsdoc to add the result.
##' @param param a list of parameters to pass to \code{goUpDownTables()} and extract \code{doZip} from.
##' @param goResult the GO result to get the up-down results from. Can be obtained by \code{twoGroupsGO()}.
##' @seealso \code{\link{twoGroupsGO}}
##' @template roxygen-template
addGoUpDownResult = function(doc, param, goResult){
  udt = goUpDownTables(param, goResult)
  doc = addParagraph(doc, paste("Maximum number of terms displayed:", param$maxNumberGroupsDisplayed))
  
  doc = addParagraph(doc, "GO categories that are overrepresented among significantly upregulated genes.")
  doc = addFlexTable(doc, udt$flexTable[["up"]])
  doc = addParagraph(doc, "GO categories that are overrepresented among significantly downregulated genes.")
  doc = addFlexTable(doc, udt$flexTable[["down"]])
  doc = addParagraph(doc, "GO categories that are overrepresented among all significantly regulated genes.")
  doc = addFlexTable(doc, udt$flexTable[["both"]])
  
  if (param$doZip){
    addTxtLinksToReport(doc, udt$txtFiles, mime="application/zip")
  } else {
    addTxtLinksToReport(doc, udt$txtFiles, mime="application/txt")
  }
}

##' @describeIn addGoUpDownResult Gets the GO up-down tables.
goUpDownTables = function(param, goResult){
  tables = ezMatrix("", rows="Cats", cols=names(goResult))
  tables = list("enrichUp"=tables, "enrichBoth"=tables, "enrichDown"=tables)
  txtFiles = character()
  for (onto in names(goResult)){
    x = goResult[[onto]]
    for (sub in names(x)){ #c("enrichUp", "enrichDown", "enrichBoth")){
      tables[[sub]][1, onto] = .getGoTermsAsTd(x[[sub]], param$pValThreshFisher,
                                             param$minCountFisher, onto=onto, maxNumberOfTerms=param$maxNumberGroupsDisplayed)
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
  ftUp = ezFlexTable(tables[["enrichUp"]], border = 2, header.columns = TRUE)
  ftBoth = ezFlexTable(tables[["enrichBoth"]], border = 2, header.columns = TRUE)
  ftDown = ezFlexTable(tables[["enrichDown"]], border = 2, header.columns = TRUE)
  return(list(flexTables=list("up"=ftUp, "both"=ftBoth, "down"=ftDown), txtFiles=txtFiles))
}

##' @describeIn goClusterTable Gets the GO terms and pastes them into a table.
.getGoTermsAsTd = function(x, pThreshGo, minCount, onto=NA, maxNumberOfTerms=40){
  
  require(GO.db)
  
  rows = character(0)
  if (!is.data.frame(x)){
    return(rows)
  }
  x = x[x$Count >= minCount & x$Pvalue < pThreshGo, ]
  x = x[order(x$Pvalue), ]
  if (nrow(x) > maxNumberOfTerms){
    x = x[1:maxNumberOfTerms, ]
  }
  if (nrow(x) == 0){
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
    terms[i] = getGOTerm(childTerms)[[1]]  ## some entries might be too long
    pValues[i] = signif(x[childTerms, "Pvalue"], 3)
    counts[i] = paste(x[childTerms, "Count"], x[childTerms, "Size"], sep="/")
  }
  return(ezFrame("Term"=terms, "p"=pValues, "N"=counts))
}
