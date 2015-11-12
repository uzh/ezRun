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
      tables[i, onto] = goResultToHtmlTable(x, param$pValThreshFisher, param$minCountFisher, onto=onto)
    }
  }
  ft = ezFlexTable(tables, border = 2, header.columns = TRUE)
  bgColors = rep(gsub("FF$", "", clusterResult$clusterColors), each=ncol(tables))
  ft = setFlexTableBackgroundColors(ft, colors=bgColors)
  return(ft)
}

## TODOP: REFAC
##' @describeIn goClusterTable Parses the GO result into a html table.
goResultToHtmlTable = function(x, pThreshGo, minCount, onto=NA, maxNumberOfTerms=40){
  
  terms = .getGoTermsAsTd(x, pThreshGo, minCount, onto, maxNumberOfTerms=maxNumberOfTerms)
  if (length(terms) >= 1){
    rows = paste(paste("<tr>", terms, "</tr>"), collapse="\n")
    table = paste("<table border='1'><tr><th>Term</th><th>p</th><th>N</th></tr>", rows, "</table>")
  } else {
    table = " "
  }
  return(table)
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

##' @describeIn addGoUpDownResult Gets the GO up-down tables.
goUpDownTables = function(param, goResult){
  tables = ezMatrix("", rows="Cats", cols=names(goResult))
  txtFiles = character()
  for (onto in names(goResult)){
    x = goResult[[onto]]
    tables[1, onto] = goResultToHtmlTable2(x, param$pValThreshFisher,
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

## TODOP: REFAC, obsolete after refactor 1 table 3 colors into 3 seperate tables.
##' @describeIn addGoUpDownResult Parses the GO result into a html table.
goResultToHtmlTable2 = function(goResults, pThreshGo, minCount, onto=NA, maxNumberOfTerms=40) {
  
  upTerms = .getGoTermsAsTd(goResults[["enrichUp"]], pThreshGo, minCount, onto, color="red", maxNumberOfTerms=maxNumberOfTerms)
  bothTerms = .getGoTermsAsTd(goResults[["enrichBoth"]], pThreshGo, minCount, onto, color="black", maxNumberOfTerms=maxNumberOfTerms)
  downTerms = .getGoTermsAsTd(goResults[["enrichDown"]], pThreshGo, minCount, onto, color="blue", maxNumberOfTerms=maxNumberOfTerms)
  terms = c(upTerms, bothTerms, downTerms)
#   table = {"Term"=, "p"=, "N"=} # only feasably refactorable after .getGoTermsAsTd has been refactored
#   ft = ezFlexTable(table, header.columns = TRUE)
  lines = paste(paste("<tr>", terms, "</tr>"), collapse="\n")
  table = paste("<table border='1'><tr><th>Term</th><th>p</th><th>N</th></tr>", lines, "</table>")
  return(table)
}

## TODOP: REFAC and doc
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
  # return(ezFrame("Term"=terms, "p"=pValues, "N"=counts))
  rows
}
