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
##' @param seqAnno the sequence annotation.
##' @template roxygen-template
##' @seealso \code{\link{goClusterResults}}
##' @return Returns a flex table containing the GO information of the cluster result.

goClusterTableRmd = function(param, clusterResult, seqAnno){
  ontologies = names(clusterResult$GO)
  ktables = list()
  linkTable = ezMatrix("", rows = 1:clusterResult$nClusters, cols = ontologies)
  enrichrTable = ezMatrix("", rows = 1:clusterResult$nClusters, cols = "Enrichr")
  for (i in 1:clusterResult$nClusters){
    genesToUse = rownames(seqAnno) %in% names(clusterResult$clusterNumbers)[clusterResult$clusterNumbers==i]
    genesList = paste(seqAnno$gene_name[genesToUse], collapse="\\n")
    jsCall = paste0('enrich({list: "', genesList, '", popup: true});')
    enrichrTable[i, 1] = paste0("<a href='javascript:void(0)' onClick='", jsCall, 
                                "'>Enrichr</a>")
    ## Prepare the table for kable
    ktableCluster <- list()
    for (onto in ontologies){
      x = clusterResult$GO[[onto]][[i]]
      goFrame = .getGoTermsAsTd(x, param$pValThreshFisher, param$minCountFisher, 
                                onto=onto)
      ktableCluster[[onto]] <- goFrame
      if (nrow(goFrame)==0)
        next
      linkTable[i, onto] = paste0("Cluster-", onto, "-", i, ".html")
      ezInteractiveTable(goFrame, tableLink=linkTable[i, onto], digits=3,
                         title=paste("GO categories of cluster", i, "and ontology", onto))
      linkTable[i, onto] = ezLink(linkTable[i, onto], target="_blank")
      goFrame$Term = substr(goFrame$Term, 1, 30)
    }
    ## This is some ugly code to append some "" cell, so they can used in kable
    maxNrow <- max(sapply(ktableCluster, nrow))
    ktableCluster <- lapply(ktableCluster, 
                            function(x){rbind(as.matrix(x), 
                                              ezMatrix("", rows=seq_len(maxNrow-nrow(x)), 
                                                       cols=seq_len(ncol(x))))}
                            )
    ktableCluster <- do.call(cbind, ktableCluster)
    if(nrow(ktableCluster) == 0L){
      ## for later grouping in cluster kables, we need empty cells.
      ktableCluster <- ezMatrix("", rows=1, cols=colnames(ktableCluster))
    }
    ktables[[i]] <- ktableCluster
  }
  return(list(ktables=ktables, linkTable=linkTable, enrichrTable=enrichrTable))
}


revigoUpDownTables <- function(param, goResult){
  revigoLinks = ezMatrix("", rows=c('enrichUp', 'enrichDown', 'enrichBoth'), 
                         cols=c('BP', 'MF', 'CC'))
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
      revigoLinks[row, col] <- kableExtra::text_spec("ReViGO", format="html", 
                                                     link = revigoLinks[row, col])
    }
  }
  return(t(revigoLinks))
}

##' @describeIn addGoUpDownResult Gets the GO up-down tables.
goUpDownTables = function(param, goResult){
  goTable <- list()
  ktables = list("enrichUp"=goTable, "enrichDown"=goTable, "enrichBoth"=goTable)
  txtFiles = character() ## TODO make a list of list; similar to resultList
  linkTable = ezMatrix("", rows = names(goResult), 
                       cols = c("enrichUp", "enrichDown", "enrichBoth"))
  for (onto in names(goResult)){ ## BP, MF , CC
    x = goResult[[onto]]
    for (sub in names(x)){ #c("enrichUp", "enrichDown", "enrichBoth")){
      message("sub: ", sub)
      xSub = x[[sub]]
      if (is.data.frame(xSub)){
        ## We always output the goseq results files
        name = paste0(onto, "-", param$comparison, "-", sub)
        if (!is.null(xSub$Pvalue)){
          xSub = xSub[order(xSub$Pvalue), ]
          xSub = cbind("GO ID"=rownames(xSub), xSub)
        }
        txtFile = ezValidFilename(paste0(name, ".txt"), replace="-")
        txtFiles <- append(txtFiles, txtFile)
        ezWrite.table(xSub, file=txtFile, row.names=FALSE)
      }
      goFrame = .getGoTermsAsTd(xSub, param$pValThreshFisher,
                                param$minCountFisher, onto=onto,
                                maxNumberOfTerms=param$maxNumberGroupsDisplayed)
      ktables[[sub]][[onto]] = goFrame
      if (nrow(goFrame)==0)
        next
      linkTable[onto, sub] = paste0("Cluster-", onto, "-", sub, ".html")
      ezInteractiveTable(goFrame, tableLink=linkTable[onto, sub], digits=3,
                         title=paste(sub("enrich", "", sub), 
                                     "enriched GO categories of ontology", onto))
      linkTable[onto, sub] = ezLink(linkTable[onto, sub], 
                                            target = "_blank")
    }
  }
  for(sub in names(ktables)){
    ### Add the ""
    maxNrow <- max(sapply(ktables[[sub]], nrow))
    ktables[[sub]] <- lapply(ktables[[sub]],
                             function(x){rbind(as.matrix(x),
                                               ezMatrix("", rows=seq_len(maxNrow-nrow(x)),
                                                        cols=seq_len(ncol(x))))}
                             )
    ktables[[sub]] <- do.call(cbind, ktables[[sub]])
  }
  return(list(ktables=ktables, txtFiles=txtFiles, linkTable=linkTable))
}

##' @describeIn goClusterTable Gets the GO terms and pastes them into a table.
.getGoTermsAsTd = function(x, pThreshGo, minCount, onto=NA, maxNumberOfTerms=40){
  
  require("GO.db")
  require(AnnotationDbi)
  
  if (!is.data.frame(x)){
    message("got no data frame")
    return(ezFrame("Term"=character(0), "ID"=character(0), 
                   "p"=numeric(0), "N"=integer(0)))
  }
  x = x[x$Count >= minCount & x$Pvalue < pThreshGo, ]
  x = x[order(x$Pvalue), ]
  if (nrow(x) > maxNumberOfTerms){
    x = x[1:maxNumberOfTerms, ]
  }
  if (nrow(x) == 0){
    return(ezFrame("Term"=character(0), "ID"=character(0),
                   "p"=numeric(0), "N"=integer(0)))
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
  goAncestorList = AnnotationDbi::as.list(ANCESTOR[goIds]) ## without the explicit choice of AnnotationDbi:: this fails in RnaBamStats ..... no idea why
  
  goRoots = character()
  for (goId in goIds){
    if (length(intersect(goIds, goAncestorList[[goId]])) == 0){
      goRoots[goId] = goId
    }
  }
  goOffsprings = unique(AnnotationDbi::as.list(OFFSPRING[goIds]))[[1]]
  goAncestors = unique(unlist(goAncestorList))
  goRelatives = union(intersect(goAncestors, goOffsprings), goIds)
  
  terms = character()
  ids = character()
  pValues = numeric()
  counts = character()
  for (i in 1:length(goRoots)){
    childTerms = getChildTerms(goRoots[i], goIds, goRelatives, indent="", CHILDREN)
    for (term in childTerms){
      terms = append(terms, names(childTerms)[childTerms==term])
      ids = append(ids, term)
      pValues = append(pValues, x[term, "Pvalue"])
      counts = append(counts, paste(x[term, "Count"], x[term, "Size"], sep="/"))
    }
  }
  return(ezFrame("Term"=terms, "ID"=ids,"p"=pValues, "N"=counts))
}
