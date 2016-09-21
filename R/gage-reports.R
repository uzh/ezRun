###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Adds the gage tables
##' @description Adds the gage tables to an html file.
##' @template doc-template
##' @templateVar object tables
##' @param param a list of parameters, possibly passed to other functions as well:
##' \itemize{
##'  \item{gageThreshold}{ the threshold for significant pathways.}
##'  \item{gageGeneThreshold}{ the threshold for genes selected within a pathway.}
##'  \item{pathview}{ a logical indicating whether to add links to the kegg pathview.}
##' }
##' @param gageResults the results of the gage analysis.
##' @seealso \code{\link{gageAnalysis}}
##' @seealso \code{\link{twoGroupsGO}}
##' @template roxygen-template
addGageTables = function(doc, param = NULL, gageResults = NULL) {
  addParagraph(doc, paste("Gene sets used:", paste(names(gageResults[['all']]), collapse=", ")))
  if(any(grepl('kg', names(gageResults[['all']])))) {
    addParagraph(doc, pot("kg = KEGG pathways: dise (disease pathways), sigmet (signaling or metabolism pathways)",
                                hyperlink="http://www.genome.jp/kegg/pathway.html"))
  }
  if(any(grepl('msigdb', names(gageResults[['all']])))) {
    addParagraph(doc, pot("msigdb = MSigDB pathway", hyperlink="http://www.broadinstitute.org/gsea/msigdb/index.jsp"))
  }
  addParagraph(doc, paste0("Significance threshold pathways: ", param[['gageThreshold']],
                                 ". Only pathways below this treshold are represented."))
  addParagraph(doc, paste0("Significance threshold genes selected within a pathway: ", param[['gageGeneThreshold']],
                                 ". Only genes below this treshold are represented"))
  addParagraph(doc, "Warning : only pathways with at least one gene significant will be displayed. Only top 30 pathways are represented")
  
  gene.pValue=param[['gageGeneThreshold']]
  
  outerTable = ezFrame(Heatmaps=character(0), "Significant Pathways"=character(0))
  for (i in names(gageResults[['significant']])) {
    checkpoint = nrow(gageResults[['significant']][[i]][['combined']]) > 0 | nrow(gageResults[['significant']][[i]][['both']]) > 0
    if(!checkpoint) next
    #tableRows = list()
    for (signal in c("combined", "both")) {
      lab.expr = paste(signal, 'expr', sep='.')
      lab.sigGenes = paste(signal, 'sigGenes', sep='.')
      lab.pValue = paste(signal, 'pValue', sep='.')
      lab.png = paste(signal, 'png', sep='.')
      lab.pathCol = paste(signal, "pathColors", sep=".")
      
      x = gageResults[['significant']][[i]]
      res = x[[signal]]
      
      # Exit if no sigGenes
      if(nrow(x[[lab.sigGenes]])==0) next
      
      # Change numbers
      formatCol = c("p.geomean","stat.mean","p.val","q.val","exp1")
      res[,formatCol] = formatC(as.numeric(res[,formatCol]), digits=3)
      
      # Add links
      links = x[["links"]][rownames(res),]
      res = cbind(res,links)
      sigGenes = x[[lab.sigGenes]][["Set"]]
      
      # Add number of significant genes per pathway
      # Warning : only pathways with at least one gene significant will be represented!
      SigGenes = table(sigGenes)
      res = res[row.names(res) %in% names(SigGenes),,drop=F]
      res = cbind(res, SigGenes = SigGenes[row.names(res)])
      res = cbind(res, SigProp = formatC(as.numeric(res[,"SigGenes"])/as.numeric(res[,"set.size"]), digits=3))
      
      # Add links to kegg pathview
      gset.name = x$gset.name
      if (param[['pathview']] & grepl("^kg", x$gset.name)) {
        pathways = rownames(res)
        SpeciesName = getSpeciesName(param)
        kegg.id = getKeggId(SpeciesName, param)
        pattern = paste0(kegg.id,"[[:digit:]]+")
        kegg.pathId = unlist(regmatches(pathways, gregexpr(pattern, pathways)))
        pngFiles = paste0(kegg.pathId,".",x$gset.name,"-",signal,".png")
        pngLinks = paste0("<A  HREF='",pngFiles,"'>click here</A>")
        res = cbind(res, PathView = pngLinks)
      }
      
      # Add links to pathway names
      rownames(res) = paste0("<A HREF=\"", res[,"links"],"\">",rownames(res),"</A>")
      res <- res[, colnames(res) != "links", drop=F]
      
      # Writing plot and table
      resTable = ezFlexTable(res, talign="right", header.columns = TRUE, add.rownames=TRUE)
      bgColors = gsub("FF$", "", unique(x[[lab.pathCol]]))
      resTable = setFlexTableBackgroundColors(resTable, j=1, colors=bgColors)
      
      outerTable[signal, ] = c(x[[lab.png]], as.html(resTable))
      # imgLink = x[[lab.png]]
      # pathColors = unique(x[[lab.pathCol]])
      # tbl = ezAddTableWhite(res, bgcolors=matrix(gsub("FF$", "", unique(pathColors)), nrow=length(unique(pathColors)), ncol=1))
      #tableRows[[signal]] = cbind(imgLink, tbl)
    }
    if (nrow(outerTable) > 0){
      addFlexTable(doc, ezGrid(outerTable, add.rownames=TRUE))
    }
    # table = ezGrid(rbind(unlist(tableRows)), header.columns=TRUE)
    # table = addHeaderRow(table, cbind(paste("Heatmap Plot logRatio Signal for", i), paste(i, "significant pathways")))
    # addFlexTable(doc, table)
  }
}
