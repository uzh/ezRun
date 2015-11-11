###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


addGageTables = function(doc, param = NULL, gageResults = NULL) {
  doc = addParagraph(doc, paste("Gene sets used:", paste(names(gageResults[['all']]), collapse=", ")))
  if(any(grepl('kg', names(gageResults[['all']])))) {
    doc = addParagraph(doc, as.html(pot("<span style='margin-left:2em'>kg = <A HREF='http://www.genome.jp/kegg/pathway.html'>KEGG</A>
                                        pathways: dise (disease pathways) , sigmet (signaling or metabolism pathways)")))
  }
  if(any(grepl('msigdb', names(gageResults[['all']])))) {
    doc = addParagraph(doc, as.html(pot("<span style='margin-left:2em'>msigdb =
                                        <A HREF='http://www.broadinstitute.org/gsea/msigdb/index.jsp'>MSigDB</A> pathway")))
  }
  doc = addParagraph(doc, paste0("Significance threshold pathways: ", param[['gageThreshold']],
                                 ". Only pathways below this treshold are represented."))
  doc = addParagraph(doc, paste0("Significance threshold genes selected within a pathway: ", param[['gageGeneThreshold']],
                                 ". Only genes below this treshold are represented"))
  doc = addParagraph(doc, "Warning : only pathways with at least one gene significant will be displayed. Only top 30 pathways are represented")
  
  gene.pValue=param[['gageGeneThreshold']]
  
  for (i in names(gageResults[['significant']])) {
    #     use = paste0("<table border=0><tr><th>Heatmap Plot logRatio Signal for ",i,"</th><th>",i," significant pathways</th></tr>")
    checkpoint = nrow(gageResults[['significant']][[i]][['combined']]) > 0 | nrow(gageResults[['significant']][[i]][['both']]) > 0
    if(!checkpoint) next
    #     ezWrite(use, con=html)
    tableRows = list()
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
      res <- res[,!colnames(res) %in% "links", drop=F]
      
      
      # Writing plot and table
      #       ezWrite("<tr valign=top><td>", con=html)
      #       writeImageRowToHtml(x[[lab.png]], con=html)
      #       ezWrite("</td><td>", con=html)
      imgrow = x[[lab.png]]
      pathColors = unique(x[[lab.pathCol]])
      #       writeTableToHtmlWhite(res, con=html, 
      #                             bgcolors=matrix(gsub("FF$", "", unique(pathColors)), nrow=length(unique(pathColors)), ncol=1))
      #       ezWrite("</td></tr>", con=html)
      tbl = ezAddTableWhite(res, bgcolors=matrix(gsub("FF$", "", unique(pathColors)), nrow=length(unique(pathColors)), ncol=1))
      tableRows[[signal]] = cbind(imgrow, tbl)
    }
    #     ezWrite("</table>", con=html)
    table = ezGrid(rbind(unlist(tableRows)), header.columns=TRUE)
    table = addHeaderRow(table, cbind(paste("Heatmap Plot logRatio Signal for", i), paste(i, "significant pathways")))
    doc = addFlexTable(doc, table)
  }
}
