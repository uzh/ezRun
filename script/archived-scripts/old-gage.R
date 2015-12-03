


writeGageTables = function(html, param = NULL, gageResults = NULL) {
  ezWrite("<h3>GAGE Enrichment Analysis</h3>", con=html)
  use = paste0("<p>Gene sets used: ", paste(names(gageResults[['all']]), collapse=", "), "<br>")
  ezWrite(use, con=html)
  if(any(grepl('kg', names(gageResults[['all']])))) {
    use = "<p><span style='margin-left:2em'>kg = <A HREF='http://www.genome.jp/kegg/pathway.html'>KEGG</A> pathways: dise (disease pathways) , sigmet (signaling or metabolism pathways)<br>"
    ezWrite(use, con=html)
  }
  if(any(grepl('msigdb', names(gageResults[['all']])))) {
    use = "<p><span style='margin-left:2em'>msigdb = <A HREF='http://www.broadinstitute.org/gsea/msigdb/index.jsp'>MSigDB</A> pathway<br>"
    ezWrite(use, con=html)
  }
  use = paste0("<p>Significance threshold pathways: ", param[['gageThreshold']], ". Only pathways below this treshold are represented<br>")
  ezWrite(use, con=html)
  use = paste0("<p>Significance threshold genes selected within a pathway: ", param[['gageGeneThreshold']], ". Only genes below this treshold are represented<br>")
  ezWrite(use, con=html)
  ezWrite("<p>Warning : only pathways with at least one gene significant will be displayed. Only top 30 pathways are represented<br>", con=html)  
  
  
  gene.pValue=param[['gageGeneThreshold']]
  
  for (i in names(gageResults[['significant']])) {
    use = paste0("<table border=0><tr><th>Heatmap Plot logRatio Signal for ",i,"</th><th>",i," significant pathways</th></tr>")
    checkpoint = nrow(gageResults[['significant']][[i]][['combined']]) > 0 | nrow(gageResults[['significant']][[i]][['both']]) > 0
    if(!checkpoint) next
    ezWrite(use, con=html)
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
      if(param[['pathview']] & grepl("^kg", x$gset.name)) {
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
      ezWrite("<tr valign=top><td>", con=html)
      writeImageRowToHtml(x[[lab.png]], con=html)
      ezWrite("</td><td>", con=html)
      pathColors = unique(x[[lab.pathCol]])
      writeTableToHtmlWhite(res, con=html, 
                            bgcolors=matrix(gsub("FF$", "", unique(pathColors)), 
                                            nrow=length(unique(pathColors)), ncol=1))
      ezWrite("</td></tr>", con=html)
    }
    ezWrite("</table>", con=html)
  }
}


writeTableToHtmlWhite = function(x, con=stdout(), bgcolors=matrix("#ffffff", nrow=nrow(x), ncol=ncol(x)),
                                 valign="middle", border=1, head=""){
  
  ezWrite("<table border='", border, "'><tr>", con=con)
  ezWrite(paste0("<th>", c(head, colnames(x)), "</th>", collapse="\n"), con=con)
  ezWrite("</tr>", con=con)
  if (nrow(x) > 0){
    for (i in 1:nrow(x)){
      ezWrite("<tr><th>", rownames(x)[i], "</th>", con=con)
      ezWrite(paste0("<td valign='", valign, "' bgcolor='", bgcolors[i,], "'><font color='white'>", x[i,], "</font></td>", collapse="\n"), con=con)
      ezWrite("</tr>", con=con)
    }
  }
  ezWrite("</table>", con=con)
}
