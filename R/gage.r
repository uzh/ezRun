###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


runGageAnalysis = function(result, param=NULL, output=NULL, rawData=NULL, gene.pValue=param[['gageGeneThreshold']]) {
  
  #require(gage, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  #require(pathview, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  #require(reshape2, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  
  require(gage)
  require(pathview)
  require(reshape2)
  
  stopifnot(param$featureLevel == "gene")
  
  # Load Gene sets (KEGG and Human MSigDb)
  geneSets = getGeneSets(param)
  # Run GAGE
  gageResults = gageAnalysis(result, param=param, rawData=rawData, geneSets=geneSets)
  # Get normalized expression and P-values from genes in all pathways
  gageResults[['all']] = getExpressionGage(gageResults[['all']], result, rawData, 
                                           param, signal = c("greater", "less", "both"))
  
  # Get Significant Genes for each signal for significant pathways
  gageResults[['all']] = lapply(gageResults[['all']], gageSigGenes, gene.pValue=gene.pValue, signal='both')
  gageResults[['all']] = lapply(gageResults[['all']], gageSigGenes, gene.pValue=gene.pValue, signal='greater')
  gageResults[['all']] = lapply(gageResults[['all']], gageSigGenes, gene.pValue=gene.pValue, signal='less')
  
  # Write results txt table
  writeGageResults(gageResults, param=param, output=output, prefix='gage')
  # Get normalized expression and P-values from genes in significant pathways for posterior data visualization
  gageResults[['significant']] = getExpressionGage(gageResults[['significant']], result, rawData, 
                                                   param, signal = c("greater", "less", "combined", "both"))

  # Get Significant Genes for each signal for significant pathways
  gageResults[['significant']] = lapply(gageResults[['significant']], gageSigGenes, gene.pValue=gene.pValue, signal='combined')
  gageResults[['significant']] = lapply(gageResults[['significant']], gageSigGenes, gene.pValue=gene.pValue, signal='both')
  gageResults[['significant']] = lapply(gageResults[['significant']], gageSigGenes, gene.pValue=gene.pValue, signal='greater')
  gageResults[['significant']] = lapply(gageResults[['significant']], gageSigGenes, gene.pValue=gene.pValue, signal='less')
  
  # Create HeatMaps Plots
  for (i in names(geneSets)) {
    for (signal in c("combined", "both")) {
      prefix =  paste0('gage-heatmap',"-", i)
      lab.pathCol = paste(signal, "pathColors", sep=".")
      lab.png = paste(signal, "png", sep=".")
      res = gageHeatmap(gageResults[['significant']][[i]],param = param, output=output, gene.pValue=gene.pValue, signal=signal, prefix=prefix)
      gageResults[['significant']][[i]][[lab.pathCol]] = res$pathwayColors
      gageResults[['significant']][[i]][[lab.png]] = basename(as.character(res[["fileName"]]))
    }
  }
  
  # PathView
  kgSets = grep("^kg", names(geneSets), value = T)
  if(param[['pathview']] & length(kgSets)>0 ) {
    SpeciesName = getSpeciesName(param)
    kegg.id = getKeggId(SpeciesName, param)
    for (i in kgSets) {
      for (signal in c("combined", "both")) {
        #message(paste(i,signal))
        gagePathview(gageResults[['significant']][[i]], param = param, 
                     output=output, signal=signal, result = result, 
                     anno = rawData$seqAnno, kegg.id = kegg.id)
      }
    }
  }
  return(gageResults)
}

readGmt = function (file) {
  f <- readLines(file)
  lst = sapply(f, function(x) unlist(strsplit(x, "\t", fixed = TRUE)))
  names(lst) = sapply(lst, function(x) x[1])
  lst = lapply(lst, function(x) x[-(1:2)])
  return(lst)
}

getSpeciesName = function(param){
  buildFields = strsplit(param[["refBuild"]], "/", fixed=TRUE)[[1]]
  return(buildFields[1])
}

getKeggId = function(x, param) {
  kegg.org.info = read.delim(file.path(param[['KEGGdb']], param[['KEGGOrgId']]))
  
  # Query is in kegg.org.info
  if(x %in% kegg.org.info$fgcz) {
    kegg.id = kegg.org.info$KeggOrg[kegg.org.info$fgcz==x]
    return(as.character(kegg.id))
    stop()
  }
  
  # If not in the original, check if the name is in the full name field using "_" as a separator
  strOrgName = gsub("_", " ", x)
  
  if(any(grepl(strOrgName, kegg.org.info$FullName, ignore.case = T))) {
    kegg.id = kegg.org.info$KeggOrg[grep(strOrgName, kegg.org.info$FullName, ignore.case = T)]
    return(as.character(kegg.id))
  } else {
    warning(paste("Genome build not find in the KEGG organism list. Please add genome manually to",param[['KEGGOrgId']]))
    return(NA)
  }
}

getGeneSets = function(param) {
  SpeciesName = getSpeciesName(param)
  kegg.id = getKeggId(SpeciesName, param)
  
  # Load KEGG sets
  if(is.na(kegg.id)) stop('Specie unknown in kegg table')
  
  lf = grep("symbol", list.files(param[['KEGGgmt']], pattern=kegg.id, full.names=T), value=T)
  
  if(grepl(",", param[['genesets']])) param[['genesets']] <- unlist(strsplit(param[['genesets']],","))

  if(length(lf)>0) {
    # Read precomputed KEGG gene sets with HGNC symbol
    geneSet = lapply(lf, readGmt)
    keggSetNames = gsub(paste(kegg.id,"_",".symbol.gmt",sep="|"),"",basename(lf))
    names(geneSet) = keggSetNames
    geneSet = geneSet[names(geneSet) %in% param[['genesets']]]
  } else {
    stop(paste('No KEGG database found at', param[['KEGGgmt']]))
    ## Load data from kegg.gsets
    #geneSet = kegg.gsets(species = kegg.id, id.type = "entrez")
    
    ## Transform to HGNC symbol
    #require(biomaRt)
    #ensembl=useMart("ensembl")
    #listDatasets(ensembl)
    #dataset <- "hsapiens_gene_ensembl"
    #ensembl <- useMart("ensembl", dataset = dataset)
    #res <- lapply(keggSet$kg.sets, 
    #              function(x) getBM(attributes='hgnc_symbol', 
    #                                filters ='entrezgene', 
    #                                values = x, 
    #                                mart = ensembl)[,'hgnc_symbol'])
    #keggSet$kg.symbol <- res
    
    
  }
  
  # Adding MSigDB if Homo Sapiens
  if(SpeciesName == 'Homo_sapiens') {
    # READ MSigDB (read gmt file)
    use = unlist(strsplit(param[['MSigDB']], ','))
    for (i in use) {
      gmtFile = file.path(param[['MSigDBPath']],i)
      msigdb = readGmt(gmtFile)
      msigName = paste0('msigdb.', i)
      geneSet[[msigName]] = msigdb
    }
  }
  return(geneSet)
}

gageAnalysis = function(result, rawData=NULL, param=NULL, geneSets=NULL ) {
  
  # Prepare gage input data
  gage.input = as.matrix(result[[param[['gageInput']]]])
  seqAnno = rawData$seqAnno
  probes = names(result$pValue)
  dimnames(gage.input)[[1]] = sapply(seqAnno[match(probes, rownames(seqAnno)), "gene_name"], as.character)
  gage.input = gage.input[result$usedInTest,]
  
  # Run gage and get all the results + add html links
  geneSets.names = names(geneSets)

#   lapply(seq_along(geneSets), function(i) paste(names(geneSets)[[i]], geneSets[[i]]))
  
  fc.gsets.all = lapply(seq_along(geneSets), function(i, data=gage.input, paramI=param) {
    x = geneSets[[i]]
    res = gage(data, gsets = x, ref = NULL, samp = NULL)
    both = gage(data, gsets = x, ref = NULL, samp = NULL, same.dir = FALSE)
    res$both = both$greater
    res$statsboth = both$stats
    res$gsets = x
    res$gset.name = names(geneSets)[[i]]
    res$links = matrix(rep(NA,length(names(x))), dimnames=list(names(x), "link"))
    htmlbase = NA
    if(grepl("^kg\\.", res$gset.name)) {
      htmlbase = paramI[['KEGGhtml']]
      keggGsetId = gsub(" .*","",names(x))
      res$links =  matrix(paste0(htmlbase,keggGsetId), dimnames=list(names(x), "link"))
    } 
    if(grepl("^msigdb", res$gset.name)) {
      htmlbase = paramI[['MSigDBhtml']]
      res$links =  matrix(file.path(htmlbase,names(x)), dimnames=list(names(x), "link"))
    }
    return(res)
  })
  names(fc.gsets.all) <- geneSets.names
  
  
  # Filter gage results using param[['gageThreshold']] q value
  fc.gsets = lapply(fc.gsets.all, function(x, cutoff=param[['gageThreshold']], qpval = "q.val") {
    greater.pval = x$greater[,qpval]
    use = greater.pval <= cutoff & !is.na(greater.pval)
    x$greater = x$greater[use, , drop=F]
    x$greater = x$greater[order(x$greater[,qpval]), , drop=F]
    x$greater = suppressWarnings(cbind(x$greater, signal='greater'))
    if(nrow(x$greater) > 30) x$greater = x$greater[1:30,]
    
    less.pval = x$less[,qpval]
    use = less.pval <= cutoff & !is.na(less.pval)
    x$less = x$less[use, , drop=F]
    x$less = x$less[order(x$less[,qpval]), , drop=F]
    x$less = suppressWarnings(cbind(x$less, signal='less'))
    if(nrow(x$less) > 30) x$less = x$less[1:30,]
    
    x$combined = rbind(x$greater,x$less)
    combined.pval = x$combined[,qpval]
    use = order(as.numeric(combined.pval))
    x$combined = x$combined[use, , drop=F]
    if(nrow(x$combined) > 30) x$combined = x$combined[1:30,]
    
    both.pval = x$both[,qpval]
    use = both.pval <= cutoff & !is.na(both.pval)
    x$both = x$both[use, , drop=F]
    x$both = x$both[order(x$both[,qpval]), , drop=F]
    x$both = suppressWarnings(cbind(x$both, signal='both'))
    if(nrow(x$both) > 30) x$both = x$both[1:30,]
    
    return(x)
  })
  
  
  return(list(all=fc.gsets.all, significant=fc.gsets))
  
}

writeGageResults = function(gageResults, param=NULL, output=NULL, prefix=NULL, signal=c("greater", "less", "both") ) {
  gageAll = gageResults[['all']]
  outDir = "."
  if(is.null(prefix)) prefix = "gage"
  
  for(gs in names(gageAll)) {
    x = gageAll[[gs]]
    for (test in signal) {
      outfile = paste(paste(param[['comparison']],prefix, gs, test, sep = "-"), "txt", sep=".")
      outfile = file.path(outDir, outfile)
      res = x[[test]]
      links = x[["links"]][rownames(res),]
      res = cbind(res,links)
      
      fSet = row.names(res)
      fGenes = lapply(fSet, 
             function(y, datf=x$less.sigGenes){
               paste(datf[datf$Set==y,"Gene"],collapse=";") 
             })
      
      res = cbind(res, unlist(fGenes))
      colnames(res)[ncol(res)] = paste0("genes(",param[['gageGeneThreshold']],")")
      ezWrite.table(res, file=outfile)
    }
  }
}

getExpressionGage = function(gageResults, result=NULL, rawData=NULL, param = NULL, signal=NULL) {
  xgetExpressionGage = function (x, gset, anno = NULL, expr = NULL, pValue = NULL, param = NULL) {
    xNames = row.names(x)
    genes_symbols = unique(unlist(gset[xNames]))
    genes_symbols = intersect(genes_symbols, anno$gene_name)
    ## TODO: currently works only for ENSEMBL genome builds
    gsEnsembl = rownames(anno)[match(genes_symbols, anno$gene_name)]
    
    exprRes = expr[gsEnsembl, ,drop=FALSE]
    rownames(exprRes) = genes_symbols
    exprRes = shiftZeros(exprRes, param[['minSignal']])
    
    pValueRes = pValue[gsEnsembl]
    names(pValueRes) = genes_symbols
    
    res = list(expr = exprRes, pValue = pValueRes)
    return(res)
  }
  
  fc.gsets = lapply(gageResults, function(x, 
                                          xanno = rawData$seqAnno, 
                                          xexpr = result$xNorm, 
                                          xpValue = result$fdr, 
                                          xparam = param, 
                                          xsignal = signal) {
    for (i in xsignal) {
      expName = paste(i, "expr", sep=".")
      pValueName = paste(i, "pValue", sep=".")
      x[[expName]] = matrix(nrow = 0, ncol = 0)
      x[[pValueName]] = matrix(nrow = 0, ncol = 0)
      if(nrow(x[[i]]) > 0) {
        ge = xgetExpressionGage(x[[i]], x$gsets, anno = xanno, expr = xexpr, pValue = xpValue, param = xparam)
        x[[expName]] = ge$expr
        x[[pValueName]] = ge$pValue
      }
    }
    return(x)
    
    
  })
}

gageSigGenes = function(x, gene.pValue=NULL, signal=NULL ){
  # Select significant genes
  pValueVar = paste(signal,'pValue',sep='.')
  use = x[[pValueVar]] <= gene.pValue & !is.na(x[[pValueVar]])
  genesUse = names(x[[pValueVar]])[use]
  
  label = paste(signal,"sigGenes",sep='.')
  x[[label]] = matrix(nrow=0, ncol=2, dimnames=list(c(),c("Gene", "Set")))
  if(!is.null(genesUse)) {
    xAnnot = melt(x$gsets[rownames(x[[signal]])])
    xAnnot = xAnnot[xAnnot$value %in% genesUse,]
    names(xAnnot) = c("Gene", "Set")
    x[[label]] = xAnnot
  }
  return(x)
  
}

gageHeatmap = function(x, param=NULL, output=NULL, gene.pValue=NULL, signal=NULL, fileName=NULL, prefix='gage-heatmap'){
  library(gplots, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  
  sampleColors = getSampleColors(param$grouping)[order(param$grouping)]
  
  # Select Significant genes
  x = gageSigGenes(x, gene.pValue, signal=signal)
  
  # Define labels for list elements
  lab.expr = paste(signal, 'expr', sep='.')
  lab.sigGenes = paste(signal, 'sigGenes', sep='.')
  lab.pValue = paste(signal, 'pValue', sep='.')
  
  # Exit if no sigGenes
  if(nrow(x[[lab.sigGenes]])==0) return(list(pathwayColors=NA, fileName=NA))
  
  # Get expression
  use = as.character(x[[lab.sigGenes]]$Gene)
  
  xlogSignal = log2(shiftZeros(x[[lab.expr]][use,,drop=F], param$minSignal))
  xCentered = (xlogSignal - rowMeans(xlogSignal))[, order(param$grouping), drop=FALSE]
  xPathway = x[[lab.sigGenes]]$Set
  
  
  # nGroups = nlevels(factor(param$grouping))
  
  # Pathways colors
  xPathway = factor(xPathway, levels=unique(xPathway))
  pathwayColsPal = rainbow(nlevels(xPathway), start=0.1, end=0.8, s=0.8, v=0.7)
  pathwayColors = pathwayColsPal[as.numeric(xPathway)]
  
  # HeatMap
  if(is.null(fileName) & nrow(xCentered) >= 2 & ncol(xCentered) >= 2) {
    fileName = paste0(param[['comparison']],"-", prefix,"-", signal, ".png")
    
    png(fileName, width=max(800, 400 + 10 * ncol(xCentered)), height=1000)  ## REFAC
    rowDendro = FALSE
    colDendro = T
    showDendro = "column"
    heatmap.2(xCentered,
              col=ezRedBlueScale(256),
              ColSideColors=sampleColors, RowSideColors=pathwayColors,
              scale="none",
              Colv=colDendro, Rowv=rowDendro, dendrogram=showDendro,
              key=TRUE, density.info="none", trace="none",
              keysize=1, cexCol=1.5,
              margins=c(14,9), cexRow = 0.00001) ## TODO: why????
    
    dev.off()
  }  
  return(list(pathwayColors=pathwayColors, fileName=fileName))
}

gagePathview = function(x, param=NULL, output=NULL, signal=NULL, kegg.id=NULL, gene.pValue=NULL, result = result, anno = rawData$seqAnno){
  library(pathview, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  
  # Define labels for list elements
  lab.sigGenes = paste(signal, 'sigGenes', sep='.')
  
  # Extract significant genes ID
  genes = x[[lab.sigGenes]]
  
  if(nrow(genes)==0) return(NULL)
  
  # Add LogRatio to significant genes
  gsEnsembl = rownames(anno)[match(genes$Gene, anno$gene_name)]
  use = match(gsEnsembl, names(result$pValue))
  log2Ratio = result$log2Ratio[use]
  names(log2Ratio) = genes$Gene
  genes$log2Ratio = log2Ratio
  
  # Convert to list
  genes = split(genes, factor(genes$Set))
  genes = lapply(genes, function(x){res=x$log2Ratio; names(res)=x$Gene; return(res)})
  
  # Convert pathway ID
  pathways = names(genes)
  pattern = paste0(kegg.id,"[[:digit:]]+")
  kegg.pathId = unlist(regmatches(pathways, gregexpr(pattern, pathways)))
  names(genes) = kegg.pathId
  
  suffix = paste(x$gset.name, signal, sep="-")
  
  for (kegg.path.id in names(genes)) {
    #message(cat(signal,kegg.path.id,"\n"))
    cat(signal,kegg.id,kegg.path.id,"\n")
    pv.out = try(suppressMessages(pathview(gene.data = genes[[kegg.path.id]], 
                                           pathway.id = kegg.path.id, 
                                           species = kegg.id, 
                                           out.suffix = suffix, 
                                           kegg.native = T,
                                           kegg.dir = param[['KEGGxml']],
                                           gene.idtype = "SYMBOL",
                                           limit = list(gene = 3, cpd = 3))), silent = TRUE)
    
  }  
}

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



