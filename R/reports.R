###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Wrapper for \code{FlexTable()}
##' @description Wraps \code{FlexTable()} with defaults to remove the cell header and cell borders.
##' @param x a matrix or data.frame to turn into an object of the class FlexTable.
##' @param header a logical indicating whether to use a header for the table.
##' @template addargs-template
##' @templateVar fun FlexTable
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{FlexTable}}
##' @return Returns an object of the class FlexTable.
##' @examples
##' ezFlexTable(cbind(a=1:5,b=11:15))
ezFlexTable = function(x, header=FALSE, ...){
  FlexTable(x, body.cell.props = cellProperties(border.width = 0),
            header.cell.props = cellProperties(border.width = 0),
            header.columns = header, ...)
}

# how to add help text? for each plot seperately or not?
##' @title Gets an image link as html
##' @description Gets an image link as html. Also plots and creates the image.
##' @param ezPlotter an object of the class EzPlotter or inheriting from it.
##' @param file a character specifying the name of the image with a .png suffix.
##' @param mouseOverText a character specifying the text being displayed when mousing over the image.
##' @param addPdfLink a logical indicating whether to add a link on the image to a pdf version of itself.
##' @param width an integer specifying the width of each plot to create an image from.
##' @param height an integer specifying the height of each plot to create an image from.
##' @template addargs-template
##' @templateVar fun ezPlotter
##' @template roxygen-template
##' @seealso \code{\link{EzPlotter}}
##' @return Returns a character specifying a link to an image in html.
##' @examples
##' imageLink = ezImageFileLink(EzPlotterIris$new())
##' theDoc = bsdoc(title = 'My document')
##' theDoc = addParagraph(theDoc, imageLink)
##' writeDoc(theDoc, "example.html")
ezImageFileLink = function(ezPlotter, file=NULL, mouseOverText=ezPlotter$mouseOverText,
                           addPdfLink=TRUE, width=480, height=480, ...){
  pngName = ezPlotter$plotPng(file=file, width=width, height=height, ...)
  if (addPdfLink) {
    pdfName = ezPlotter$plotPdf(file=sub(".png$", ".pdf", file), width=width, height=height, ...)
    imgFilePot = pot(paste('<img src="', pngName, '" title="', mouseOverText, '"/>'),
                     hyperlink = pdfName)
  } else {
    imgFilePot = pot(paste('<img src="', pngName, '" title="', mouseOverText, '"/>'))
  }
  return(as.html(imgFilePot))
}

# important: creating the report like this will not use a connection, so this causes many adjustments
# currently not possible to use from old report opener:
# writeLines(readLines(ezCSSFile()), con=html) ## no custom css for reporteRs
# writeJavaScriptIgvStarter(htmlFile, param$projectId, html) ## perhaps refactorable with addJavascript(), but how to put in html head?
##' @title Opens an html report
##' @description Opens an html report using \code{bsdoc()} from the ReporteRs package. Also adds some introductory elements.
##' @param title a character specifying the title of the html report.
##' @param dataset usually a data.frame from the meta field of an EzDataset.
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{bsdoc}}
##' @seealso \code{\link[ReporteRs]{writeDoc}}
##' @return Returns an object of the class bsdoc to add further elements.
##' @examples
##' theDoc = openBsdocReport(title="My html report")
##' ezAddBootstrapMenu(theDoc)
##' closeBsdocReport(doc=theDoc, file="example.html")
openBsdocReport = function(title="", dataset=NULL){
  doc = bsdoc(title = title)
  pot1 = pot(paste("Started on", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "--&#160;"))
  pot2 = as.html(pot("Documentation", hyperlink = "http://fgcz-sushi.uzh.ch/doc/methods-20140422.html"))
  doc = addFlexTable(doc, ezFlexTable(cbind(pot1, pot2)))
  doc = addTitleWithAnchor(doc, title)
  if (!is.null(dataset)){
    ezWrite.table(dataset, file="dataset.tsv", head="Name")
    doc = addParagraph(doc, pot("dataset.tsv", hyperlink = "dataset.tsv"))
  }
  return(doc)
}

##' @describeIn openBsdocReport Adds a paragraph showing the finishing time and writes the document. \code{file} must have a .html suffix.
closeBsdocReport = function(doc, file){
  doc = addParagraph(doc, paste("Finished", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  writeDoc(doc, file=file)
}

##' @describeIn openBsdocReport Adds a bootstrapmenu with the FGCZ link and optional navigation. Section links can be added to the Navigation with \code{#Anchorname}.
ezAddBootstrapMenu = function(doc, titles=NULL){
  bootStrap = BootstrapMenu("Functional Genomics Center Zurich", link = "http://www.fgcz.ethz.ch")
  if (ezIsSpecified(titles)){
    ddMenu = DropDownMenu("Navigation")
    for (each in titles){
      ddMenu = addLinkItem(ddMenu, label=sub("#", "", each), link=each)
    }
    bootStrap = addLinkItem(bootStrap, dd=ddMenu)
  }
  doc = addBootstrapMenu(doc, bootStrap)
}

##' @title Adds a title with an anchor
##' @description Adds a title with an anchor by using \code{addTitle()} and \code{addCodeBlock()} from the ReporteRs package.
##' @param doc an object of the class bsdoc to add the anchored title to.
##' @param title a character specifying the title.
##' @param level an integer specifying the heading level.
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{addTitle}}
##' @seealso \code{\link[ReporteRs]{addCodeBlock}}
##' @examples
##' theDoc = openBsdocReport(title="My html report")
##' title = "My title"
##' ezAddBootstrapMenu(theDoc, paste0("#", title))
##' addTitleWithAnchor(theDoc, title)
##' closeBsdocReport(doc=theDoc, file="example.html")
addTitleWithAnchor = function(doc, title, level=1){
  doc = addParagraph(doc, as.html(pot(paste0("<a name='", title, "'></a>"))))
  doc = addTitle(doc, value=title, level=level)
}

##' @title Writes an error report
##' @description Writes an error report to an html file. Also creates the file and closes it.
##' @param htmlFile a character representing the path to write the file in. Must have a .html suffix.
##' @param param a list of parameters to extract the \code{name} from.
##' @param dataset usually a data.frame from the meta field of an EzDataset.
##' @param error a character vector representing the error message(s).
##' @template roxygen-template
##' @seealso \code{\link{openBsdocReport}}
##' @seealso \code{\link{closeBsdocReport}}
##' @examples
##' param = ezParam()
##' htmlFile = "example.html"
##' writeErrorReport(htmlFile, param)
writeErrorReport = function(htmlFile, param=param, dataset=NULL, error="Unknown Error"){
  html = openBsdocReport(title=paste("Error:", param$name), dataset=dataset)
  html = ezAddBootstrapMenu(html)
  html = addTitle(html, "Error message", level=2)
  for (i in 1:length(error)){
    html = addParagraph(html, error[i])
  }
  closeBsdocReport(html, htmlFile)
}

##' @title Writes a text link
##' @description Writes a text link to an html file.
##' @param txtNames a character representing the link name.
##' @param mime a character representing the type of the link.
##' @template connection-template
##' @template roxygen-template
##' @examples
##' param = ezParam()
##' htmlFile = "example_html"
##' html =  openHtmlReport(htmlFile,param)
##' writeTxtLinksToHtml("dataset.tsv", con=html)
addTxtLinksToReport = function(txtNames, mime="text/plain", doc){
  for (each in txtNames){
    doc = addParagraph(doc, pot(paste("<a href='", each, "' type='", mime, "'>", each, "</a>")))
  }
}

##' @title Adds a table
##' @description Adds a table to a bsdoc object.
##' @param x a matrix or data.frame to paste a table from.
##' @param doc an object of the class bsdoc to add the table to.
##' @param bgcolors a matrix specifying the background colors.
##' @param valign a character specifying where to align the table elements vertically. Use either "top", "middle" or "bottom".
##' @param border an integer specifying the border width.
##' @param head a character specifying the contents of the upperleft corner of the table.
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{addFlexTable}}
##' @examples
##' x = matrix(1:25,5)
##' rownames(x) = letters[1:5]
##' colnames(x) = LETTERS[1:5]
##' html = openBsdocReport()
##' addTableToReport(x, html, head="Example", bgcolors="red")
##' closeBsdocReport(html, "example.html")
addTableToReport = function(x, doc, bgcolors=NULL, valign="middle", border=1, head=""){
  if (is.null(bgcolors)){
    bgcolors = matrix("#ffffff", nrow=nrow(x), ncol=ncol(x))
  }
  x = cbind(rownames(x),x)
  bodyCells = cellProperties(border.width=border, vertical.align=valign)
  table = FlexTable(x, header.columns = FALSE, body.cell.props=bodyCells,
                    header.cell.props=cellProperties(border.width = border))
  table = setFlexTableBackgroundColors(table, j=2:length(colnames(x)), colors=bgcolors)
  table = addHeaderRow(table, c(head, colnames(x)[2:length(colnames(x))]))
  doc = addFlexTable(doc, table)
}

##' @describeIn addTableToReport Does the same with a white font and returning the table instead of adding it to the document.
addTableToReportWhite = function(x, doc, bgcolors=NULL, valign="middle", border=1, head=""){
  if (is.null(bgcolors)){
    bgcolors = matrix("#ffffff", nrow=nrow(x), ncol=ncol(x))
  }
  x = cbind(rownames(x),x)
  x = as.html(pot(paste('<font color="white">', x, '</font>')))
  bodyCells = cellProperties(border.width=border, vertical.align=valign)
  table = FlexTable(x, header.columns = FALSE, body.cell.props=bodyCells,
                    header.cell.props=cellProperties(border.width = border))
  table = setFlexTableBackgroundColors(table, j=2:length(colnames(x)), colors=bgcolors)
  table = addHeaderRow(table, c(head, colnames(x)[2:length(colnames(x))]))
  return(table)
}

##' @title Adds a summary of the count result
##' @description Adds a summary of the count result to a bsdoc object.
##' @param doc an object of the class bsdoc to add the table to.
##' @param param a list of parameters to influence the output:
##' \itemize{
##'  \item{batch}{ a logical indicating whether the second factor was used.}
##'  \item{comparison}{ which comparison was used.}
##'  \item{normMethod}{ the normalization method.}
##'  \item{sigThresh}{ the threshold...}
##'  \item{useSigThresh}{ ...and whether it should be used.}
##' }
##' @param result
##' \itemize{
##'  \item{analysis}{ which analysis was used.}
##'  \item{featureLevel}{ which feature level was used.}
##'  \item{countName}{ which data column was used.}
##'  \item{method}{ which method was used.}
##'  \item{pValue}{ counts the number of features.}
##'  \item{isPresentProbe}{ counts the number of features with counts above threshold.}
##' }
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{addFlexTable}}
addCountResultSummary = function(doc, param, result){
  doc = addTitle(doc, "Result Summary", level=2)
  settings = c("Analysis:"=result$analysis)
  settings = append(settings, c("Feature level:"=result$featureLevel))
  settings = append(settings, c("Data Column Used:"=result$countName))
  settings = append(settings, c("Method:"=result$method))
  if (ezIsSpecified(param$batch)){
    settings = append(settings, c("Statistical Model:"="used provided second factor"))
  }
  settings = append(settings, c("Comparison:"=param$comparison))
  if (!is.null(param$normMethod)){
    settings = append(settings, c("Normalization:"=param$normMethod))
  }
  settings = append(settings, c("Number of features:"=length(result$pValue)))
  if (!is.null(result$isPresentProbe)){
    settings = append(settings, c("Number of features with counts above threshold:"=sum(result$isPresentProbe)))
  }
  if (param$useSigThresh){
    settings = append(settings, c("Log2 signal threshold:"=signif(log2(param$sigThresh), digits=4)))
    settings = append(settings, c("Linear signal threshold:"=signif(param$sigThresh, digits=4)))
  }
  doc = addFlexTable(doc, ezFlexTable(as.data.frame(settings), add.rownames=TRUE))
}

##' @title Adds a result file
##' @description Adds a result file in text format or zipped.
##' @param doc an object of the class bsdoc to add the results to.
##' @param param a list of parameters that pastes the \code{comparison} into the file name and does a zip file if \code{doZip} is true.
##' @param result a list of results.
##' @param rawData a list of raw data. Usually obtained from \code{loadCountDataset()}.
##' @param useInOutput a logical specifying whether to use most of the result information.
##' @param file a character representing the name of the result file.
##' @template roxygen-template
##' @return Returns the name of the result file.
##' @seealso \code{\link{writeTxtLinksToHtml}}
addResultFile = function(doc, param, result, rawData, useInOutput=TRUE,
                           file=paste0("result--", param$comparison, ".txt")){
  seqAnno = rawData$seqAnno
  probes = names(result$pValue)[useInOutput]
  y = data.frame(row.names=probes, stringsAsFactors=FALSE, check.names=FALSE)
  y[ , colnames(seqAnno)] = sapply(seqAnno[match(probes, rownames(seqAnno)), ], as.character)
  y$"log2 Signal" = result$log2Expr[useInOutput]
  y$"isPresent" = result$isPresentProbe[useInOutput]
  y$"log2 Ratio" = result$log2Ratio[useInOutput]
  y$"gfold (log2 Change)" = result$gfold[useInOutput]
  y$"log2 Effect" = result$log2Effect[useInOutput]
  y$"probesetCount" = result$nProbes[useInOutput]
  y$"presentProbesetCount" = result$nPresentProbes[useInOutput]
  y$ratio = result$ratio[useInOutput]
  y$pValue = result$pValue[useInOutput]
  y$fdr = result$fdr[useInOutput]
  for (nm in grep("Tukey pValue", names(result), value=TRUE)){
    y[[nm]] = result[[nm]][useInOutput]
  }
  if (!is.null(result$groupMeans)){
    groupMeans = result$groupMeans[useInOutput, ]
    colnames(groupMeans) = paste("log2 Avg of", colnames(groupMeans))
    y = data.frame(y, groupMeans, check.names=FALSE, stringsAsFactors=FALSE)
  }

  if (!is.null(result$xNorm)){
    yy = result$xNorm[useInOutput, ]
    colnames(yy) = paste(colnames(yy), "[normalized count]")
    y = cbind(y, yy)
  }
  yy = getRpkm(rawData)[useInOutput, ]
  if (!is.null(yy)){
    colnames(yy) = paste(colnames(yy), "[FPKM]")
    y = cbind(y, yy)
  }
  y = y[order(y$fdr, y$pValue), ]
  ezWrite.table(y, file=file, head="Identifier", digits=4)
  if (param$doZip){
    zipLink = zipFile(file)
    if (!is.null(doc)){
      doc = addParagraph(doc, pot(paste("<a href='", zipLink, "' type='application/zip'>", zipLink, "</a>")))
    }
  } else {
    if (!is.null(doc)){
      doc = addParagraph(doc, pot(paste("<a href='", file, "' type='application/txt'>", file, "</a>")))
    }
  }
  return(list(resultFile=file))
}

##' @title Adds QC scatter plots
##' @description Adds QC scatter plots to an html file.
##' @param doc an object of the class bsdoc to add the plots to.
##' @param param a list of parameters. If \code{writeScatterPlots} is false, the function returns NULL.
##' @param design a data.frame containing the factorial design.
##' @param conds a named character vector containing the conditions of the factorial design.
##' @param rawData a list of raw data. Usually obtained from \code{loadCountDataset()}.
##' @param signalCond a set of values containing signals averaged by the conditions.
##' @param isPresentCond Either NULL or a data.frame containing coloring information.
##' @param colors a character vector containing rgb codes. The default is a scale from blue to red.
##' @param types a character vector containing the types.
##' @template roxygen-template
addQcScatterPlots = function(doc, param, design, conds, rawData, signalCond, isPresentCond, types=NULL){
  samples = rownames(design)
  nConds = length(unique(conds))
  signal = getSignal(rawData)
  signal[signal <= 0] = NA
  isPresent = ezPresentFlags(signal, presentFlag=rawData$presentFlag, param=param, isLog=rawData$isLog)
  signalRange = range(signal, na.rm=TRUE)
  doc = addTitle(doc, "Scatter Plots by Conditions", level=2)
  if (!is.null(rawData$seqAnno$gc)){
    gcTypes = data.frame("GC < 0.4"=as.numeric(rawData$seqAnno$gc) < 0.4,
                         "GC > 0.6"=as.numeric(rawData$seqAnno$gc) > 0.6,
                         check.names=FALSE)
  } else {
    gcTypes = NULL
  }
  if (!is.null(rawData$seqAnno$width)){
    widthTypes = data.frame("width < 500nt"=as.numeric(rawData$seqAnno$width) < 500, 
                            "width > 5000nt"=as.numeric(rawData$seqAnno$width) > 5000,
                            check.names=FALSE)
  } else {
    widthTypes = NULL
  }
  if (nConds > 1 & nConds <=  param$allPairsMaxCondNumber){
    plotter = EzPlotterAllPairScatter$new(x=signalCond)
    defLink = ezImageFileLink(plotter, file="allPairs-scatter.png", isPresent=isPresentCond, types=types)
    if (!is.null(gcTypes)){
      plotter = EzPlotterAllPairScatter$new(x=signalCond)
      gcLink = ezImageFileLink(plotter, file="allPairs-scatter-byGc.png", main="color by GC", isPresent=isPresentCond, types=gcTypes) 
    }
    if (!is.null(widthTypes)){
      plotter = EzPlotterAllPairScatter$new(x=signalCond)
      widthLink = ezImageFileLink(plotter, file="allPairs-scatter-byWidth.png", main="color by width", isPresent=isPresentCond, types=widthTypes)
    }
    doc = addFlexTable(doc, ezFlexTable(cbind(defLink, gcLink, widthLink)))
  }
  for (i in 1:min(4, ncol(design))){
    for (cond in unique(design[,i])){
      idx = which(cond == design[,i])
      if (length(idx) > 1){
        idx = idx[order(samples[idx])] ## order alphabetically
        condName = paste(colnames(design)[i], cond)
        doc = addTitle(doc, condName, level=3)
        pngName = ezValidFilename(paste0(condName, "-scatter.png"))
        plotter = EzPlotterScatter(y=signal[ ,idx])
        doc = addParagraph(doc, ezImageFileLink(plotter, file=pngName, isPresent=isPresent[ ,idx], types=types,
                                                lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL))
        if (!is.null(gcTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByGcScatter.png"))
          plotter = EzPlotterScatter(y=signal[ ,idx])
          doc = addParagraph(doc, ezImageFileLink(plotter, file=pngName, isPresent=isPresent[ ,idx], types=gcTypes,
                                                  lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL))
        }
        if (!is.null(widthTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByWidthScatter.png"))
          plotter = EzPlotterScatter(y=signal[ ,idx])
          doc = addParagraph(doc, ezImageFileLink(plotter, file=pngName, isPresent=isPresent[ ,idx], types=widthTypes,
                                                  lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL))
        }
      }
    }
  }
}


addTestScatterPlots = function(doc, param, x, result, seqAnno, types=NULL){
  if (is.null(types)){
    types = data.frame(row.names=rownames(x))
    if ("IsControl" %in% colnames(seqAnno)){
      types$Controls = seqAnno[ rownames(x), "IsControl"]
    }
  }
  msg = "Highlighting significants with: "
  if (!is.null(param$pValueHighlightThresh)){
    significants = result$pValue <= param$pValueHighlightThresh & result$usedInTest
    types$Significants = result$pValue <= param$pValueHighlightThresh & result$usedInTest
    msg = paste(msg, "p <= ", param$pValueHighlightThresh)
    if (!is.null(param$log2RatioHighlightThresh)){
      msg = paste(msg, "and log ratio >= ", param$log2RatioHighlightThresh)
      if (!is.null(result$log2Ratio)){
        types$Significants = types$Significants & abs(result$log2Ratio) >= param$log2RatioHighlightThresh
      } else {
        types$Significants = types$Significants & result$log2Effect >= param$log2RatioHighlightThresh
      }
    }
  }
  doc = addTitle(doc, "Scatter Plots", level=2)
  doc = addParagraph(doc, msg)
  doc = addTitle(doc, "Between-group Comparison", level=3)
  links = character()
  if (ncol(result$groupMeans) == 2 & !is.null(param$sampleGroup) & !is.null(param$refGroup)){
    sampleValues = 2^result$groupMeans[ , param$sampleGroup]
    refValues = 2^result$groupMeans[ , param$refGroup]
    plotter = EzPlotterScatter$new(x=refValues, y=sampleValues)
    links["scatter"] = ezImageFileLink(plotter, file=paste0(param$comparison, "-scatter.png"),
                                  isPresent=result$usedInTest, types=types,
                                  xlab=param$refGroup, ylab=param$sampleGroup)
    plotter = EzPlotterVolcano$new(log2Ratio=result$log2Ratio, pValue=result$pValue)
    links["volcano"] = ezImageFileLink(plotter, file=paste0(param$comparison, "-volcano.png"),
                                  isPresent=result$usedInTest, types=types, main=param$comparison)
    plotter = EzPlotterVolcano$new(log2Ratio=result$log2Ratio, pValue=result$fdr)
    links["volcanoFdr"] = ezImageFileLink(plotter, file=paste0(param$comparison, "-FDR-volcano.png"),
                                  isPresent=result$usedInTest, types=types, main=param$comparison, yType="FDR")
  } else {
    plotter = EzPlotterAllPairScatter$new(x=2^result$groupMeans)
    links["allPair"] = ezImageFileLink(plotter, file=paste0(param$comparison, "-scatter.png"),
                                     isPresent=result$usedInTest, types=types)
  }
  
  myBreaks = seq(0, 1, by=0.002)
  links["pValueHist"] = paste0(param$comparison, "-pValueHist.png")
  png(file=links["pValueHist"], height=400, width=800)
  histUsed = hist(result$pValue[result$usedInTest], breaks=myBreaks, plot=FALSE)
  histAbs = hist(result$pValue[!result$usedInTest], breaks=myBreaks, plot=FALSE)
  xx = rbind(used=histUsed$counts, absent=histAbs$counts)
  xx = shrinkToRange(xx, c(0, max(xx["used", ])))
  #colnames(x) = histUsed$mids
  barplot(xx, space=0, border=NA, col=c("blue", "darkorange"), 
          xlab="p-value", ylab="counts", ylim=c(0, max(xx["used", ])),
          main="p-value histogram")
  abline(h=sum(result$usedInTest)/ncol(xx))
  at = c(0.01, 0.1, 0.25, 0.5, 0.75, 1)
  axis(1, at=at*ncol(xx), labels = at)
  legend("top", c("used", "absent"), col=c("blue", "darkorange"), pch=20, cex=1)
  links["pValueHist"] = as.html(pot(paste('<img src="', links["pValueHist"], '"/>')))
  dev.off()
  doc = addFlexTable(doc, ezFlexTable(rbind(links)))
  
  if (is.null(x)){
    return()
  }
  
  theRange = 2^(range(x, na.rm=TRUE))
  if (!ezIsSpecified(param$batch)){ ## TODO: we no longer use pairing, we now use batch which is more general; however these plots only work if batch is a real pairing
    for (group in unique(c(param$refGroup, colnames(result$groupMeans)))){
      idx = which(group == param$grouping)
      if (length(idx) > 1){
        doc = addTitle(doc, paste("Intra-group Comparison:", group), level=3)
        pngName = paste0(group, "-scatter.png")
        xlab = paste("Avg of", group)
        refValue = result$groupMeans[ , group]
        plotter = EzPlotterScatter$new(x=2^refValues, y=2^x[, idx, drop=FALSE])
        doc = addParagraph(doc, ezImageFileLink(plotter, file=pngName, isPresent=result$isPresent[, idx, drop=FALSE],
                                                types=types, lim=theRange, xlab=xlab))
        if (ncol(result$groupMeans) == 2){
          otherGroup = setdiff(colnames(result$groupMeans), group)
          pngName = paste0(group, "-over-", otherGroup, "-scatter.png")
          xlab = paste("Avg of", otherGroup)
          refValue = result$groupMeans[ , otherGroup]
          doc = addParagraph(doc, ezImageFileLink(plotter, file=pngName, isPresent=result$isPresent[, idx, drop=FALSE],
                                                  types=types, lim=theRange, xlab=xlab))
        }
      }
    }
  } else {
    doc = addTitle(doc, paste("Pairs:", param$sampleGroup, "over", param$refGroup), level=3)
    use = param$grouping %in% c(param$sampleGroup, param$refGroup)
    if (all(table(param$batch[use], param$grouping[use]) == 1)){
      groups = paste(param$grouping, param$batch, sep="--")
      sampleGroups = sort(unique(groups[param$grouping == param$sampleGroup]))
      refGroups = sort(unique(groups[param$grouping == param$refGroup]))
      avgValues = averageColumns(x[ ,use], groups[use], mean)
      avgPresent= averageColumns(x[ ,use], groups[use], function(x){mean(x) > 0.5})
      sampleValues = avgValues[ , sampleGroups, drop=FALSE]
      refValues = avgValues[ , refGroups, drop=FALSE]
      samplePresent = avgPresent[ ,sampleGroups, drop=FALSE]
      refPresent = avgPresent[ , refGroups, drop=FALSE]
      pngName = paste0(param$sampleGroup, "-over-", param$refGroup, "-pairs.png")
      plotter = EzPlotterScatter$new(x=2^refValues, y=2^sampleValues)
      doc = addParagraph(doc, ezImageFileLink(plotter, file=pngName, isPresent=samplePresent | refPresent,
                                              types=types, lim=theRange, xlab=colnames(refValues)))
    }
  }
}

#######################################################################################################################################
#######################################################################################################################################
#####################################################  GO ANALYSIS STUFF  #############################################################
#######################################################################################################################################
#######################################################################################################################################


addGOClusterResult = function(doc, param, clusterResult){
  ontologies = names(clusterResult$GO)
  tables = ezMatrix("", rows=paste("Cluster", 1:clusterResult$nClusters), cols=ontologies)
  for (onto in ontologies){
    for (i in 1:clusterResult$nClusters){
      x = clusterResult$GO[[onto]][[i]]
      tables[i, onto] = goResultToHtmlTable(x, param$pValThreshFisher, param$minCountFisher, onto=onto);
    }
  }
  addTableToReport(tables, doc, border=2,
                   bgcolors=matrix(gsub("FF$", "", clusterResult$clusterColors), nrow=clusterResult$nClusters, ncol=1))
}


addGOTables = function(doc, param, goResult){
  doc = addTitle(doc, "GO Enrichment Analysis", level=3)
  doc = addParagraph(doc, "Red GO categories are overrepresented among the significantly upregulated genes")
  doc = addParagraph(doc, "Blue GO categories are overrepresented among the significantly downregulated genes")
  doc = addParagraph(doc, "Black GO categories are overrepresented among all signifcantly regulated genes")
  doc = addParagraph(doc, paste("Maximum number of terms displayed:", param$maxNumberGroupsDisplayed))
  tables = ezMatrix("", rows="Cats", cols=names(goResult))
  txtFiles = character()
  for (onto in names(goResult)){
    x = goResult[[onto]]
    tables[1,onto] = goResultToHtmlTable2(x, param$pValThreshFisher, 
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
  addTableToReport(tables, doc, border=2, valign="top")
  if (param$doZip){
    addTxtLinksToReport(txtFiles, mime="application/zip", doc=doc)
  } else {
    addTxtLinksToReport(txtFiles, mime="application/txt", doc=doc)
  }
}


#######################################################################################################################################
#######################################################################################################################################
#########################################################  GAGE STUFF  ################################################################
#######################################################################################################################################
#######################################################################################################################################


addGageTables = function(doc, param = NULL, gageResults = NULL) {
  doc = addTitle(doc, "GAGE Enrichment Analysis", level=3)
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
      tbl = addTableToReportWhite(res, doc,
                               bgcolors=matrix(gsub("FF$", "", unique(pathColors)), nrow=length(unique(pathColors)), ncol=1))
      tableRows[[signal]] = cbind(imgrow, tbl)
    }
#     ezWrite("</table>", con=html)
    table = ezFlexTable(rbind(unlist(tableRows)), header = TRUE)
    table = addHeaderRow(table, cbind(paste("Heatmap Plot logRatio Signal for", i), paste(i, "significant pathways")))
    doc = addFlexTable(doc, table)
  }
}



