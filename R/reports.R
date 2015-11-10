###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Wrapper for \code{FlexTable()}
##' @description Wraps \code{FlexTable()} with defaults to remove the cell header and cell borders.
##' @param x a matrix or data.frame to turn into an object of the class FlexTable.
##' @param border an integer specifying the width of the table borders.
##' @param valign "bottom", "middle" or "top" specifying the position of table cell contents.
##' @param header.columns a logical indicating whether to use a header for the table.
##' @template addargs-template
##' @templateVar fun FlexTable
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{FlexTable}}
##' @return Returns an object of the class FlexTable.
##' @examples
##' ezFlexTable(data.frame(a=1:5,b=11:15))
ezFlexTable = function(x, border = 1, valign = "top", header.columns = FALSE,  ...){
  if (!is.data.frame(x) & !is.matrix(x)){
    x = ezFrame(x)
  }
  bodyCells = cellProperties(border.width = border, vertical.align=valign)
  headerCells = cellProperties(border.width = border)
  FlexTable(x, body.cell.props = bodyCells,
            header.cell.props = headerCells,
            header.columns = header.columns, ...)
}

##' @describeIn ezFlexTable A flex table without borders.
ezGrid = function(x, header.columns = FALSE,  valign = "top", ...){
  if (!is.data.frame(x) & !is.matrix(x)){
    x = ezFrame(x)
  }
  FlexTable(x, body.cell.props = cellProperties(border.width = 0, vertical.align = valign),
            header.cell.props = cellProperties(border.width = 0),
            header.columns = header.columns, ...)
}

# how to add help text? for each plot seperately or not?
##' @title Gets an image link as html
##' @description Gets an image link as html. Also plots and creates the image.
##' @param plotCmd an expression of plot commands.
##' @param file a character specifying the name of the image with a .png suffix.
##' @param name a character specifying the name of the image together with \code{plotType}, if \code{file} is null.
##' @param plotType a character specifying the name of the image together with \code{name}, if \code{file} is null.
##' @param mouseOverText a character specifying the text being displayed when mousing over the image.
##' @param addPdfLink a logical indicating whether to add a link on the image to a pdf version of itself.
##' @param width an integer specifying the width of each plot to create an image from.
##' @param height an integer specifying the height of each plot to create an image from.
##' @param ppi an integer specifying points per inch.
##' @param envir the environment to evaluate \code{plotCmd} in.
##' @template roxygen-template
##' @return Returns a character specifying a link to an image in html.
##' @examples
##' x = 1:10
##' plotCmd = expression({
##'   plot(x)
##'   text(2,1, "my Text")
##' })
##' ezImageFileLink(plotCmd)
ezImageFileLink = function(plotCmd, file=NULL, name="imagePlot", plotType="plot", mouseOverText="my mouse over",
                           addPdfLink=TRUE, width=480, height=480, ppi=72, envir=parent.frame()){
  if (is.null(file)){
    file = paste0(name, "-", plotType, ".png")
  }
  png(file, width=width, height=height)
  eval(plotCmd, envir = envir)
  dev.off()
  if (addPdfLink) {
    pdfName = sub(".png$", ".pdf", file)
    pdf(file=pdfName, width=width/ppi, height=height/ppi)
    eval(plotCmd, envir=envir)
    dev.off()
    imgFilePot = pot(paste("<img src='", file, "' title='", mouseOverText, "'/>"),
                     hyperlink = pdfName)
  } else {
    imgFilePot = pot(paste("<img src='", pngName, "' title='", mouseOverText, "'/>"))
  }
  return(as.html(imgFilePot))
}

##' @title Pastes image links as html
##' @description A simple wrapper that pastes image links as html.
##' @param image character(s) specifying links to image files.
##' @template roxygen-template
##' @return Returns html link(s).
##' @examples
##' imgLinks("link.png")
imgLinks = function(image){
  links = character()
  for (each in image){
    links[each] = as.html(pot(paste0("<img src='", each, "'/>")))
  }
  return(links)
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
  require(ReporteRs, warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  doc = bsdoc(title = title)
  pot1 = pot(paste("Started on", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "--&#160;"))
  pot2 = as.html(pot("Documentation", hyperlink = "http://fgcz-sushi.uzh.ch/doc/methods-20140422.html"))
  doc = addFlexTable(doc, ezGrid(cbind(pot1, pot2)))
  addTitleWithAnchor(doc, title)
  if (!is.null(dataset)){
    ezWrite.table(dataset, file="dataset.tsv", head="Name")
    doc = addParagraph(doc, pot("dataset.tsv", hyperlink = "dataset.tsv"))
  }
  return(doc)
}

##' @describeIn openBsdocReport Adds a paragraph showing the finishing time and writes the document. \code{file} must have a .html suffix.
closeBsdocReport = function(doc, file, titles=NULL){
  doc = addParagraph(doc, paste("Finished", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  bootStrap = BootstrapMenu("Functional Genomics Center Zurich", link = "http://www.fgcz.ethz.ch")
  if (ezIsSpecified(titles)){
    ddMenu = DropDownMenu("Navigation")
    for (each in titles){
      ddMenu = addLinkItem(ddMenu, label=each, link=paste0("#", each))
    }
    bootStrap = addLinkItem(bootStrap, dd=ddMenu)
  }
  doc = addBootstrapMenu(doc, bootStrap)
  writeDoc(doc, file=file)
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
##' @param doc an object of the class bsdoc to add the table to.
##' @param x a matrix or data.frame to paste a table from.
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
##' ezAddTable(html, x, head="Example", bgcolors="red")
##' closeBsdocReport(html, "example.html")
ezAddTable = function(doc, x, bgcolors=NULL, valign="middle", border=1, head=""){
  bodyCells = cellProperties(border.width=border, vertical.align=valign)
  table = FlexTable(x, header.columns = FALSE, body.cell.props=bodyCells,
                    header.cell.props=cellProperties(border.width = border))
  if (!is.null(bgcolors)){
    table = setFlexTableBackgroundColors(table, j=1:ncol(x), colors=bgcolors)
  }
  table = addHeaderRow(table, colnames(x))
  addFlexTable(doc, table)
}

##' @describeIn ezAddTable Does the same with a white font and returning the table instead of adding it to the document.
ezAddTableWhite = function(x, bgcolors=NULL, valign="middle", border=1, head=""){
  if (is.null(bgcolors)){
    bgcolors = matrix("#ffffff", nrow=nrow(x), ncol=ncol(x))
  }
  ##x = cbind(rownames(x),x)
  x = as.html(pot(paste('<font color="white">', x, '</font>')))
  bodyCells = cellProperties(border.width=border, vertical.align=valign)
  table = FlexTable(x, header.columns = FALSE, body.cell.props=bodyCells,
                    header.cell.props=cellProperties(border.width = border))
  table = setFlexTableBackgroundColors(table, j=1:ncol(x), colors=bgcolors)
  table = addHeaderRow(table, colnames(x))
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
  settings = character()
  settings["Analysis:"] = result$analysis
  settings["Genome Build:"] = param$ezRef@refBuild
  settings["Feature level:"] = result$featureLevel
  settings["Data Column Used:"] = result$countName
  settings["Method:"] = result$method
  if (ezIsSpecified(param$batch)){
    settings["Statistical Model:"] = "used provided second factor"
  }
  settings["Comparison:"] = param$comparison
  if (!is.null(param$normMethod)){
    settings["Normalization:"] = param$normMethod
  }
  settings["Number of features:"] = length(result$pValue)
  if (!is.null(result$isPresentProbe)){
    settings["Number of features with counts above threshold:"] = sum(result$isPresentProbe)
  }
  if (param$useSigThresh){
    settings["Log2 signal threshold:"] = signif(log2(param$sigThresh), digits=4)
    settings["Linear signal threshold:"] = signif(param$sigThresh, digits=4)
  }
  doc = addFlexTable(doc, ezGrid(settings, add.rownames=TRUE))
}


addSignificantCounts = function(doc, result, pThresh=c(0.1, 0.05, 1/10^(2:5))){
  
  sigTable = ezFlexTable(getSignificantCountsTable(result, pThresh=pThresh), header.columns = TRUE, add.rownames = TRUE)
  sigFcTable = ezFlexTable(getSignificantFoldChangeCountsTable(result, pThresh=pThresh), header.columns = TRUE, add.rownames=TRUE)
  doc = addTitle(doc, "Significant Counts", level=3)
  tbl = ezGrid(cbind(as.html(sigTable), as.html(sigFcTable)))
  doc = addFlexTable(doc, tbl)
  return(doc)
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
    plotCmd = expression({
      ezAllPairScatter(signalCond, isPresent=isPresentCond, types=types)
    })
    defLink = ezImageFileLink(plotCmd, file="allPairs-scatter.png",
                              width=min(max(ncol(signalCond) * 200, 480), 2000),
                              height=min(max(ncol(signalCond) * 200, 480), 2000))
    if (!is.null(gcTypes)){
      plotCmd = expression({
        ezAllPairScatter(signalCond, main="color by GC", isPresent=isPresentCond, types=gcTypes)
      })
      gcLink = ezImageFileLink(plotCmd, file="allPairs-scatter-byGc.png",
                               width=min(max(ncol(signalCond) * 200, 480), 2000),
                               height=min(max(ncol(signalCond) * 200, 480), 2000))
    }
    if (!is.null(widthTypes)){
      plotCmd = expression({
        ezAllPairScatter(signalCond, main="color by width", isPresent=isPresentCond, types=widthTypes)
      })
      widthLink = ezImageFileLink(plotCmd, file="allPairs-scatter-byWidth.png",
                                  width=min(max(ncol(signalCond) * 200, 480), 2000),
                                  height=min(max(ncol(signalCond) * 200, 480), 2000))
    }
    doc = addFlexTable(doc, ezGrid(cbind(defLink, gcLink, widthLink)))
  }
  for (i in 1:min(4, ncol(design))){
    for (cond in unique(design[,i])){
      idx = which(cond == design[,i])
      if (length(idx) > 1){
        idx = idx[order(samples[idx])] ## order alphabetically
        condName = paste(colnames(design)[i], cond)
        doc = addTitle(doc, condName, level=3)
        pngName = ezValidFilename(paste0(condName, "-scatter.png"))
        plotCmd = expression({
          ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=types, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
        })
        doc = addParagraph(doc, ezImageFileLink(plotCmd, file=pngName,
                                                width=min(ncol(as.matrix(signal)), 6) * 480,
                                                height=ceiling(ncol(as.matrix(signal))/6) * 480))
        if (!is.null(gcTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByGcScatter.png"))
          plotCmd = expression({
            ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=gcTypes, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
          })
          doc = addParagraph(doc, ezImageFileLink(plotCmd, file=pngName,
                                                  width=min(ncol(as.matrix(signal)), 6) * 480,
                                                  height=ceiling(ncol(as.matrix(signal))/6) * 480))
        }
        if (!is.null(widthTypes)){
          pngName = ezValidFilename(paste0(condName, "-ByWidthScatter.png"))
          plotCmd = expression({
            ezScatter(y=signal[ ,idx], isPresent=isPresent[ ,idx], types=widthTypes, lim=signalRange, xlab=paste("Avg of", cond), ylab=NULL)
          })
          doc = addParagraph(doc, ezImageFileLink(plotCmd, file=pngName,
                                                  width=min(ncol(as.matrix(signal)), 6) * 480,
                                                  height=ceiling(ncol(as.matrix(signal))/6) * 480))
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
    plotCmd = expression({
      ezScatter(x=refValues, y=sampleValues, isPresent=result$usedInTest, types=types, xlab=param$refGroup, ylab=param$sampleGroup)
    })
    links["scatter"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-scatter.png"),
                                       width=min(ncol(as.matrix(sampleValues)), 6) * 480,
                                       height=ceiling(ncol(as.matrix(sampleValues))/6) * 480)
    plotCmd = expression({
      ezVolcano(log2Ratio=result$log2Ratio, pValue=result$pValue, isPresent=result$usedInTest, types=types, main=param$comparison)
    })
    links["volcano"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-volcano.png"))
    plotCmd = expression({
      ezVolcano(log2Ratio=result$log2Ratio, pValue=result$fdr, isPresent=result$usedInTest, types=types, main=param$comparison, yType="FDR")
    })
    links["volcanoFdr"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-FDR-volcano.png"))
  } else {
    plotCmd = expression({
      ezAllPairScatter(x=2^result$groupMeans, isPresent=result$usedInTest, types=types)
    })
    links["allPair"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-scatter.png"),
                                       width=min(max(ncol(result$groupMeans) * 200, 480), 2000),
                                       height=min(max(ncol(result$groupMeans) * 200, 480), 2000))
  }
  
  plotCmd = expression({
    myBreaks = seq(0, 1, by=0.002)
    histUsed = hist(result$pValue[result$usedInTest], breaks=myBreaks, plot=FALSE)
    histAbs = hist(result$pValue[!result$usedInTest], breaks=myBreaks, plot=FALSE)
    xx = rbind(used=histUsed$counts, absent=histAbs$counts)
    xx = shrinkToRange(xx, c(0, max(xx["used", ])))
    barplot(xx, space=0, border=NA, col=c("blue", "darkorange"), 
            xlab="p-value", ylab="counts", ylim=c(0, max(xx["used", ])),
            main="p-value histogram")
    abline(h=sum(result$usedInTest)/ncol(xx))
    at = c(0.01, 0.1, 0.25, 0.5, 0.75, 1)
    axis(1, at=at*ncol(xx), labels = at)
    legend("top", c("used", "absent"), col=c("blue", "darkorange"), pch=20, cex=1)
  })
  links["pValueHist"] = ezImageFileLink(plotCmd, file=paste0(param$comparison, "-pValueHist.png"), height=400, width=800)
  doc = addFlexTable(doc, ezGrid(rbind(links)))
  
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
        refValues = result$groupMeans[ , group]
        plotCmd = expression({
          ezScatter(x=2^refValues, y=2^x[, idx, drop=FALSE], isPresent=result$isPresent[, idx, drop=FALSE], types=types, lim=theRange, xlab=xlab)
        })
        doc = addParagraph(doc, ezImageFileLink(plotCmd, file=pngName,
                                                width=min(ncol(as.matrix(x[, idx, drop=FALSE])), 6) * 480,
                                                height=ceiling(ncol(as.matrix(x[, idx, drop=FALSE]))/6) * 480))
        if (ncol(result$groupMeans) == 2){
          otherGroup = setdiff(colnames(result$groupMeans), group)
          pngName = paste0(group, "-over-", otherGroup, "-scatter.png")
          xlab = paste("Avg of", otherGroup)
          refValues = result$groupMeans[ , otherGroup]
          plotCmd = expression({
            ezScatter(x=2^refValues, y=2^x[, idx, drop=FALSE], isPresent=result$isPresent[, idx, drop=FALSE], types=types, lim=theRange, xlab=xlab)
          })
          doc = addParagraph(doc, ezImageFileLink(plotCmd, file=pngName,
                                                  width=min(ncol(as.matrix(x[, idx, drop=FALSE])), 6) * 480,
                                                  height=ceiling(ncol(as.matrix(x[, idx, drop=FALSE]))/6) * 480))
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
      plotCmd = expression({
        ezScatter(x=2^refValues, y=2^sampleValues, isPresent=samplePresent | refPresent, types=types, lim=theRange, xlab=colnames(refValues))
      })
      doc = addParagraph(doc, ezImageFileLink(plotCmd, file=pngName,
                                              width=min(ncol(as.matrix(sampleValues)), 6) * 480,
                                              height=ceiling(ncol(as.matrix(sampleValues))/6) * 480))
    }
  }
}

#######################################################################################################################################
#######################################################################################################################################
#####################################################  GO ANALYSIS STUFF  #############################################################
#######################################################################################################################################
#######################################################################################################################################


# addGOClusterResult = function(doc, param, clusterResult){
#   ontologies = names(clusterResult$GO)
#   tables = ezMatrix("", rows=paste("Cluster", 1:clusterResult$nClusters), cols=ontologies)
#   for (onto in ontologies){
#     for (i in 1:clusterResult$nClusters){
#       x = clusterResult$GO[[onto]][[i]]
#       tables[i, onto] = goResultToHtmlTable(x, param$pValThreshFisher, param$minCountFisher, onto=onto);
#     }
#   }
#   ezAddTable(doc, tables, border=2,
#                    bgcolors=matrix(gsub("FF$", "", clusterResult$clusterColors), nrow=clusterResult$nClusters, ncol=ncol(tables)))
# }


## tables within a table within a table works flawlessly with this code, I never experienced the <br></br> problem
## perhaps tables in the functions goClusterTable() and goUpDownTables() should be done as a list with each flex table being passed as.html()
# param = ezParam(userParam = list('refBuild' = 'Schizosaccharomyces_pombe/Ensembl/EF2/Annotation/Version-2013-03-07'))
# clusterResult = list("GO"=list("A"=letters[1:3], "B"=letters[4:6], "C"=letters[7:9]), nClusters=3, clusterColors=rainbow(3))
# names(clusterResult$GO) = c("CC", "BP", "MF")
# doc = openBsdocReport("My title")
# # goLink = as.html(ezGrid(c("Background color corresponds to the row colors in the heatmap plot.",
# #                           as.html(goClusterTable(param, clusterResult)))))
# #tbl = ezGrid(ezFrame("Cluster Plot"=imgLinks("Rplot001.png"), "GO categories of feature clusters"=goLink), header.columns = TRUE)
# innerTables = list()
# for (i in seq(2,24,2)){
#   innerTables[[i]] = as.html(ezFlexTable(matrix(rep(i, i), ifelse(i %% 4 == 0, 4, 2)), header.columns = TRUE))
# }
# table1 = ezMatrix(unlist(innerTables), rows=1:3, cols=letters[1:4])
# tableLink = as.html(ezFlexTable(table1, header.columns = TRUE))
# file.remove("Rplot001.png")
# png()
# plot(10:1)
# dev.off()
# tbl = ezGrid(ezFrame("Cluster Plot"=imgLinks("Rplot001.png"), "tables"=tableLink), header.columns = TRUE)
# doc = addFlexTable(doc, tbl)
# closeBsdocReport(doc, "example.html")


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


addGoUpDownResult = function(doc, param, goResult){
  
  udt = goUpDownTables(param, goResult)
  
  doc = addTitle(doc, "GO Enrichment Analysis", level=3)
  doc = addParagraph(doc, "Red GO categories are overrepresented among the significantly upregulated genes")
  doc = addParagraph(doc, "Blue GO categories are overrepresented among the significantly downregulated genes")
  doc = addParagraph(doc, "Black GO categories are overrepresented among all signifcantly regulated genes")
  doc = addParagraph(doc, paste("Maximum number of terms displayed:", param$maxNumberGroupsDisplayed))
  doc = addFlexTable(doc, udt$flexTable)
  
  if (param$doZip){
    addTxtLinksToReport(udt$txtFiles, mime="application/zip", doc=doc)
  } else {
    addTxtLinksToReport(udt$txtFiles, mime="application/txt", doc=doc)
  }  
  return(doc) #
}

goUpDownTables = function(param, goResult){
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
  ft = ezFlexTable(tables, border = 2, header.columns = TRUE)
  return(list(flexTable=ft, txtFiles=txtFiles))
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
  ezAddTable(doc, tables, border=2, valign="top")
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
      tbl = ezAddTableWhite(res, bgcolors=matrix(gsub("FF$", "", unique(pathColors)), nrow=length(unique(pathColors)), ncol=1))
      tableRows[[signal]] = cbind(imgrow, tbl)
    }
    #     ezWrite("</table>", con=html)
    table = ezGrid(rbind(unlist(tableRows)), header.columns=TRUE)
    table = addHeaderRow(table, cbind(paste("Heatmap Plot logRatio Signal for", i), paste(i, "significant pathways")))
    doc = addFlexTable(doc, table)
  }
}
