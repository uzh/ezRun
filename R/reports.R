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
##' @param talign "left", "middle" or "right" specifying the position of text within the cells.
##' @param header.columns a logical indicating whether to use a header for the table.
##' @template addargs-template
##' @templateVar fun FlexTable()
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{FlexTable}}
##' @return Returns an object of the class FlexTable.
##' @examples
##' ezFlexTable(data.frame(a=1:5,b=11:15))
ezFlexTable = function(x, border = 1, valign = "top", talign = "left", header.columns = FALSE,  ...){
  if (!is.data.frame(x) & !is.matrix(x)){
    x = ezFrame(x)
  }
  bodyCells = cellProperties(border.width=border, padding=2, vertical.align=valign)
  bodyPars = parProperties(text.align = talign)
  headerCells = cellProperties(border.width=border, padding=2)
  FlexTable(x, body.cell.props = bodyCells, body.par.props = bodyPars,
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
    imgFilePot = pot(paste("<img src='", file, "' title='", mouseOverText, "'/>"))
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

## currently not possible to use from old report opener:
## writeLines(readLines(ezCSSFile()), con=html) ## no custom css for reporteRs
## TODOP: get addJavaScriptIgvStarter() and addIgvSessionLink() to work.
##' @title Opens an html report
##' @description Opens an html report using \code{bsdoc()} from the ReporteRs package. Also adds some introductory elements.
##' @param title a character specifying the title of the html report.
##' @template dataset-template
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{bsdoc}}
##' @seealso \code{\link[ReporteRs]{writeDoc}}
##' @return Returns an object of the class bsdoc to add further elements.
##' @examples
##' theDoc = openBsdocReport(title="My html report")
##' closeBsdocReport(doc=theDoc, file="example.html")
openBsdocReport = function(title="", dataset=NULL, param=NULL){
  doc = bsdoc(title = title)
  if (!is.null(param)){
    addJavaScriptIgvStarter(htmlFile, param$projectId, doc)
  }
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

##' @describeIn openBsdocReport Adds a java function to start Igv from Jnlp.
addJavaScriptIgvStarter = function(htmlFile, projectId, doc){
  jnlpLines1 = paste('<jnlp spec="6.0+" codebase="http://www.broadinstitute.org/igv/projects/current">',
                     '<information>',
                     '<title>IGV 2.3</title>',
                     '<vendor>The Broad Institute</vendor>',
                     '<homepage href="http://www.broadinstitute.org/igv"/>',
                     '<description>IGV Software</description>',
                     '<description kind="short">IGV</description>',
                     '<icon href="IGV_64.png"/>',
                     '<icon kind="splash" href="IGV_64.png"/>',
                     '<offline-allowed/>',
                     '</information>',
                     '<security>',
                     '<all-permissions/>',
                     '</security>',
                     '<update check="background" />', #check="never" policy="never"/>',
                     '<resources>',
                     '<java version="1.6+" initial-heap-size="256m" max-heap-size="1100m" />',
                     '<jar href="igv.jar" download="eager" main="true"/>',
                     '<jar href="batik-codec__V1.7.jar" download="eager"/>',
                     '<jar href="goby-io-igv__V1.0.jar" download="eager"/>',  
                     '<property name="apple.laf.useScreenMenuBar" value="true"/>',
                     '<property name="com.apple.mrj.application.growbox.intrudes" value="false"/>',
                     '<property name="com.apple.mrj.application.live-resize" value="true"/>',
                     '<property name="com.apple.macos.smallTabs" value="true"/>',
                     '<property name="http.agent" value="IGV"/>',
                     '<property name="development" value="false"/>',
                     '</resources>',
                     '<application-desc  main-class="org.broad.igv.ui.Main">',
                     '<argument>--genomeServerURL=http://fgcz-gstore.uzh.ch/reference/igv_genomes.txt</argument>',
                     paste0('<argument>--dataServerURL=', "http://fgcz-gstore.uzh.ch/list_registries/", projectId, '</argument>'),
                     '<argument>',
                     sep="\n")
  jnlpLines2 = paste("</argument>",
                     '</application-desc>',
                     '</jnlp>',
                     sep="\n")
  javaScript = paste("function startIgvFromJnlp(label, locus){",
                     "var theSession = document.location.href.replace('", htmlFile, "', 'igvSession.xml');",
                     "var igvLink = 'data:application/x-java-jnlp-file;charset=utf-8,';",
                     "igvLink += '", RCurl::curlEscape(jnlpLines1), "';",
                     "igvLink += theSession;",
                     "igvLink += '", RCurl::curlEscape(jnlpLines2), "';",
                     "document.write(label.link(igvLink))",
                     "}")
  doc = addJavascript(doc, text=javaScript)
}

##' @describeIn openBsdocReport Adds the session info, the bootstrap menu, a paragraph showing the finishing time and writes the document. \code{file} must have a .html suffix.
closeBsdocReport = function(doc, file, titles=NULL){
  ezSessionInfo()
  doc = addParagraph(doc, pot("sessionInfo.txt", hyperlink = "sessionInfo.txt"))
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
##' @template doc-template
##' @templateVar object anchored title
##' @param title a character specifying the title.
##' @param level an integer specifying the heading level.
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{addTitle}}
##' @seealso \code{\link[ReporteRs]{addCodeBlock}}
##' @examples
##' theDoc = openBsdocReport(title="My html report")
##' title = "My title"
##' addTitleWithAnchor(theDoc, title)
##' closeBsdocReport(doc=theDoc, file="example.html")
addTitleWithAnchor = function(doc, title, level=1){
  doc = addParagraph(doc, as.html(pot(paste0("<a name='", title, "'></a>"))))
  doc = addTitle(doc, value=title, level=level)
}

##' @title Writes an error report
##' @description Writes an error report to an html file. Also creates the file and closes it.
##' @template htmlFile-template
##' @param param a list of parameters to extract the \code{name} from.
##' @template dataset-template
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

##' @title Adds text links
##' @description Adds texts link to an html file.
##' @template doc-template
##' @templateVar object text links
##' @param txtNames a character representing the file names.
##' @param doZip a logical indicating whether to zip the files.
##' @param mime a character representing the type of the links.
##' @template roxygen-template
##' @examples
##' param = ezParam()
##' htmlFile = "example.html"
##' doc =  openBsdocReport(title="My report")
##' addTxtLinksToReport(doc, "dataset.tsv")
##' closeBsdocReport(doc, htmlFile)
addTxtLinksToReport = function(doc, txtNames, doZip=FALSE, mime=ifelse(doZip, "application/zip", "application/text")){
  for (each in txtNames){
    if (doZip){
      each = zipFile(each)
    }
    doc = addParagraph(doc, pot(paste("<a href='", each, "' type='", mime, "'>", each, "</a>")))
  }
}

## NOTE: perhaps obsolete now with ReporteRs
##' @title Adds a table
##' @description Adds a table to a bsdoc object.
##' @template doc-template
##' @templateVar object table
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
##' @template doc-template
##' @templateVar object table
##' @param param a list of parameters to influence the output:
##' \itemize{
##'  \item{batch}{ a logical indicating whether the second factor was used.}
##'  \item{comparison}{ which comparison was used.}
##'  \item{normMethod}{ the normalization method.}
##'  \item{sigThresh}{ the threshold...}
##'  \item{useSigThresh}{ ...and whether it should be used.}
##' }
##' @template result-template
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{addFlexTable}}
addCountResultSummary = function(doc, param, result){
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

##' @title Adds tables of the significant counts
##' @description Adds tables of the significant counts.
##' @template doc-template
##' @templateVar object table
##' @template result-template
##' @param pThresh a numerical indicating the p-value threshold. genes
##' @template roxygen-template
addSignificantCounts = function(doc, result, pThresh=c(0.1, 0.05, 1/10^(2:5))){
  sigTable = ezFlexTable(getSignificantCountsTable(result, pThresh=pThresh),
                         header.columns = TRUE, add.rownames = TRUE, talign = "right")
  sigFcTable = ezFlexTable(getSignificantFoldChangeCountsTable(result, pThresh=pThresh),
                           header.columns = TRUE, add.rownames = TRUE, talign = "right")
  tbl = ezGrid(cbind(as.html(sigTable), as.html(sigFcTable)))
  doc = addFlexTable(doc, tbl)
  return(doc)
}

##' @describeIn addSignificantCounts Gets the table containing the significant counts.
getSignificantCountsTable = function(result, pThresh=1/10^(1:5), genes=NULL){
  sigTable = ezMatrix(NA, rows=paste("p <", pThresh), cols=c("#significants", "FDR"))
  for (i in 1:length(pThresh)){
    isSig = result$pValue < pThresh[i] & result$usedInTest == 1
    if (is.null(genes)){
      sigTable[i, "#significants"] = sum(isSig, na.rm=TRUE)
    } else {
      sigTable[i, "#significants"] = length(na.omit(unique(genes[isSig])))
    }
    if ( sigTable[i, "#significants"] > 0){
      sigTable[i, "FDR"] = signif(max(result$fdr[isSig], na.rm=TRUE), digits=4)
    }
  }
  sigTable
}

##' @describeIn addSignificantCounts Gets the table containing the significant fold change counts.
getSignificantFoldChangeCountsTable = function(result, pThresh=1/10^(1:5), fcThresh = c(1, 1.5, 2, 3, 4, 8, 10), genes=NULL){
  
  ## counts the significant entries
  ## if genes is given counts the number of different genes that are significant
  if (!is.null(result$log2Ratio)){
    fc = 2^abs(result$log2Ratio)
  } else {
    stopifnot(!is.null(result$log2Effect))
    fc = 2^abs(result$log2Effect)
  }
  
  sigFcTable = ezMatrix(NA, rows=paste("p <", pThresh), cols=paste("fc >=", fcThresh))
  for (i in 1:length(pThresh)){
    for (j in 1:length(fcThresh)){
      isSig = result$pValue < pThresh[i] & result$usedInTest == 1 & fc >= fcThresh[j]
      if (is.null(genes)){
        sigFcTable[i, j] = sum(isSig, na.rm=TRUE)
      } else {
        sigFcTable[i, j] = length(unique(na.omit(genes[isSig])))
      }
    }
  }
  sigFcTable
}

##' @title Adds a result file
##' @description Adds a result file in text format or zipped.
##' @template doc-template
##' @templateVar object result
##' @param param a list of parameters that pastes the \code{comparison} into the file name and does a zip file if \code{doZip} is true.
##' @template result-template
##' @template rawData-template
##' @param useInOutput a logical specifying whether to use most of the result information.
##' @param file a character representing the name of the result file.
##' @template roxygen-template
##' @return Returns the name of the result file.
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
  doc = addTxtLinksToReport(doc, file, param$doZip)
  return(list(resultFile=file))
}
