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
##' ezFlexTable(data.frame(a=1:5, b=11:15))
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

# how to add help text? for each plot separately or not?
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
  eval(plotCmd, envir=envir)
  dev.off()
  imgFilePot = paste0("<span><img src='", file, "' title='", mouseOverText, "'/></span>")
  if (addPdfLink) {
    pdfName = sub(".png$", ".pdf", file)
    pdf(file=pdfName, width=width/ppi, height=height/ppi)
    eval(plotCmd, envir=envir)
    dev.off()
    imgFilePot = paste0("<a href='", pdfName, "'>", imgFilePot, "</a>")
  }
  return(imgFilePot)
}

## currently not possible to use from old report opener:
##' @title Opens an html report
##' @description Opens an html report using \code{bsdoc()} from the ReporteRs package. Also adds some introductory elements.
##' @param title a character specifying the title of the html report.
##' @param doc an object of the class bsdoc to write and close.
##' @param file a character specifying the path to write the report in.
##' @param titles a character vector containing the titles of the report. Used to create the bootstrapmenu.
##' @template roxygen-template
##' @seealso \code{\link[ReporteRs]{bsdoc}}
##' @seealso \code{\link[ReporteRs]{writeDoc}}
##' @return Returns an object of the class bsdoc to add further elements.
##' @examples
##' theDoc = openBsdocReport(title="My html report")
##' closeBsdocReport(doc=theDoc, file="example.html")
openBsdocReport = function(title=""){
  require(ReporteRs)
  doc = bsdoc(title = title)
  pot1 = paste("Started on", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "--&#160;")
  pot2 = as.html(pot("Documentation", hyperlink = "http://fgcz-sushi.uzh.ch/doc/methods-20140422.html"))
  addFlexTable(doc, ezGrid(cbind(pot1, pot2)))
  addTitle(doc, title, id=title)
  return(doc)
}

##' @describeIn openBsdocReport Adds the session info, the bootstrap menu, a paragraph showing the finishing time and writes the document. \code{file} must have a .html suffix.
closeBsdocReport = function(doc, file, titles=NULL){
  ezSessionInfo()
  addParagraph(doc, pot("sessionInfo.txt", hyperlink = "sessionInfo.txt"))
  addParagraph(doc, paste("Finished", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  bootStrap = BootstrapMenu("Functional Genomics Center Zurich", link = "http://www.fgcz.ethz.ch")
  if (ezIsSpecified(titles)){
    ddMenu = DropDownMenu("Navigation")
    for (each in titles){
      ddMenu = addLinkItem(ddMenu, label=each, link=paste0("#", each))
    }
    bootStrap = addLinkItem(bootStrap, dd=ddMenu)
  }
  addBootstrapMenu(doc, bootStrap)
  writeDoc(doc, file=file)
}

##' @title Adds a dataset
##' @description Adds a dataset to a bsdoc object and the reference build.
##' @template doc-template
##' @templateVar object dataset
##' @template dataset-template
##' @param param a list of parameters to extract \code{refBuild} from.
##' @template roxygen-template
addDataset = function(doc, dataset, param){
  ezWrite.table(dataset, file="input_dataset.tsv", head="Name")
  tableLink = "InputDataset.html"
  ezInteractiveTable(dataset, tableLink=tableLink, title="Input Dataset")
  addParagraph(doc, ezLink(tableLink, "Input Dataset", target="_blank"))
}

### -----------------------------------------------------------------
### create error report with Rmd
###
writeErrorReport <- function(htmlFile, param=param, error="Unknown Error"){
  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "ErrorReport.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="ErrorReport.Rmd", envir=new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE)
}

##' @title Adds text links
##' @description Adds texts link to a bsdoc object.
##' @template doc-template
##' @templateVar object text links
##' @param txtNames a character vector representing the link names.
##' @param linkName a single character representing the link (with .html).
##' @param txtName a single character representing the link name.
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
    addParagraph(doc, ezLink(each, type=mime))
    #addParagraph(doc, pot(paste("<a href='", each, "' type='", mime, "'>", each, "</a>")))
  }
}

##' @describeIn addTxtLinksToReport Gets the link, its name and returns an html link that will open new windows/tabs.
newWindowLink = function(linkName, txtName=NULL){
  .Deprecated("use ezLink")
  if (is.null(txtName)){
    title = sub(".html", "", linkName)
  } else {
    title = txtName
  }
  jsCall = paste0('popup({linkName: "', linkName, '"});')
  return(pot(paste0("<a href='javascript:void(0)' onClick='", jsCall, "'>", title, "</a>")))
}

ezLink = function(link, label=link, target="", type=""){
  linkTag = paste0("<a href='", link, "'")
  if (target != ""){
    linkTag = paste0(linkTag, " target='", target, "'")
  }
  if (type != ""){
    linkTag = paste0(linkTag, " type='", type, "'")
  }  
  linkTag = paste0(linkTag, ">")
  paste0(linkTag, label, "</a>")
}

##' @title Adds a summary of the count result
##' @description Adds a summary of the count result to a bsdoc object.
##' @template doc-template
##' @templateVar object table
##' @param param a list of parameters to influence the output:
##' \itemize{
##'  \item{grouping2}{ indicates whether a second factor was used.}
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
  settings["Feature level:"] = result$featureLevel
  settings["Data Column Used:"] = result$countName
  settings["Method:"] = result$method
  if (ezIsSpecified(param$grouping2)){
    settings["Statistical Model:"] = "used provided second factor"
  }
  settings["Comparison:"] = param$comparison
  if (!is.null(param$normMethod)){
    settings["Normalization:"] = param$normMethod
  }
  if(!is.null(param$deTest)){
    settings["Differential expression test:"] <- param$deTest
  }
  if (param$useSigThresh){
    settings["Log2 signal threshold:"] = signif(log2(param$sigThresh), digits=4)
    settings["Linear signal threshold:"] = signif(param$sigThresh, digits=4)
  }
  addFlexTable(doc, ezGrid(settings, add.rownames=TRUE))
}

addCountResultSummarySE = function(doc, param, se){
  settings = character()
  settings["Analysis:"] = metadata(se)$analysis
  settings["Feature level:"] = metadata(se)$featureLeve
  settings["Data Column Used:"] = metadata(se)$countName
  settings["Method:"] = metadata(se)$method
  if (ezIsSpecified(param$grouping2)){
    settings["Statistical Model:"] = "used provided second factor"
  }
  settings["Comparison:"] = param$comparison
  if (!is.null(param$normMethod)){
    settings["Normalization:"] = param$normMethod
  }
  if(!is.null(param$deTest)){
    settings["Differential expression test:"] <- param$deTest
  }
  if (param$useSigThresh){
    settings["Log2 signal threshold:"] = signif(log2(param$sigThresh), digits=4)
    settings["Linear signal threshold:"] = signif(param$sigThresh, digits=4)
  }
  addFlexTable(doc, ezGrid(settings, add.rownames=TRUE))
}

makeCountResultSummary = function(param, se){
  settings = character()
  settings["Analysis:"] = metadata(se)$analysis
  settings["Feature level:"] = metadata(se)$featureLevel
  settings["Data Column Used:"] = metadata(se)$countName
  settings["Method:"] = metadata(se)$method
  if (ezIsSpecified(param$grouping2)){
    settings["Statistical Model:"] = "used provided second factor"
  }
  settings["Comparison:"] = param$comparison
  if (!is.null(param$normMethod)){
    settings["Normalization:"] = param$normMethod
  }
  if(!is.null(param$deTest)){
    settings["Differential expression test:"] <- param$deTest
  }
  if (param$useSigThresh){
    settings["Log2 signal threshold:"] = signif(log2(param$sigThresh), digits=4)
    settings["Linear signal threshold:"] = signif(param$sigThresh, digits=4)
  }
  return(as.data.frame(settings))
}

makeSignificantCounts = function(se, pThresh=c(0.1, 0.05, 1/10^(2:5))){
  sigTable = getSignificantCountsTableSE(se, pThresh=pThresh)
  sigFcTable = getSignificantFoldChangeCountsTableSE(se, pThresh=pThresh)
  return(as.data.frame(cbind(sigTable, sigFcTable)))
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

getSignificantCountsTableSE = function(se, pThresh=1/10^(1:5), genes=NULL){
  sigTable = ezMatrix(NA, rows=paste("p <", pThresh), cols=c("#significants", "FDR"))
  for (i in 1:length(pThresh)){
    isSig = rowData(se)$pValue < pThresh[i] & rowData(se)$usedInTest == 1
    if (is.null(genes)){
      sigTable[i, "#significants"] = sum(isSig, na.rm=TRUE)
    } else {
      sigTable[i, "#significants"] = length(na.omit(unique(genes[isSig])))
    }
    if ( sigTable[i, "#significants"] > 0){
      sigTable[i, "FDR"] = signif(max(rowData(se)$fdr[isSig], na.rm=TRUE), 
                                  digits=4)
    }
  }
  sigTable
}

##' @describeIn addSignificantCounts Gets the table containing the significant fold change counts.
getSignificantFoldChangeCountsTable = function(result, pThresh=1/10^(1:5), 
                                               fcThresh = c(1, 1.5, 2, 3, 4, 8, 10), 
                                               genes=NULL){
  
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

getSignificantFoldChangeCountsTableSE = function(se, pThresh=1/10^(1:5), 
                                                 fcThresh = c(1, 1.5, 2, 3, 4, 8, 10), 
                                                 genes=NULL){
  
  ## counts the significant entries
  ## if genes is given counts the number of different genes that are significant
  if (!is.null(rowData(se)$log2Ratio)){
    fc = 2^abs(rowData(se)$log2Ratio)
  } else {
    stopifnot(!is.null(rowData(se)$log2Effect))
    fc = 2^abs(rowData(se)$log2Effect)
  }
  
  sigFcTable = ezMatrix(NA, rows=paste("p <", pThresh), cols=paste("fc >=", fcThresh))
  for (i in 1:length(pThresh)){
    for (j in 1:length(fcThresh)){
      isSig = rowData(se)$pValue < pThresh[i] & rowData(se)$usedInTest == 1 & 
        fc >= fcThresh[j]
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
  if (!is.null(y$featWidth)){
    y$featWidth = as.integer(y$featWidth)
  }
  if (!is.null(y$gc)){
    y$gc = as.numeric(y$gc)
  }
  ezWrite.table(y, file=file, head="Identifier", digits=4)
  addParagraph(doc, paste("Full result table for opening with a spreadsheet program (e.g. Excel: when",
                          "opening with Excel, make sure that the Gene symbols are loaded into a",
                          "column formatted as 'text' that prevents conversion of the symbols to dates):"))
  addTxtLinksToReport(doc, file, param$doZip)
  useInInteractiveTable = c("gene_name", "type", "description", "width", "gc", "isPresent", "log2 Ratio", "pValue", "fdr")
  useInInteractiveTable = intersect(useInInteractiveTable, colnames(y))
  tableLink = sub(".txt", "-viewTopSignificantGenes.html", file)
  ezInteractiveTable(head(y[, useInInteractiveTable, drop=FALSE], param$maxTableRows), tableLink=tableLink, digits=3,
                     title=paste("Showing the", param$maxTableRows, "most significant genes"))
  return(list(resultFile=file))
}

addResultFileSE = function(doc, param, se, useInOutput=TRUE,
                           file=paste0("result--", param$comparison, ".txt")){
  
  se <- se[useInOutput, ]
  y = data.frame(rowData(se), row.names=rownames(se),
                 stringsAsFactors=FALSE, check.names=FALSE)
  y$"isPresent" = y$isPresentProbe
  y$isPresentProbe <- NULL
  y$"log2 Ratio" = y$log2Ratio
  y$log2Ratio <- NULL
  y$"gfold (log2 Change)" = y$gfold
  y$gfold <- NULL
  y$usedInTest <- NULL ## don't output usedInTest.
  
  # We don't export this groupMeans to result file
  #if (!is.null(result$groupMeans)){
  #  groupMeans = result$groupMeans[useInOutput, ]
  #  colnames(groupMeans) = paste("log2 Avg of", colnames(groupMeans))
  #  y = data.frame(y, groupMeans, check.names=FALSE, stringsAsFactors=FALSE)
  #}
  
  if (!is.null(assays(se)$xNorm)){
    yy = assays(se)$xNorm
    colnames(yy) = paste(colnames(yy), "[normalized count]")
    y = cbind(y, yy)
  }
  yy = getRpkm(se)
  if (!is.null(yy)){
    colnames(yy) = paste(colnames(yy), "[FPKM]")
    y = cbind(y, yy)
  }
  y = y[order(y$fdr, y$pValue), ]
   if (!is.null(y$featWidth)){
     ## This is to round the with after averaging the transcript lengths
     y$featWidth = as.integer(y$featWidth)
  }

  ezWrite.table(y, file=file, head="Identifier", digits=4)
  addParagraph(doc, paste("Full result table for opening with a spreadsheet program (e.g. Excel: when",
                          "opening with Excel, make sure that the Gene symbols are loaded into a",
                          "column formatted as 'text' that prevents conversion of the symbols to dates):"))
  addTxtLinksToReport(doc, file, param$doZip)
  useInInteractiveTable = c("gene_name", "type", "description", "width", "gc", "isPresent", "log2 Ratio", "pValue", "fdr")
  useInInteractiveTable = intersect(useInInteractiveTable, colnames(y))
  tableLink = sub(".txt", "-viewTopSignificantGenes.html", file)
  ezInteractiveTable(head(y[, useInInteractiveTable, drop=FALSE], param$maxTableRows), tableLink=tableLink, digits=3,
                     title=paste("Showing the", param$maxTableRows, "most significant genes"))
  return(list(resultFile=file))
}

makeResultFile = function(param, se, useInOutput=TRUE,
                          file=paste0("result--", param$comparison, ".txt")){
  require(tools)
  require(DT, quietly=TRUE)
  se <- se[useInOutput, ]
  y = data.frame(rowData(se), row.names=rownames(se),
                 stringsAsFactors=FALSE, check.names=FALSE)
  y <- cbind(y, as.data.frame(granges(rowRanges(se))))
  y$"isPresent" = y$isPresentProbe
  y$isPresentProbe <- NULL
  y$"log2 Ratio" = y$log2Ratio
  y$log2Ratio <- NULL
  y$"gfold (log2 Change)" = y$gfold
  y$gfold <- NULL
  y$usedInTest <- NULL ## don't output usedInTest.
  
  if (!is.null(assays(se)$xNorm)){
    yy = assays(se)$xNorm
    colnames(yy) = paste(colnames(yy), "[normalized count]")
    y = cbind(y, yy)
  }
  yy = getRpkm(se)
  if (!is.null(yy)){
    colnames(yy) = paste(colnames(yy), "[FPKM]")
    y = cbind(y, yy)
  }
  y = y[order(y$fdr, y$pValue), ]
  if (!is.null(y$featWidth)){
    ### This is to round the with after averaging the transcript lengths
    y$featWidth = as.integer(y$featWidth)
  }
  
  ezWrite.table(y, file=file, digits=4, row.names=FALSE)
  ans <- list()
  ans$resultFile <- file
  if(isTRUE(param$doZip)){
    zipFile <- sub(file_ext(file), "zip", file)
    zip(zipfile=zipFile, files=file)
    ans$resultZip <- zipFile
  }
  
  ## Interactive gene tables
  useInInteractiveTable = c("gene_name", "type", "description", "featWidth", "gc", 
                            "isPresent", "log2 Ratio", "pValue", "fdr")
  useInInteractiveTable = intersect(useInInteractiveTable, colnames(y))
  tableLink = sub(".txt", "-viewTopSignificantGenes.html", file)
  tableDT <- ezInteractiveTableRmd(head(y[, useInInteractiveTable, drop=FALSE], 
                                        param$maxTableRows),
                                   digits=3,
                                   title=paste("Showing the", param$maxTableRows, "most significant genes"))
  DT::saveWidget(tableDT, tableLink)
  ans$resultHtml <- tableLink
  return(ans)
}


makeWebgestaltFiles <- function(param, resultFile){
    fileNames = list()
    result = ezRead.table(resultFile)
    setwdNew('Webgestalt')
    comparison = sub('.txt', '', sub('result--', '', basename(resultFile)))
    background = rownames(result[result$isPresent, ])
    ezWrite.table(background, paste0('ORA_BG_Webgestalt_',comparison, '.txt'), row.names = FALSE, col.names = FALSE)
    
    GSEA = cbind(rownames(result[result$isPresent, ]), result[result$isPresent, 'log2 Ratio'])
    ezWrite.table(GSEA, paste0('GSEA_Input_log2FC_Webgestalt_',comparison, '.rnk'), row.names = FALSE, col.names = FALSE)
    
    GSEA_pVal = cbind(rownames(result[result$isPresent, ]), sign(result[result$isPresent, 'log2 Ratio']) * -log10(result[result$isPresent, 'pValue']))
    ezWrite.table(GSEA_pVal, paste0('GSEA_Input_pVal_Webgestalt_',comparison, '.rnk'), row.names = FALSE, col.names = FALSE)
    
    ORA_Up = rownames(result[result$isPresent & result$pValue < param[['pValueHighlightThresh']] & result[['log2 Ratio']] >= param[['log2RatioHighlightThresh']], ])
    if(length(ORA_Up) > 0){
        ezWrite.table(ORA_Up, paste0('ORA_Up_Webgestalt_', comparison, '.txt') ,row.names = FALSE, col.names = FALSE)
    }
    
    ORA_Both = rownames(result[result$isPresent & result$pValue < param[['pValueHighlightThresh']] & abs(result[['log2 Ratio']]) >= param[['log2RatioHighlightThresh']], ])
    if(length(ORA_Both) > 0){
        ezWrite.table(ORA_Both, paste0('ORA_Both_Webgestalt_', comparison, '.txt') ,row.names = FALSE, col.names = FALSE)
    }
    
    ORA_Down = rownames(result[result$isPresent & result$pValue < param[['pValueHighlightThresh']] & result[['log2 Ratio']] <= (-1*param[['log2RatioHighlightThresh']]), ])
    if(length(ORA_Down) > 0){
        ezWrite.table(ORA_Down, paste0('ORA_Down_Webgestalt_', comparison, '.txt') ,row.names = FALSE, col.names = FALSE)
    }
    setwd('..')
    return('success')
}

runWebgestaltGSEA <- function(param, rnkFile){
  require(WebGestaltR)
  outputDirectory <- file.path(getwd(), 'Webgestalt/GSEA_Results')
  system(paste('mkdir -p', outputDirectory))
  speciesName = limma::strsplit2(param$ezRef['refBuild'],'/')[1]
  organism = paste0(tolower(substr(speciesName,1,1)), sub('.*_', '' ,speciesName))
  myEnrichDatabases = c('geneontology_Biological_Process_noRedundant',
                        'geneontology_Cellular_Component_noRedundant',
                        'geneontology_Molecular_Function_noRedundant',
                        'pathway_KEGG',
                        'pathway_Panther',
                        'pathway_Reactome',
                        'pathway_Wikipathway')
  myEnrichDatabases = intersect(listGeneSet(organism = organism)$name, myEnrichDatabases)
  if(length(intersect(organism, listOrganism())) == 1){
    for (i in 1:length(myEnrichDatabases)){
      projectName = paste(myEnrichDatabases[i], sub('^GSEA_Input_', '', sub('.rnk', '', basename(rnkFile))), sep = '_')
      enrichResult <- WebGestaltR(enrichMethod = "GSEA", organism = organism,
                                  enrichDatabase = myEnrichDatabases[i], interestGeneFile = rnkFile,
                                  interestGeneType = "ensembl_gene_id", collapseMethod = "mean",reportNum = 30, fdrThr = 0.1, topThr = 30,
                                  isOutput = TRUE, outputDirectory = outputDirectory, projectName = projectName, nThreads = param$cores)
    } 
  }
  return('success')
}

############################################################
### probably not needed:

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

##' @title Adds a java function
##' @description Adds a java function to start Igv from Jnlp.
##' @template htmlFile-template
##' @param projectId a character representing the project ID.
##' @template doc-template
##' @templateVar object JS starter
##' @template roxygen-template
addJavaScriptIgvStarter = function(htmlFile, projectId, doc){
  jnlpLines1 = paste('<jnlp spec="6.0+" codebase="http://data.broadinstitute.org/igv/projects/current">',
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
  addJavascript(doc, text=javaScript)
}
