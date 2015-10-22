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
##' @param file a character specifying the name of the image.
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
##' closeBsdocReport(doc=theDoc, file="example.html")
openBsdocReport = function(title="", dataset=NULL){
  html = bsdoc(title = title)
  file.copy(ezBannerFile(), ".")
  html = addImage(html, basename(ezBannerFile()))
  html = addParagraph(html, pot("Functional Genomics Center Zurich", hyperlink = "http://www.fgcz.ethz.ch"))
  pot1 = pot(paste("Started on", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "--&#160;"))
  pot2 = as.html(pot("Documentation", hyperlink = "http://fgcz-sushi.uzh.ch/doc/methods-20140422.html"))
  html = addFlexTable(html, ezFlexTable(cbind(pot1, pot2)))
  html = addTitle(html, title)
  if (!is.null(dataset)){
    ezWrite.table(dataset, file="dataset.tsv", head="Name")
    html = addParagraph(html, pot("dataset.tsv", hyperlink = "dataset.tsv"))
  }
  return(html)
}

##' @describeIn openBsdocReport Adds a paragraph showing the finishing time and writes the document.
closeBsdocReport = function(doc, file){
  doc = addParagraph(doc, paste("Finished", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  writeDoc(doc, file=file)
}
