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
##' @templateVar fun \code{FlexTable()}
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


## how to add help text? for each plot seperately or not?
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

closeBsdocReport = function(doc, file){
  doc = addParagraph(doc, paste("Finished", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  writeDoc(doc, file=file)
}

# writeDoc(html, "bla.html")

