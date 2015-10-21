###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezImageTable = function(x, header=FALSE, ...) {
  FlexTable(x, body.cell.props = cellProperties(border.width = 0),
            header.cell.props = cellProperties(border.width = 0),
            header.columns = header, ...)
}

## how to add help text? for each plot seperately or not?
ezImageFileLink = function(ezPlotter, file=NULL, mouseOverText=ezPlotter$mouseOverText, helpText=ezPlotter$helpText,
                           addPdfLink=TRUE, width=480, height=480, ...) {
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
