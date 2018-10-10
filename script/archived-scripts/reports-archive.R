##' @title Adds a table
##' @description Adds a table to a bsdoc object.
##' @template doc-template
##' @templateVar object table
##' @param x a matrix or data.frame to paste a table from.
##' @param bgcolors a matrix specifying the background colors.
##' @param valign a character specifying where to align the table elements vertically. Use either "top", "middle" or "bottom".
##' @param border an integer specifying the border width.
##' @param head a character specifying the contents of the upper-left corner of the table.
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

# ## NOTEP: used once in gage-reports.R, but that needs to be refactored anyway.
# ##' @describeIn ezAddTable Does the same with a white font and returning the table instead of adding it to the document.
# ezAddTableWhite = function(x, bgcolors=NULL, valign="middle", border=1, head=""){
#   if (is.null(bgcolors)){
#     bgcolors = matrix("#ffffff", nrow=nrow(x), ncol=ncol(x))
#   }
#   ##x = cbind(rownames(x),x)
#   x = as.html(pot(paste('<font color="white">', x, '</font>')))
#   bodyCells = cellProperties(border.width=border, vertical.align=valign)
#   table = FlexTable(x, header.columns = FALSE, body.cell.props=bodyCells,
#                     header.cell.props=cellProperties(border.width = border))
#   table = setFlexTableBackgroundColors(table, j=1:ncol(x), colors=bgcolors)
#   table = addHeaderRow(table, colnames(x))
#   return(table)
# }
