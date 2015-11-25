
# together with ReporteRs:
# maybe not possible (at least in no easy way) as each package has their own classes
# but it's no problem if the table is on a separate html page

require(ReporteRs)
require(DT)
doc = openBsdocReport("testing DT")
x = datatable(iris)
link = "exampleTable.html"
saveWidget(x, link)
doc = addParagraph(doc, pot(link, hyperlink=link))
closeBsdocReport(doc, "example.html")

datatable(iris, extensions=c("ColReorder"), options = list(dom = 'Rlfrtip'))
datatable(iris, extensions=c("ColReorder", "ColVis"), options = list(dom = 'C<"clear">lfrtip'))
## but how to add both ColReorder and ColVis? they need a different character for dom in the options list...
datatable(iris, extensions=c("ColReorder", "ColVis", "TableTools"), options = list(dom = 'T<"clear">lfrtip', tableTools = list(sSwfPath = copySWF())))
## I think these three extensions could be useful, but all of them need a different dom in options...

## example with modified start nRows, modified column alignment, sorting and filtering each column.
datatable(iris, extensions="ColVis", filter="top", options = list(dom = 'C<"clear">lfrtip',
                                                                  pageLength = 25,
                                                                  columnDefs = list(list(className = 'dt-left', targets=1:2)),
                                                                  order = list(list(1, 'asc'), list(2, 'desc'))))



