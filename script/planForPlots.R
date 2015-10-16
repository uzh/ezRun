
require(ReporteRs)

# htmlFile = "example_html"
# param = ezParam()
# title = "Der Titel"
# dataset = iris
# writeHtmlReportWithReporteRs(htmlFile, param, title, dataset)
writeHtmlReportWithReporteRs = function(htmlFile, param=param, title="", dataset=NULL){
  html = bsdoc(title = htmlFile)
  html = addTitle(html, pot(paste("<center>", title, "</center>", sep = "")))
  file.copy(ezBannerFile(), ".")
  html = addImage(html, "banner.png", par.properties = parCenter())
  html = addParagraph(html, "test", par.properties = parRight())
  html = addParagraph(html, pot("<center>Test</center>"))
  mymenu = BootstrapMenu( title = 'FGCZ', link = 'http://www.fgcz.ch/')
  mydd = DropDownMenu( label = 'Mon menu' )
  mydd = addLinkItem( mydd, label = 'GitHub', 'http://github.com/')
  mydd = addLinkItem( mydd, separator.after = TRUE)
  mydd = addLinkItem( mydd, label = 'Wikipedia', 'http://www.wikipedia.de')
  mymenu = addLinkItem( mymenu, label = 'ReporteRs', 'http://github.com/davidgohel/ReporteRs')
  mymenu = addLinkItem( mymenu, dd = mydd )
  html = addBootstrapMenu( html, mymenu )
  writeDoc(html, "test.html")
}

## all plots would be generated from a plotter class
EzPlotter =
  setRefClass("EzPlotter",
              fields = c("name", "data", "helpText", "mouseOverText"),
              methods = list(
                plot = function(...)
                {
                  "Plots \\code{data} with the default plot function from the graphics package."
                  graphics::plot(data, ...)
                },
                plotPng = function(width=NULL, height=NULL, ...)
                {
                  "Creates a .png file of a plot."
                  png(file = paste(name, ".png", sep=""), width, height)
                  .self$plot(...)
                  dev.off()
                },
                plotPdf = function(width=NULL, height=NULL, ...)
                {
                  "Creates a .pdf file of a plot."
                  width = width/72
                  height = height/72
                  pdf(file = paste(name, ".pdf", sep=""), width, height)
                  .self$plot(...)
                  dev.off()
                },
                writeData = function()
                {
                  
                }
              )
  )

EzPlotterIris =
  setRefClass("EzPlotterIris",
              contains="EzPlotter",
              methods=list(
                initialize = function(name=NULL)
                {
                  if (ezIsSpecified(name)){
                    name <<- name
                  } else {
                    name <<- "EzPlotterIris"
                  }
                  data <<- iris
                  helpText <<- "Iris is a flower dataset."
                  mouseOverText <<- "Showing mouseOver text."
                },
                plot=function(...)
                {
                  graphics::plot(.self$data$Sepal.Length, .self$data$Sepal.Width, ...)
                }
              )
  )

# testplot = EzPlotterIris$new()
# testplot$plot()
# testplot$plotPng()
# testplot$plotPdf()

addEzImage = function(theDoc, ezPlotter, mouseOverText=ezPlotter$mouseOverText, helpText=ezPlotter$helpText,
                      addPdfLink=TRUE, width=480, height=480, ...) {
  ## add the mouse over text: there seems to be no way to do this without writing html code ourselves.
  ## put it into a 2x1 table: FlexTable doesn't support images.
  ## put in the second row the pdf link and the help text: Tried to put them on the same line, but it doesn't seem to work.
  ezPlotter$plotPng(width, height, ...)
  # theDoc = addImage(theDoc, paste(ezPlotter$name, ".png", sep="")) ## directly the image, but no mouseover
  theDoc = addParagraph(theDoc, pot(paste('<img src="', paste(ezPlotter$name, ".png", sep=""),
                                          '" title="', mouseOverText, '"/>')), par.properties=parCenter())
  theDoc = addParagraph(theDoc, helpText, par.properties=parCenter())
  if (addPdfLink) {
    ezPlotter$plotPdf(width, height, ...)
    theDoc = addParagraph(theDoc, pot("pdf", hyperlink=paste(ezPlotter$name, ".pdf", sep="")),
                          par.properties=parCenter())
  }
}

## in the report generating scripts I want to write
## if mouseOverText and helpText is null the default help text will be used
theDoc = bsdoc(title = 'My document')
theDoc = addTitle(theDoc, "A title")
theDoc = addEzImage(theDoc, EzPlotterIris$new(name="myIrisPlot"),
                    main="Iris sepal shape", xlab="Iris sepal length", ylab="Iris sepal width",
                    width=800,height=800)
writeDoc(theDoc, "my.html")

## alternative implementation (does the same, but looks uglier)
# bla = pot_img("EzPlotterIris.png")
# sopar = set_of_paragraphs(bla, pot(helpText), pot("pdf", hyperlink="EzPlotterIris.pdf"))
# theDoc = addParagraph(theDoc, sopar, par.properties=parCenter())

## doesn't work
# ft = FlexTable(numrow=2,numcol=1)
# ft[[1]] = addImage(theDoc, "myirisPlot.png")
# ft[[2]] = addParagraph(theDoc, "testing")
# theDoc = addFlexTable(theDoc, ft)



