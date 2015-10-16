
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
                  png(filename = paste(name, ".png", sep=""), width, height)
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
                  graphics::plot(data$Sepal.Length, data$Sepal.Width, ...)
                }
              )
  )

EzPlotterVolcano =
  setRefClass("EzPlotterVolcano",
              contains="EzPlotter",
              methods=list(
                initialize = function(name=NULL, log2Ratio, pValue)
                {
                  if (ezIsSpecified(name)){
                    name <<- name
                  } else {
                    name <<- "EzPlotterVolcano"
                  }
                  data <<- list(log2Ratio=log2Ratio, pValue=pValue)
                  helpText <<- "Volcanos are hot."
                  mouseOverText <<- "Showing mouseOver text."
                },
                plot=function(xlim=NULL, ylim=NULL, pch=16, isPresent=NULL, types=NULL,
                              colors=rainbow(ncol(types)), legendPos="bottomright",
                              cex.main=1.0, cex=0.8, yType="p-value", ...)
                {
                  ## These are here to prevent changing the data field, but perhaps that wouldn't matter?
                  xValues = data$log2Ratio 
                  yValues = -log10(data$pValue)
                  if (is.null(xlim)){
                    lrMaxAbs = max(abs(xValues), na.rm=TRUE)
                    xm = min(lrMaxAbs, 5)
                    xlim = c(-xm, xm)
                    xValues = shrinkToRange(xValues, theRange = xlim)
                  }
                  if (is.null(ylim)){
                    ym = min(max(yValues, na.rm=TRUE), 10)
                    ylim = c(0, ym)
                    yValues = shrinkToRange(yValues, theRange=ylim)
                  }
                  par(pty="s")
                  par(cex.main=cex.main, cex=cex)
                  graphics::plot(xValues, yValues, pch=pch, xlim=xlim, ylim=ylim,
                       col="gray", xlab="log2 ratio", ylab=paste("-log10(", yType, ")" ,sep=""),
                       ...)
                  if (!is.null(isPresent)){
                    points(xValues[isPresent], yValues[isPresent], col="black", pch=pch)
                  }
                  if (!is.null(types)){
                    for (j in 1:ncol(types)){
                      points(xValues[types[,j]], yValues[types[,j]], col=colors[j], pch=pch)
                    }
                    if (!is.null(legendPos)){
                      legend(legendPos, colnames(types), col=colors, cex=1.2, pch=20, bty="o", pt.bg="white")
                    }
                  }
                }
              )
  )

EzPlotterXYScatterScatter =
  setRefClass("EzPlotterXYScatterScatter",
              contains="EzPlotter",
              methods=list(
                initialize = function(name=NULL, xVec, yVec)
                {
                  if (ezIsSpecified(name)){
                    name <<- name
                  } else {
                    name <<- "EzPlotterXYScatterScatter"
                  }
                  data <<- list(xVec=xVec, yVec=yVec)
                  helpText <<- "ScatterScatter."
                  mouseOverText <<- "Showing mouseOver text."
                },
                plot=function(xlim=range(data$xVec, data$yVec, na.rm=TRUE), ylim=xlim,
                              isPresent=NULL, types=NULL, absentColor="gray",
                              frame=TRUE, axes=TRUE, shrink=FALSE,
                              pch=16, colors=rainbow(ncol(types)), legendPos="bottomright", ...)
                {
                  ## These are here to prevent changing the data field, but perhaps that wouldn't matter?
                  xValues = data$xVec
                  yValues = data$yVec
                  if (shrink){
                    xValues = shrinkToRange(xValues, xlim)
                    yValues = shrinkToRange(yValues, ylim)
                  }
                  par(pty="s")
                  if (is.null(isPresent)){
                    graphics::plot(xValues, yValues, log="xy", pch=pch, xlim=xlim, ylim=ylim,
                         col="black", frame=frame, axes=axes, ...)
                  } else {
                    graphics::plot(xValues, yValues, log="xy", pch=pch, xlim=xlim, ylim=ylim,
                         col=absentColor, frame=frame, axes=axes, ...)
                    points(xValues[isPresent], yValues[isPresent], col="black", pch=pch, ...)
                  }
                  if (!is.null(types) && ncol(types) > 0){
                    for (j in 1:ncol(types)){
                      points(xValues[types[,j]], yValues[types[,j]], col=colors[j], pch=pch, ...)
                    }
                    if (!is.null(legendPos)){
                      legend(legendPos, colnames(types), col=colors, cex=1.2, pt.cex=1.5, pch=20, bty="o", pt.bg="white")
                    }
                  }
                  abline(0, 1, col="blue")
                  abline(log10(2), 1, col="blue", lty=2);
                  abline(-log10(2), 1, col="blue", lty=2);
                }
              )
  )

## add the mouse over text: there seems to be no way to do this without writing html code ourselves.
## put it into a 2x1 table: FlexTable doesn't support images.
## put in the second row the pdf link and the help text: Tried to put them on the same line, but it doesn't seem to work.
addEzImage = function(theDoc, ezPlotter, mouseOverText=ezPlotter$mouseOverText, helpText=ezPlotter$helpText,
                      addPdfLink=TRUE, width=480, height=480, ...) {
  ezPlotter$plotPng(width=width, height=height, ...)
  # theDoc = addImage(theDoc, paste(ezPlotter$name, ".png", sep="")) ## directly the image, but no mouseover
  theDoc = addParagraph(theDoc, pot(paste('<img src="', paste(ezPlotter$name, ".png", sep=""),
                                          '" title="', mouseOverText, '"/>')), par.properties=parCenter())
  theDoc = addParagraph(theDoc, helpText, par.properties=parCenter())
  if (addPdfLink) {
    ezPlotter$plotPdf(width=width, height=height, ...)
    theDoc = addParagraph(theDoc, pot("pdf", hyperlink=paste(ezPlotter$name, ".pdf", sep="")),
                          par.properties=parCenter())
  }
}

## in the report generating scripts I want to write
## if mouseOverText and helpText is null the default help text will be used
theDoc = bsdoc(title = 'My document')
theDoc = addTitle(theDoc, "A title")
theDoc = addEzImage(theDoc, EzPlotterIris$new(name="Iris"),
                    main="Iris sepal shape", xlab="Iris sepal length", ylab="Iris sepal width",
                    width=600,height=600)
theDoc = addEzImage(theDoc, EzPlotterVolcano$new(name="Volcano", log2Ratio=1:10, pValue=10:1))
theDoc = addParagraph(theDoc, pot('<center><table><tr><td><img src="Iris.png"></td>'))
theDoc = addParagraph(theDoc, pot('<td><img src="Iris.png"></td></tr></table></center>'))
types = data.frame(matrix(rep(1:10, each=10), 10))
theDoc = addEzImage(theDoc, EzPlotterVolcano$new(name="Volcano2", log2Ratio=1:100, pValue=rep(10^(-4:5), each=10)),
                    xlim=c(0,100), ylim=c(-5,4), pch=16, isPresent=1:50, types=types,
                    colors=rainbow(ncol(types)), legendPos="bottomleft",
                    main="Volcano2")
theDoc = addEzImage(theDoc, EzPlotterXYScatterScatter$new(name="ScatterScatter", xVec=1:10, yVec=1:10),
                    main="ScatterScatter")
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


