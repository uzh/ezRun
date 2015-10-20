
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
                plotPng = function(file=NULL, width=480, height=480, ...)
                {
                  "Creates a .png file of a plot."
                  if (ezIsSpecified(file)) {
                    filename = file
                  } else {
                    filename = paste(name, ".png", sep="")
                  }
                  png(filename = filename, width, height)
                  .self$plot(...)
                  dev.off()
                  return(filename)
                },
                plotPdf = function(file=NULL, width=480, height=480, ...)
                {
                  "Creates a .pdf file of a plot."
                  if (ezIsSpecified(file)) {
                    filename = file
                  } else {
                    filename = paste(name, ".pdf", sep="")
                  }
                  width = round(width/72, digits=2)
                  height = round(height/72, digits=2)
                  pdf(file = filename, width, height)
                  .self$plot(...)
                  dev.off()
                  return(filename)
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

EzPlotterSmoothScatter =
  setRefClass("EzPlotterSmoothScatter",
              contains="EzPlotter",
              methods=list(
                initialize = function(name=NULL, x=NULL, y)
                {
                  if (ezIsSpecified(name)){
                    name <<- name
                  } else {
                    name <<- "EzPlotterSmoothScatter"
                  }
                  data <<- list(x=x, y=y)
                  helpText <<- "SmoothScatter."
                  mouseOverText <<- "Showing mouseOver text."
                },
                plot = function(isPresent=NULL, cex=0.8, lim=range(x, y, na.rm=TRUE),
                                xlab=NULL, ylab=NULL, nPlotsPerRow=6, plotWidth=350,
                                plotHeight=400, cex.main=1.0, ...)
                {
                  yValues = as.matrix(y)
                  if (is.null(ylab)){
                    ylab=colnames(y)
                  }
                  
                  # treat the special case when the reference is not given but there are only two plots
                  if (ncol(y) == 2 & is.null(x)){
                    par(cex.main=cex.main, cex=cex)
                    smoothScatter(log2(y[ ,1]), log2(y[ ,2]), xlim=log2(lim), ylim=log2(lim),
                                  xlab=ylab[1], ylab=ylab[2], ...)
                    abline(0, 1, col="blue")
                    return()
                  }
                  
                  ## all other cases
                  nPlots = ncol(y)
                  nImgRow <- ceiling(nPlots / nPlotsPerRow)
                  nImgCol <- min(nPlots, nPlotsPerRow)
                  
                  ###### TODO: get change for plotPng of EzPlotter
                  if (!is.null(file)){
                    png(filename=file, height=nImgRow * plotHeight, width=nImgCol * plotWidth);
                    on.exit(dev.off())
                  }
                  ###### or: implement additional function, one which calls plotPng with adjusted width and height
                  ###### and then calls this plot
                  
                  
                  par(mfrow=c(nImgRow, nImgCol))
                  par(cex.main=cex.main, cex=cex)
                  if (nPlots == 1){
                    main = ""
                  } else {
                    main = ylab
                    ylab[] = ""
                  }
                  for (i in 1:nPlots){
                    if (is.null(x)){
                      xVal = apply(y[ , -i, drop=FALSE], 1, ezGeomean, na.rm=TRUE)
                      if (is.null(xlab)){
                        xlab="Average"
                      }
                    } else {
                      if (is.null(dim(x))){
                        xVal = x
                      } else {
                        xVal = x[ ,i]
                        xlab = colnames(x)[i]
                      }
                    }
                    smoothScatter(log2(xVal), log2(y[ ,i]), xlim=log2(lim), ylim=log2(lim),
                                  main=main[i], xlab=ylab[1], ylab=ylab[2], ...)
                    abline(0, 1, col="blue")
                  }
                }
              )
  )

## obsolete? can be done with ezImageFileLink and addParagraph as well
addEzImage = function(theDoc, ezPlotter, file=NULL, mouseOverText=ezPlotter$mouseOverText, helpText=ezPlotter$helpText,
                      addPdfLink=TRUE, width=480, height=480, ...) {
  pngName = ezPlotter$plotPng(file=file, width=width, height=height, ...)
  if (addPdfLink) {
    pdfName = ezPlotter$plotPdf(file=sub(".png$", ".pdf", file), width=width, height=height, ...)
    theDoc = addParagraph(theDoc, pot(paste('<img src="', pngName, '" title="', mouseOverText, '"/>'),
                                            hyperlink = pdfName), par.properties=parCenter())
  } else {
    theDoc = addParagraph(theDoc, pot(paste('<img src="', pngName,
                                            '" title="', mouseOverText, '"/>')), par.properties=parCenter())
  }
  theDoc = addParagraph(theDoc, helpText, par.properties=parCenter())
  return(theDoc)
}

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

cd = getwd()
setwdNew("./scratch")

theDoc = bsdoc(title = 'My document')
theDoc = addTitle(theDoc, "A title")
theDoc = addEzImage(theDoc, EzPlotterIris$new(name="Iris"),
                    main="Iris sepal shape", xlab="Iris sepal length", ylab="Iris sepal width",
                    width=600,height=600)
bla = ezImageFileLink(EzPlotterIris$new(name="Iris"),
                      main="Iris sepal shape", xlab="Iris sepal length", ylab="Iris sepal width",
                      width=600,height=600)
theDoc = addParagraph(theDoc, bla)
theDoc = addEzImage(theDoc, EzPlotterVolcano$new(name="Volcano", log2Ratio=1:10, pValue=10:1))
types = data.frame(matrix(rep(1:10, each=10), 10))
theDoc = addEzImage(theDoc, EzPlotterVolcano$new(name="Volcano2", log2Ratio=1:100, pValue=rep(10^(-4:5), each=10)),
                    xlim=c(0,100), ylim=c(-5,4), pch=16, isPresent=1:50, types=types,
                    colors=rainbow(ncol(types)), legendPos="bottomleft",
                    main="Volcano2")
theDoc = addEzImage(theDoc, EzPlotterXYScatterScatter$new(name="ScatterScatter", xVec=1:10, yVec=1:10),
                    main="ScatterScatter")

x = ezMatrix("", rows=1, cols=c("iris", "xy", "iris2"))
x[1, "iris"] = ezImageFileLink(EzPlotterIris$new(name="Iris2"))
x[1, "xy"] = ezImageFileLink(EzPlotterXYScatterScatter$new(name="ScatterScatterScatter", xVec=1:10, yVec=1:10))
x[1, "iris2"] = ezImageFileLink(EzPlotterIris$new(name="Iris3"))
theDoc = addFlexTable(theDoc, ezImageTable(x))

# x = ezMatrix("", rows=1, cols=sampleNames)
# for (nm in sampleNames){
#   x[1, nm] = ezImageFileLink(EzPlotterIris$new(name="Iris2", data[sm]))
# }
# theDoc = addFlexTable(theDoc, ezImageTable(x))

writeDoc(theDoc, "my.html")
setwd(cd)

