
require(ezRun)
require(ReporteRs)

ezImageFileLink3 = function(plotCmd, file=NULL, name="imagePlot", plotType="plot", mouseOverText="my mouse over",
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
    imgFilePot = pot(paste('<img src="', file, '" title="', mouseOverText, '"/>'),
                     hyperlink = pdfName)
  } else {
    imgFilePot = pot(paste('<img src="', pngName, '" title="', mouseOverText, '"/>'))
  }
  return(as.html(imgFilePot))
}

## normal suage
x = 1:10
plotCmd = expression({
  plot(x)
  text(2,1, "my Text")
})

ezImageFileLink3(plotCmd)


## example for a nested call to ezImageFileLink3
fileLinkCaller = function(plotCmd){
  ezImageFileLink3(plotCmd, envir=parent.frame())
}
myFunc = function(){
  x = 1:10
  plotCmd = expression({
    plot(x)
    text(2,1, "my Text")
  })
  fileLinkCaller(plotCmd)
}
myFunc()




cd = getwd()
setwdNew("~/scratch")

theDoc = bsdoc(title = 'My document')
theDoc = addTitle(theDoc, "A title")
x = 1:10
volcanoCmd = expression({ezVolcano(x, 1/x)})
volcanolink = ezImageFileLink3(volcanoCmd, name = "volcano", mouseOverText = "volcano mouse over", plotType = "volcano")
theDoc = addParagraph(theDoc, volcanolink)
writeDoc(theDoc, "myPlot.html")
setwd(cd)

