
require(ReporteRs)

# cd = getwd()
# setwdNew("./scratch")
# 
# theDoc = bsdoc(title = 'My document')
# theDoc = addTitle(theDoc, "A title")
# irislink = ezImageFileLink(EzPlotterIris$new(name="Iris"),
#                     main="Iris sepal shape", xlab="Iris sepal length", ylab="Iris sepal width",
#                     width=600,height=600)
# theDoc = addParagraph(theDoc, irislink)
# types = data.frame(matrix(rep(1:10, each=10), 10))
# volcanolink = ezImageFileLink(EzPlotterVolcano$new(name="Volcano2", log2Ratio=1:100, pValue=rep(10^(-4:5), each=10)),
#                     xlim=c(0,100), ylim=c(-5,4), pch=16, isPresent=1:50, types=types,
#                     colors=rainbow(ncol(types)), legendPos="bottomleft",
#                     main="Volcano2")
# theDoc = addParagraph(theDoc, volcanolink)
# scascalink = ezImageFileLink(EzPlotterXYScatterScatter$new(name="ScatterScatter", xVec=1:10, yVec=1:10),
#                     main="ScatterScatter")
# theDoc = addParagraph(theDoc, scascalink)
# x = ezMatrix("", rows=1, cols=c("iris", "xy", "iris2"))
# x[1, "iris"] = ezImageFileLink(EzPlotterIris$new(name="Iris2"))
# x[1, "xy"] = ezImageFileLink(EzPlotterXYScatterScatter$new(name="ScatterScatterScatter", xVec=1:10, yVec=1:10))
# x[1, "iris2"] = ezImageFileLink(EzPlotterIris$new(name="Iris3"))
# theDoc = addFlexTable(theDoc, ezFlexTable(x))
# vier = ezImageFileLink(EzPlotterIris$new(name="Iris4"))
# fuenf = ezImageFileLink(EzPlotterIris$new(name="Iris5"))
# theDoc = addFlexTable(theDoc, ezFlexTable(cbind(vier,fuenf)))
# smoscafilink = ezImageFileLink(EzPlotterSmoothScatter$new(name="SmoothScatter", y=data.frame(a=1:10,b=21:30)), file="SmoothScatterplot.png", cex=2)
# theDoc = addParagraph(theDoc, smoscafilink)
# smoscafilink2 = ezImageFileLink(EzPlotterSmoothScatter$new(name="SmoothScatterMultiple", y=data.frame(a=1:10,b=21:30,c=90:81), x=51:60))
# theDoc = addParagraph(theDoc, smoscafilink2)
# scatterlink = ezImageFileLink(EzPlotterScatter$new(name="Scatter", y=data.frame(a=1:10,b=21:30)))
# theDoc = addParagraph(theDoc, scatterlink)
# scatterlink2 = ezImageFileLink(EzPlotterScatter$new(name="ScatterMultiple", y=data.frame(a=1:10,b=21:30,c=90:81), x=51:60))
# theDoc = addParagraph(theDoc, scatterlink2)
# allscatterlink = ezImageFileLink(EzPlotterAllPairScatter$new(x=matrix(1:10,5)))
# theDoc = addParagraph(theDoc, allscatterlink)
# allscatterlink2 = ezImageFileLink(EzPlotterAllPairScatter$new(name="allshallscatter", x=matrix(1:25,5)))
# theDoc = addParagraph(theDoc, allscatterlink2)
# 
# theList = list()
# for (i in 1:4){
#   theList[[i]] = ezImageFileLink(EzPlotterIris$new(name=paste("Iris",i+10,sep="")))
# }
# theDoc = addFlexTable(theDoc, ezFlexTable(rbind(theList[1:4])))
# 
# testLink = "testLink.txt"
# theDoc = addParagraph(theDoc, pot(paste("<a href='", testLink, "' type='application/text'>", testLink, "</a>")))
# 
# # x = ezMatrix("", rows=1, cols=sampleNames)
# # for (nm in sampleNames){
# #   x[1, nm] = ezImageFileLink(EzPlotterIris$new(name="Iris2", data[sm]))
# # }
# # theDoc = addFlexTable(theDoc, ezFlexTable(x))
# 
# writeDoc(theDoc, "my.html")
# setwd(cd)
# 
# 





EzPlotter2 =
  setRefClass("EzPlotter2",
              fields = c("name", "ggplot", "helpText", "mouseOverText"),
              methods = list(
                initialize = function(ggplot, name=NULL)
                {
                  if (ezIsSpecified(name)){
                    name <<- name
                  } else {
                    name <<- "EzPlotter2"
                  }
                  ggplot <<- ggplot
                  helpText <<- "Default helpText"
                  mouseOverText <<- "Default mouseOverText"
                },
                plotPng = function(file=NULL, width=480, height=480)
                {
                  "Creates a .png file of a plot."
                  if (ezIsSpecified(file)) {
                    filename = file
                  } else {
                    filename = paste(name, ".png", sep="")
                  }
                  png(filename = filename, width, height)
                  plot(ggplot)
                  dev.off()
                  return(filename)
                },
                plotPdf = function(file=NULL, width=480, height=480)
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
                  plot(ggplot)
                  dev.off()
                  return(filename)
                }
              )
  )

ezImageFileLink2 = function(ggplot, file=NULL, mouseOverText=ezPlotter$mouseOverText,
                              addPdfLink=TRUE, width=480, height=480){
  
  ezPlotter = EzPlotter2$new(ggplot)
  
  pngName = ezPlotter$plotPng(file=file, width=width, height=height)
  if (addPdfLink) {
    pdfName = ezPlotter$plotPdf(file=sub(".png$", ".pdf", file), width=width, height=height)
    imgFilePot = pot(paste('<img src="', pngName, '" title="', mouseOverText, '"/>'),
                     hyperlink = pdfName)
  } else {
    imgFilePot = pot(paste('<img src="', pngName, '" title="', mouseOverText, '"/>'))
  }
  return(as.html(imgFilePot))
}

cd2 = getwd()
setwdNew("./scratch")
theDoc = bsdoc(title = 'My document')
theDoc = addTitle(theDoc, "A title")
myPlot = qplot(iris$Sepal.Length, iris$Sepal.Width)
myLink = ezImageFileLink2(myPlot, file="test.png")
theDoc = addParagraph(theDoc, myLink)
writeDoc(theDoc, "my.html")
setwd(cd2)


# EzPlotterVolcano$new(name="Volcano2", log2Ratio=1:100, pValue=rep(10^(-4:5), each=10)),
#                     xlim=c(0,100), ylim=c(-5,4), pch=16, isPresent=1:50, types=types,
#                     colors=rainbow(ncol(types)), legendPos="bottomleft",
#                     main="Volcano2")



# volcanoPlot = qplot(x=1:100, y=)




xlim=NULL
ylim=NULL
pch=16
isPresent=NULL
types=NULL
colors=rainbow(ncol(types))
legendPos="bottomright"
cex.main=1.0
cex=0.8
yType="p-value"

xValues = 1:100
yValues = -log10(rep(10^(-4:5), each=10))
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
               col="gray", xlab="log2 ratio", ylab=paste("-log10(", yType, ")" ,sep=""))
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

