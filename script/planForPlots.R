
require(ReporteRs)

cd = getwd()
setwdNew("./scratch")

theDoc = bsdoc(title = 'My document')
theDoc = addTitle(theDoc, "A title")
irislink = ezImageFileLink(EzPlotterIris$new(name="Iris"),
                    main="Iris sepal shape", xlab="Iris sepal length", ylab="Iris sepal width",
                    width=600,height=600)
theDoc = addParagraph(theDoc, irislink)
types = data.frame(matrix(rep(1:10, each=10), 10))
volcanolink = ezImageFileLink(EzPlotterVolcano$new(name="Volcano2", log2Ratio=1:100, pValue=rep(10^(-4:5), each=10)),
                    xlim=c(0,100), ylim=c(-5,4), pch=16, isPresent=1:50, types=types,
                    colors=rainbow(ncol(types)), legendPos="bottomleft",
                    main="Volcano2")
theDoc = addParagraph(theDoc, volcanolink)
scascalink = ezImageFileLink(EzPlotterXYScatterScatter$new(name="ScatterScatter", xVec=1:10, yVec=1:10),
                    main="ScatterScatter")
theDoc = addParagraph(theDoc, scascalink)
x = ezMatrix("", rows=1, cols=c("iris", "xy", "iris2"))
x[1, "iris"] = ezImageFileLink(EzPlotterIris$new(name="Iris2"))
x[1, "xy"] = ezImageFileLink(EzPlotterXYScatterScatter$new(name="ScatterScatterScatter", xVec=1:10, yVec=1:10))
x[1, "iris2"] = ezImageFileLink(EzPlotterIris$new(name="Iris3"))
theDoc = addFlexTable(theDoc, ezFlexTable(x))
vier = ezImageFileLink(EzPlotterIris$new(name="Iris4"))
fuenf = ezImageFileLink(EzPlotterIris$new(name="Iris5"))
theDoc = addFlexTable(theDoc, ezFlexTable(cbind(vier,fuenf)))
smoscafilink = ezImageFileLink(EzPlotterSmoothScatter$new(name="SmoothScatter", y=data.frame(a=1:10,b=21:30)), file="SmoothScatterplot.png", cex=2)
theDoc = addParagraph(theDoc, smoscafilink)
smoscafilink2 = ezImageFileLink(EzPlotterSmoothScatter$new(name="SmoothScatterMultiple", y=data.frame(a=1:10,b=21:30,c=90:81), x=51:60))
theDoc = addParagraph(theDoc, smoscafilink2)
scatterlink = ezImageFileLink(EzPlotterScatter$new(name="Scatter", y=data.frame(a=1:10,b=21:30)))
theDoc = addParagraph(theDoc, scatterlink)
scatterlink2 = ezImageFileLink(EzPlotterScatter$new(name="ScatterMultiple", y=data.frame(a=1:10,b=21:30,c=90:81), x=51:60))
theDoc = addParagraph(theDoc, scatterlink2)
allscatterlink = ezImageFileLink(EzPlotterAllPairScatter$new(x=matrix(1:10,5)))
theDoc = addParagraph(theDoc, allscatterlink)
allscatterlink2 = ezImageFileLink(EzPlotterAllPairScatter$new(name="allshallscatter", x=matrix(1:25,5)))
theDoc = addParagraph(theDoc, allscatterlink2)

theList = list()
for (i in 1:4){
  theList[[i]] = ezImageFileLink(EzPlotterIris$new(name=paste("Iris",i+10,sep="")))
}
theDoc = addFlexTable(theDoc, ezFlexTable(rbind(theList[1:4])))

# x = ezMatrix("", rows=1, cols=sampleNames)
# for (nm in sampleNames){
#   x[1, nm] = ezImageFileLink(EzPlotterIris$new(name="Iris2", data[sm]))
# }
# theDoc = addFlexTable(theDoc, ezFlexTable(x))

writeDoc(theDoc, "my.html")
setwd(cd)
