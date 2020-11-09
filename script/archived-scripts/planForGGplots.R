


EzPlotter2 =
  setRefClass("EzPlotter2",
              fields = c("name", "ggPlot", "helpText", "mouseOverText"),
              methods = list(
                initialize = function(ggPlot, name=NULL)
                {
                  if (ezIsSpecified(name)){
                    name <<- name
                  } else {
                    name <<- "EzPlotter2"
                  }
                  ggPlot <<- ggPlot
                  helpText <<- "Default helpText"
                  mouseOverText <<- "Default mouseOverText"
                },
                plotPng = function(file=NULL, width=480, height=480)
                {
                  "Creates a .png file of a plot."
                  if (ezIsSpecified(file)) {
                    filename = file
                  } else {
                    filename = paste0(name, ".png")
                  }
                  png(filename = filename, width, height)
                  plot(ggPlot)
                  dev.off()
                  return(filename)
                },
                plotPdf = function(file=NULL, width=480, height=480)
                {
                  "Creates a .pdf file of a plot."
                  if (ezIsSpecified(file)) {
                    filename = file
                  } else {
                    filename = paste0(name, ".pdf")
                  }
                  width = round(width/72, digits=2)
                  height = round(height/72, digits=2)
                  pdf(file = filename, width, height)
                  plot(ggPlot)
                  dev.off()
                  return(filename)
                }
              )
  )

ezImageFileLink2 = function(ggPlot, file=NULL, mouseOverText=ezPlotter$mouseOverText,
                            addPdfLink=TRUE, width=480, height=480){
  ezPlotter = EzPlotter2$new(ggPlot)
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

ezGGplot = function(...){
  ggPlot = ggplot(...) + theme_bw(base_size=14) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  return(ggPlot)
}

require(ReporteRs)
cd2 = getwd()
#setwdNew("~/tmp")
#setwdNew("./scratch")
theDoc = bsdoc(title = 'My document')
theDoc = addTitle(theDoc, "A title")
myPlot = qplot(iris$Sepal.Length, iris$Sepal.Width)
myLink = ezImageFileLink2(myPlot, file="test.png")
theDoc = addParagraph(theDoc, myLink)
writeDoc(theDoc, "my.html")
#setwd(cd2)


ggVolcano = function(log2Ratio, pValue, xlim=NULL, ylim=NULL,
                     isPresent=NULL, types=NULL, pch=16,
                     colors=rainbow(ncol(types)),
                     legendPos="bottomright", yType="p-value"){
  
  yValues = -log10(pValue)
  if (is.null(xlim)){
    lrMaxAbs = max(abs(log2Ratio), na.rm=TRUE)
    xm = min(lrMaxAbs, 5)
    xlim = c(-xm, xm)
    log2Ratio = shrinkToRange(log2Ratio, theRange = xlim)
  }
  if (is.null(ylim)){
    ym = min(max(yValues, na.rm=TRUE), 10)
    ylim = c(0, ym)
    yValues = shrinkToRange(yValues, theRange=ylim)
  }
  
  ggPlot = ezGGplot(data=df, aes(log2Ratio, yValues))
  ggPlot = ggPlot + theme(aspect.ratio=1) +
    xlab("log2 ratio") + ylab(paste0("-log10(", yType, ")")) +
    xlim(xlim) + ylim(ylim)
  ggPlot = ggPlot + geom_point(color="grey", shape=pch)
  
  if (!is.null(isPresent)){
    ggPlot = ggPlot + geom_point(color="black", shape=pch, aes(log2Ratio[isPresent], yValues[isPresent]))
  }
  if (!is.null(types)){
    for (j in 1:ncol(types)){
      loopHelp = paste0("geom_point(shape=", pch, ",color='", colors[j], "',aes(", log2Ratio[types[,j]], ",", yValues[types[,j]], "))")
      # loopHelp = paste0("geom_point(shape=", pch, ",aes(", log2Ratio[types[,j]], ",", yValues[types[,j]], ",color='", colors[j], "'))") ## for the legend
      ggPlot = ggPlot + eval(parse(text=loopHelp))
    }
    ggPlot = ggPlot + scale_color_manual(values=colors)
    if (!is.null(legendPos)){
      # legend(legendPos, colnames(types), col=colors, cex=1.2, pch=20, bty="o", pt.bg="white")
      #       ggPlot = ggPlot + theme(legend.position=legendPos) + geom_point(aes())
      #       ggPlot
    }
  }
  return(ggPlot)
}

ggVolcano(log2Ratio = log2Ratio, pValue=yValues)



ggplot(data=ezFrame(log2Ratio=log2Ratio, yValues=yValues)) + geom_point(aes(x=log2Ratio, y=yValues))  + 
  theme_bw(base_size=14) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

ggplot(data=ezFrame(log2Ratio=log2Ratio, yValues=yValues)) + geom_point(aes(x=log2Ratio, y=yValues))

df = ezFrame(log2Ratio=log2Ratio, yValues=yValues)

ggplot() + geom_point(aes(x=log2Ratio, y=yValues))

ggplot(df, aes(x=log2Ratio,y=yValues)) + geom_point()


types = ezFrame(row.names=1:100)
for (x in unique(yValues)){
  types[[as.character(x)]] = x == yValues
}
colors=rainbow(ncol(types))

cp = ggplot() + geom_point(aes(x=log2Ratio, y=yValues))
for (j in 1:ncol(types)){
  use = types[[j]]
  cp = cp +geom_point(data=df[use, ], color=colors[j], aes(x=log2Ratio, y=yValues))
}
cp

cp = ggplot() + geom_point(aes(x=log2Ratio, y=yValues))
for (j in 1:ncol(types)){
  use = types[[j]]
  cp = cp + geom_point(data=df[use, ], color=colors[j], aes(x=log2Ratio, y=yValues))
}
cp

names(colors) = letters[1:length(colors)]
cp + scale_color_manual(name="mylegend", values=colors) + theme(legend.position="right")


cp = ggplot() + geom_point(aes(x=log2Ratio, y=yValues))
for (j in 1:ncol(types)){
  use = types[[j]]
  cp = cp +geom_point(color=colors[j], aes(x=0.5*log2Ratio[use], y=yValues[use]))
}
cp



  theme_bw(base_size=14) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())


# ggtitle("title")

# log2Ratio = 1:100
# yValues = -log10(rep(10^(-4:5), each=10))
# xlim=c(0,100)
# ylim=c(-5,4)
# pch=16
# isPresent=1:50
# types = data.frame(matrix(rep(1:10, each=10), 10))
# colors=rainbow(ncol(types))
# legendPos="bottomright"
# cex.main=1.0
# cex=0.8
# yType="p-value"

bsdoc

  
  