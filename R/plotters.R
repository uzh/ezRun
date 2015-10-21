###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


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
                plotPng = function(file=NULL, width=480, height=480,
                                   nPlotsPerRow=6, ...)
                {
                  "Creates a .png file of plots."
                  if (ezIsSpecified(file)) {
                    filename = file
                  } else {
                    filename = paste(name, ".png", sep="")
                  }
                  
                  if (ncol(data$y) == 2 & is.null(data$x)) {
                    png(filename = filename, width, height)
                    .self$plot(...)
                  } else {
                    nPlots = ncol(data$y)
                    nImgRow <- ceiling(nPlots / nPlotsPerRow)
                    nImgCol <- min(nPlots, nPlotsPerRow)
                    png(filename = filename, height=nImgRow * height, width=nImgCol * width)
                    .self$plot(...)
                  }
                  dev.off()
                  return(filename)
                },
                plotPdf = function(file=NULL, width=480, height=480,
                                   nPlotsPerRow=6, ...)
                {
                  "Creates a .pdf file of smooth scatter plots."
                  if (ezIsSpecified(file)) {
                    filename = file
                  } else {
                    filename = paste(name, ".pdf", sep="")
                  }
                  width = round(width/72, digits=2)
                  height = round(height/72, digits=2)
                  
                  if (ncol(data$y) == 2 & is.null(data$x)) {
                    pdf(file = filename, width, height)
                    .self$plot(...)
                  } else {
                    nPlots = ncol(data$y)
                    nImgRow <- ceiling(nPlots / nPlotsPerRow)
                    nImgCol <- min(nPlots, nPlotsPerRow)
                    pdf(file = filename, height=nImgRow * height, width=nImgCol * width)
                    .self$plot(...)
                  }
                  dev.off()
                  return(filename)
                },
                plot = function(cex=0.8, cex.main=1.0, lim=range(data$x, data$y, na.rm=TRUE),
                                xlab=NULL, ylab=NULL, nPlotsPerRow=6, ...)
                {
                  yValues = as.matrix(data$y)
                  if (is.null(ylab)){
                    ylab=colnames(data$y)
                  }
                  
                  # treat the special case when the reference is not given but there are only two plots
                  if (ncol(data$y) == 2 & is.null(data$x)){
                    par(cex.main=cex.main, cex=cex)
                    smoothScatter(log2(data$y[ ,1]), log2(data$y[ ,2]), xlim=log2(lim), ylim=log2(lim),
                                  xlab=ylab[1], ylab=ylab[2], ...)
                    abline(0, 1, col="blue")
                    return()
                  }
                  
                  ## all other cases
                  nPlots = ncol(data$y)
                  nImgRow <- ceiling(nPlots / nPlotsPerRow)
                  nImgCol <- min(nPlots, nPlotsPerRow)
                  par(mfrow=c(nImgRow, nImgCol))
                  par(cex.main=cex.main, cex=cex)
                  
                  if (nPlots == 1){
                    main = ""
                  } else {
                    main = ylab
                    ylab[] = ""
                  }
                  for (i in 1:nPlots){
                    if (is.null(data$x)){
                      xVal = apply(data$y[ , -i, drop=FALSE], 1, ezGeomean, na.rm=TRUE)
                      if (is.null(xlab)){
                        xlab="Average"
                      }
                    } else {
                      if (is.null(dim(data$x))){
                        xVal = data$x
                      } else {
                        xVal = data$x[ ,i]
                        xlab = colnames(data$x)[i]
                      }
                    }
                    smoothScatter(log2(xVal), log2(data$y[ ,i]), xlim=log2(lim), ylim=log2(lim),
                                  main=main[i], xlab=ylab[1], ylab=ylab[2], ...)
                    abline(0, 1, col="blue")
                  }
                }
              )
  )

EzPlotterScatter =
  setRefClass("EzPlotterScatter",
              contains="EzPlotterSmoothScatter",
              methods=list(
                initialize = function(name=NULL, x=NULL, y)
                {
                  if (ezIsSpecified(name)){
                    name <<- name
                  } else {
                    name <<- "EzPlotterScatter"
                  }
                  data <<- list(x=x, y=y)
                  helpText <<- "Scatter."
                  mouseOverText <<- "Showing mouseOver text."
                },
                plot = function(cex=0.8, cex.main=1.0, lim=range(data$x, data$y, na.rm=TRUE),
                                xlab=NULL, ylab=NULL, nPlotsPerRow=6, isPresent=NULL, types=NULL,
                                shrink=FALSE, pch=16, colors=rainbow(ncol(types)),
                                legendPos="bottomright", ...)
                {
                  yValues = as.matrix(data$y)
                  if (is.null(ylab)){
                    ylab=colnames(data$y)
                  }
                  
                  # treat the special case when the reference is not given but there are only two plots
                  if (ncol(data$y) == 2 & is.null(data$x)){
                    if (is.null(dim(isPresent))){
                      isPres = isPresent
                    } else {
                      isPres = isPresent[ ,1] | isPresent[ , 2]
                    }
                    par(cex.main=cex.main, cex=cex)
                    XYScatterScatter = EzPlotterXYScatterScatter$new(xVec=data$y[ ,1], yVec=data$y[ ,2])
                    XYScatterScatter$plot(xlim=lim, ylim=lim, isPresent=isPres, types=types, pch=pch,
                                          colors=colors, legendPos=legendPos, shrink=shrink,
                                          xlab=ylab[1], ylab=ylab[2], ...)
                    return()
                  }
                  
                  ## all other cases
                  nPlots = ncol(data$y)
                  nImgRow <- ceiling(nPlots / nPlotsPerRow)
                  nImgCol <- min(nPlots, nPlotsPerRow)
                  par(mfrow=c(nImgRow, nImgCol))
                  par(cex.main=cex.main, cex=cex)
                  
                  if (nPlots == 1){
                    main = ""
                  } else {
                    main = ylab
                    ylab[] = ""
                  }
                  for (i in 1:nPlots){
                    if (is.null(dim(isPresent))){
                      isPres = isPresent
                    } else {
                      isPres = isPresent[ ,i]
                    }
                    if (is.null(data$x)){
                      xVal = apply(data$y[ , -i, drop=FALSE], 1, ezGeomean, na.rm=TRUE)
                      if (is.null(xlab)){
                        xlab="Average"
                      }
                    } else {
                      if (is.null(dim(data$x))){
                        xVal = data$x
                      } else {
                        xVal = data$x[ ,i]
                        xlab = colnames(data$x)[i]
                      }
                    }
                    par(mar=c(4.1, 3.1, 4.1, 0.1))
                    XYScatterScatter = EzPlotterXYScatterScatter$new(xVec=xVal, yVec=data$y[ ,i])
                    XYScatterScatter$plot(xlim=lim, ylim=lim, isPresent=isPres, types=types, pch=pch,
                                          colors=colors, legendPos=legendPos, shrink=shrink,
                                          main=main[i], xlab=xlab, ylab=ylab[i], ...)
                  }
                }
              )
  )

EzPlotterFastqScreen =
  setRefClass("EzPlotterFastqScreen",
              contains="EzPlotter",
              methods=list(
                initialize = function(name=NULL, x)
                {
                  if (ezIsSpecified(name)){
                    name <<- name
                  } else {
                    name <<- "EzPlotterFastqScreen"
                  }
                  data <<- x
                  helpText <<- "Contains the methods for the fastqscreen plots."
                  mouseOverText <<- "Showing mouseOver text."
                },
                plot = function(addText=FALSE, ...)
                {
                  par(mar=c(10.1, 4.1, 4.1, 2.1))
                  bplot = barplot(data, ...)
                  if (addText) {
                    text(y=data+5, x=bplot, labels=paste(as.character(data), '%', sep=''), cex=0.7, xpd=TRUE)
                  }
                }
              )
  )





