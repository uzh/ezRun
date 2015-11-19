###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


##' @title Gets the sample colors
##' @description Gets the sample colors from the experimental conditions.
##' @param conds conditions obtained from \code{ezConditionsFromDesign()} or \code{ezConditionsFromDataset()}.
##' @param colorNames a character vector containing the names for the colors.
##' @param hueStep a numeric specifying the hue step for \code{hsv()}.
##' @template roxygen-template
##' @seealso \code{\link[grDevices]{hsv}}
##' @return Returns a character vector containing colors in hex format.
getSampleColors = function(conds, colorNames=names(conds), hueStep = 1/50){

	condSet = unique(conds)
	hueDelta = 1 / length(condSet)
	hue = (match(conds, condSet) - 1) * hueDelta
	vs = rep(1, length(conds))
	for (cond in condSet){
		idx = which(cond == conds)
		if (length(idx) <= 8){
		  vs[idx] = seq.int(from=1, by=-0.1, length.out=length(idx))
		} else {
		  vs[idx] = seq.int(from=1, to=0.2, length.out=length(idx))
		}
		if ((length(idx)+2) * hueStep > hueDelta){
			hueStep = hueDelta / (length(idx)+2)
		}
  	hue[idx] = hue[idx] + 0:(length(idx)-1) * hueStep
	}
	result = hsv(hue, vs, vs)
	names(result) = colorNames
	return(result)
}

##' @describeIn getSampleColors Gets the sample pch from the experimental conditions.
getSamplePch = function(conds, pchNames=names(conds)){
  
  cnList = split(pchNames, conds)
  idxList = lapply(cnList, function(cn){0:(length(cn)-1)})
  pch = unsplit(idxList, conds)
  names(pch) = pchNames
  return(pch)
}

##' @describeIn getSampleColors Gets the sample line types from the experimental conditions.
getSampleLty = function(conds, ltyNames=names(conds), maxLineTypes=5){
  
  cnList = split(ltyNames, conds)
  idxList = lapply(cnList, function(cn){rep(1:maxLineTypes, length.out=length(cn))})
  lty = unsplit(idxList, conds)
  names(lty) = ltyNames
  return(lty)
}

##' @title Gets a color scale from blue to red
##' @description Gets a color scale from blue to red in hex format.
##' @param n an integer specifying the amount of colors to split the scale in.
##' @param whiteLevel a numeric specifying the maximum whiteness for red, green and blue.
##' @template roxygen-template
##' @return Returns a character vector containing colors in hex format.
##' @examples
##' getBlueRedScale()
getBlueRedScale <- function(n=256, whiteLevel=0.9){

  # from blue to red going by white
  cs <- character(n)
  n1 <- ceiling(n/2)
  cs[1:n1] <- rgb(seq(from=0, to=whiteLevel, length.out=n1),
                  seq(from=0, to=whiteLevel, length.out=n1),
                  seq(from=1, to=whiteLevel, length.out=n1)
                 )
  n2 <- n - n1
  cs[(n1+1):n] <- rev(rgb(seq(from=1, to=whiteLevel, length.out=n2),
                          seq(from=0, to=whiteLevel*(n2-1)/(n1-1), length.out=n2),
                         seq(from=0, to=whiteLevel*(n2-1)/(n1-1), length.out=n2)
                      ))
  cs
}

##' @title Plots a color scale
##' @description Plots a color scale with colors derived from \code{getBlueRedScale()}.
##' @param colorRange two numerics specifying the range to plot to the axis.
##' @template colors-template
##' @param vertical a logical indicating whether to plot vertically.
##' @param at a numeric vector specifying where to put axis ticks.
##' @param labels a character vector specifying the axis labels.
##' @param by.label a numeric specifying the interval between axis labels.
##' @template roxygen-template
##' @examples
##' ezColorLegend()
ezColorLegend = function(colorRange=c(-3,3), colors=getBlueRedScale(), vertical=TRUE,
															at=seq(from=colorRange[1], to=colorRange[2], by=by.label),
															labels = as.character(at), by.label=0.5){
	pos = (at - colorRange[1])/(colorRange[2]- colorRange[1])
	#n = (colorRange[2] - colorRange[1])*2 + 1
	if (vertical){
	  par(mar=c(2,2,2,4))
	  image(t(as.matrix((1:length(colors)))), axes=FALSE, frame.plot=TRUE, col=colors)
	  axis(4, at=pos, las=2, labels=labels)
	} else {
	  par(mar=c(4,2,2,2))
	  image(as.matrix((1:length(colors))), axes=FALSE, frame.plot=TRUE, col=colors)
	  axis(1, at=pos, las=1, labels=labels)
	}
}


## still used?
.makeTailEffectPlots = function(param, signal, seqAnno, colors=NULL, ylim=c(-2,5)){

  isControl = seqAnno$IsControl
  logSignal = log2(signal[!isControl, ])
  samples = colnames(logSignal)
  sequence = seqAnno[!isControl, "Sequence"]
  gene = seqAnno[!isControl, "Accession [Agilent]"]
  xlim = range(logSignal, na.rm=TRUE)
  tails = c("A", "C", "G", "T")
  tailColors = getSampleColors(tails, colorNames=tails)

  pngNames = character()
  for(sampleName in samples){
    valueByTail = getByTail(logSignal[ , sampleName], NULL, sequence, gene)
    values = unlist(lapply(valueByTail, function(x){x[-length(x)]}))
    effectByTail = lapply(valueByTail, diff)
    effects = unlist(effectByTail)
    effects = shrinkToRange(effects, ylim)
    pngNames[sampleName] = paste0(sampleName, "-tailEffectPlot.png")
    png(pngNames[sampleName], height=400, width=400)
    for (tail in tails){
      #values = unlist(lapply(valueByTailSet[[sampleName]], function(x){x[-1]}))
      use = grep(paste0(tail, "$"), names(effects))
      if (tail == "A"){
        plot(values[use], effects[use], pch=20, main=sampleName, col=tailColors[tail], cex=0.5,
             xlim=xlim, ylim=ylim, xlab="Signal", ylab="Tail Effect")
        abline(h=0, lwd=2, col="gray")
      } else {
        points(values[use], effects[use], pch=20, col=tailColors[tail], cex=0.5)
      }
    }
    legend("bottomright", tails, col=tailColors, bty="n", cex=0.8, pt.bg="white", lty=1, lwd=3 )
    dev.off()
  }
  return(pngNames)
}


## still used?
.makeTailEffectPlotsByTail = function(param, signal, seqAnno, colors=NULL, ylim=c(-2,5)){

  valueByTailSet = list()
  effectByTailSet = list()
  isControl = seqAnno$IsControl
  logSignal = log2(signal[!isControl, ])
  samples = colnames(logSignal)
  sequence = seqAnno[!isControl, "Sequence"]
  gene = seqAnno[!isControl, "Agilent SystematicName"]
  for(sampleName in samples){
    valueByTail = getByTail(logSignal[ , sampleName], NULL, sequence, gene)
    effectByTailSet[[sampleName]] = lapply(valueByTail, diff)
    valueByTailSet[[sampleName]] = valueByTail  
  }
  xlim = range(logSignal)

  tails = c("A", "C", "G", "T")
  tail = tails[1]
  pngNames = character()
  for (tail in tails){
    pngNames[tail] = paste0("tailEffectPlot-", tail, ".png")
    png(pngNames[tail], height=500, width=500)
    for (sampleName in samples){
      values = unlist(lapply(valueByTailSet[[sampleName]], function(x){x[-length(x)]}))
      #values = unlist(lapply(valueByTailSet[[sampleName]], function(x){x[-1]}))
      effects = unlist(effectByTailSet[[sampleName]])
      effects = shrinkToRange(effects, ylim)
      use = grep(paste0(tail, "$"), names(effects))
      if (sampleName == colnames(logSignal)[1]){
        plot(values[use], effects[use], pch=20, main=paste("Tail", tail), col=colors[sampleName], 
             xlim=xlim, ylim=ylim, xlab="Signal", ylab="Tail Effect")
        abline(h=0, lwd=2, col="gray")
      } else {
        points(values[use], effects[use], pch=20, col=colors[sampleName])
      }
    }
    legend("topright", samples, col=colors, bty="n", cex=0.8, pt.bg="white", lty=1, lwd=3 )
    dev.off()
  }
  return(pngNames)
}


## still used?
.getByTail = function(values, isPresent, seq, gene, method=mean){
  names(values) = seq
  valueListByGeneAll = split(values, gene)
  seqCounts = sapply(valueListByGeneAll, function(x){length(unique(names(x)))})
  use = seqCounts > 1
  if (!is.null(isPresent)){
    names(isPresent) = seq
    presentListByGeneAll  = split(isPresent, gene)
    isPresentGene = mapply(computePresence, presentListByGeneAll)
    use[!isPresentGene] = FALSE
  }
  valueListByGene = valueListByGeneAll[use]

  sortAndName = function(x){
    x = x[order(names(x))]
    additionalNuc = sapply(strsplit(names(x), ""), tail, n=1)
    additionalNuc[1] = ""
    names(x) = additionalNuc
    return(x)
  }
  valueByTail = lapply(valueListByGene, sortAndName)
  return(valueByTail)
}

##' @title Does a volcano plot
##' @description Does a volcano plot.
##' @param log2Ratio a numeric vector containing the log2 ratios of a result.
##' @param pValue a numeric vector containing the p-Values of a result.
##' @param yType a character specifying the type of the y-value. Gets pasted onto the y-axis of the plot.
##' @template plot-template
##' @template roxygen-template
##' @examples
##' types = data.frame(matrix(rep(1:10, each=10), 10))
##' ezVolcano(log2Ratio=1:100, pValue=rep(10^(-4:5), each=10),
##' pch=16, isPresent=1:50, types=types, colors=rainbow(ncol(types)), legendPos="bottomleft")
ezVolcano = function(log2Ratio, pValue, yType="p-value",
                     lim=list(NA, NA), isPresent=NULL,
                     types=NULL, pch=16, colors=rainbow(ncol(types)), legendPos="bottomright",
                     cex.main=1.0, cex=0.8, ...){
  yValues = -log10(pValue)
  xlim = lim[[1]]
  ylim = lim[[2]]
  
  if (is.na(lim[[1]])){
    lrMaxAbs = max(abs(log2Ratio), na.rm=TRUE)
    xm = min(lrMaxAbs, 5)
    xlim = c(-xm, xm)
    log2Ratio = shrinkToRange(log2Ratio, theRange = xlim)
  }
  if (is.na(lim[[2]])){
    ym = min(max(yValues, na.rm=TRUE), 10)
    ylim = c(0, ym)
    yValues = shrinkToRange(yValues, theRange=ylim)
  }

  par(pty="s")
  par(cex.main=cex.main, cex=cex)
  plot(log2Ratio, yValues, pch=pch, xlim=xlim, ylim=ylim,
      col="gray", xlab="log2 ratio", ylab=paste0("-log10(", yType, ")"),
       ...)
  if (!is.null(isPresent)){
    points(log2Ratio[isPresent], yValues[isPresent], col="black", pch=pch)
  }
  if (!is.null(types)){
    for (j in 1:ncol(types)){
      points(log2Ratio[types[,j]], yValues[types[,j]], col=colors[j], pch=pch)
    }
    if (!is.null(legendPos)){
      legend(legendPos, colnames(types), col=colors, cex=1.2, pch=20, bty="o", pt.bg="white")
    }
  }
}

##' @title Does smooth scatter plots
##' @description Does smooth scatter plots.
##' @param x an optional reference vector or matrix.
##' @param y a matrix of values to plot.
##' @param xlab a character for labeling the reference. If \code{x} is a matrix, its colnames will be used as labels.
##' @param ylab a character vector for the labels of the plots. If \code{ylab} is NULL, the colnames from \code{y} will be used as labels.
##' @param nPlotsPerRow an integer specifying the number of plots per row.
##' @template plot-template
##' @template roxygen-template
##' @examples
##' ezSmoothScatter(y=data.frame(a=1:10, b=21:30, c=41:50))
ezSmoothScatter <- function(x=NULL, y, xlab=NULL, ylab=NULL, nPlotsPerRow=6,
                            lim=range(x, y, na.rm=TRUE), isPresent=NULL,
                            types=NULL, pch=16, colors=rainbow(ncol(types)), legendPos="bottomright",
                            cex.main=1.0, cex=0.8, ...){
  
  y = as.matrix(y)
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
                  main=main[i], xlab=xlab, ylab=ylab[i], ...)
    abline(0, 1, col="blue")
  }
}

##' @title Does scatter plots
##' @description Does scatter plots.
##' @inheritParams ezSmoothScatter
##' @param shrink a logical specifying whether to shrink the values to range.
##' @template roxygen-template
##' @examples
##' ezScatter(y=data.frame(a=1:10, b=21:30, c=41:50))
ezScatter <- function(x=NULL, y, xlab=NULL, ylab=NULL, nPlotsPerRow=6, shrink=FALSE,
                      lim=range(x, y, na.rm=TRUE), isPresent=NULL,
                      types=NULL, pch=16, colors=rainbow(ncol(types)), legendPos="bottomright", 
                      cex.main=1.0, cex=0.8, ...){
  
  y = as.matrix(y)
  if (is.null(ylab)){
    ylab=colnames(y)
  }
  
  # treat the special case when the reference is not given but there are only two plots
  if (ncol(y) == 2 & is.null(x)){
    if (is.null(dim(isPresent))){
      isPres = isPresent
    } else {
      isPres = isPresent[ ,1] | isPresent[ , 2]
    }
    par(cex.main=cex.main, cex=cex)
    ezXYScatterScatter(y[ ,1], y[ ,2], xlim=lim, ylim=lim, isPresent=isPres,
                       types=types, pch=pch, colors=colors, legendPos=legendPos, shrink=shrink,
                       xlab=ylab[1], ylab=ylab[2], ...)
    return()
  }
  
  ## all other cases
  nPlots = ncol(y)
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
    par(mar=c(4.1, 3.1, 4.1, 0.1))
    ezXYScatterScatter(xVal, y[ ,i], xlim=lim, ylim=lim, isPresent=isPres,
                       types=types, pch=pch, colors=colors, legendPos=legendPos, shrink=shrink,
                       main=main[i], xlab=xlab, ylab=ylab[i], ...)
  }
}

##' @describeIn ezScatter Does the XY scatter plot.
ezXYScatterScatter = function(xVec, yVec, absentColor="gray", shrink=FALSE, frame=TRUE, axes=TRUE,
                              xlim=range(xVec, yVec, na.rm=TRUE), ylim=xlim, isPresent=NULL,
                              types=NULL, pch=16, colors=rainbow(ncol(types)), legendPos="bottomright", ...){
  par(pty="s")
  if (shrink){
    xVec = shrinkToRange(xVec, xlim)
    yVec = shrinkToRange(yVec, ylim)
  }
  if (is.null(isPresent)){
    plot(xVec, yVec, log="xy", pch=pch, xlim=xlim, ylim=ylim,
         col="black", frame=frame, axes=axes,
         ...)
  } else {
    plot(xVec, yVec, log="xy", pch=pch, xlim=xlim, ylim=ylim,
         col=absentColor, frame=frame, axes=axes,
         ...)
    points(xVec[isPresent], yVec[isPresent], col="black", pch=pch, ...)
  }
  if (!is.null(types) && ncol(types) > 0){
    for (j in 1:ncol(types)){
      points(xVec[types[,j]], yVec[types[,j]], col=colors[j], pch=pch, ...)
    }
    if (!is.null(legendPos)){
      legend(legendPos, colnames(types), col=colors, cex=1.2, pt.cex=1.5, pch=20, bty="o", pt.bg="white")
    }
  }
  abline(0, 1, col="blue")
  abline(log10(2), 1, col="blue", lty=2);
  abline(-log10(2), 1, col="blue", lty=2);
}

##' @title Does scatter plots of all pairs
##' @description Does scatter plots of all pairs.
##' @param x a matrix containing the data to plot.
##' @param main a character specifying the main title of the plots.
##' @param shrink a logical specifying whether to shrink the values to range.
##' @param xylab a character vector containing the axis labels. If it is NULL, colnames of \code{x} will be used.
##' @template plot-template
##' @template roxygen-template
##' @examples
##' ezAllPairScatter(x=matrix(1:10,5))
ezAllPairScatter = function(x, main="", shrink=FALSE, xylab=NULL,
                            lim=range(x, na.rm=TRUE), isPresent=NULL,
                            types=NULL, pch=16, colors=rainbow(ncol(types)), legendPos="bottomright",
                            cex.main=1.0, cex=0.8, ...){
  nItems = ncol(x)
  if (is.null(xylab)){
    xylab = colnames(x)
  }
  if (nItems > 2){
   par(oma=c(2,2,2,2), mar=c(0,0,0,0), mfcol=c(nItems, nItems))
   par(pch=pch, cex=cex, pty="s", cex.main=cex.main)
   for( i in 1:nItems){
     for(j in 1:nItems){
       if (is.null(dim(isPresent))){
         isPres = isPresent ## this covers also the case where isPresent == NULL
       } else {
         isPres = isPresent[ ,i] | isPresent[ ,j]
       }
       ezXYScatterScatter(x[ ,i], x[, j], xlim=lim, ylim=lim, shrink=shrink, axes=FALSE, frame=TRUE,
                          isPresent=isPres, types=types, pch=pch, colors=colors, legendPos=NULL, ...)
       if (i == 1){
         mtext(xylab[j], 2)
       }
       if (j == nItems){
         mtext(xylab[i], 1)
       }
     }
   }
  } else {
    par(pch=pch, cex=cex, pty="s", cex.main=cex.main)
    if (is.null(dim(isPresent))){
      isPres = isPresent   ## this covers also the case where isPresent == NULL
    } else {
      isPres = isPresent[ ,1] | isPresent[ ,2]
    }
    ezXYScatterScatter(x[ ,1], x[, 2], xlim=lim, ylim=lim, shrink=shrink, xlab=xylab[1], ylab=xylab[2],
                       isPresent=isPresent, types=types, pch=pch, colors=colors, legendPos=legendPos, ...)
  }
  mtext(main, outer=TRUE, line=1)
}

##' @title Does a correlation plot
##' @description Does a correlation plot.
##' @param z the data to plot.
##' @param cond a character vector specifying the conditions.
##' @param condOrder a sorted character vector specifying the order of the conditions.
##' @param main a character specifying the main title of the plots.
##' @param labels a character vector specifying axis labels.
##' @param condLabels a character vector specifying condition labels.
##' @param plotLabels a logical specifying whether to do plot labels.
##' @template colors-template
##' @template roxygen-template
##' @examples
##' ezCorrelationPlot(z=matrix(1:100,10))
ezCorrelationPlot <- function(z, cond=NULL, condOrder=NULL, main="Correlation", 
                              labels=NULL, condLabels=NULL, plotLabels=nrow(z) < 100,
                              colors=NULL){
  
  colorScale <- gray((1:256)/256)
  condNumbers = NULL
  nRow <- nrow(z)
  
  layout(matrix(c(1,3,4,0,2,0), ncol=3, nrow=2, byrow=TRUE),
         widths=c(0.2,0.6,0.2), heights=c(0.8, 0.2))
  #par(pin=c(4,4))
  tmp <- z
  tmp[ tmp == 1] <- NA;
  if ( is.null(labels)){
    labels <- rownames(tmp)
  }
  
  rge <- range(tmp, finite=TRUE);
  #rge <- c(min(tmp, na.rm=TRUE), 1)
  breaks <- 0:256 / 256 * (rge[2] - rge[1]) + rge[1]
  
  # left and bottom plot
  if (is.null(cond)){
    par(mar=c(0, 2, 2, 0));
    image(as.matrix(1), axes=FALSE, frame.plot=FALSE, col="white")
    par(mar=c(2, 0, 0, 2));
    image(as.matrix(1), axes=FALSE, frame.plot=FALSE, col="white")
    condLabels <- NULL
  } else {
    if (is.null(condOrder)){
      condOrder <- unique(cond)
    }
    if (length(cond) != nRow){
      stop("incorrect size of cond option")
    }
    idx <- integer(nRow)
    for ( i in 1:length(condOrder)){
      idx[ condOrder[i] == cond] <- i
    }
    if (is.null(colors)){
      condCol = rainbow(length(condOrder))
      colors = condCol[idx]
    }
    idxOrdered <- order(idx)
    condNumbers <- as.matrix(idx[idxOrdered])
    
    par(mar=c(0, 2, 2, 0));
    image(t((1:nRow)[idxOrdered]), axes=FALSE, frame.plot=FALSE, col=colors)
    par(mar=c(2, 0, 0, 2));
    image(as.matrix((1:nRow)[idxOrdered]), axes=FALSE, frame.plot=FALSE, col=colors)
    
    tmp <- tmp[idxOrdered, idxOrdered]
    if (!is.null(labels)){
      labels <- labels[idxOrdered]
    }
    if (is.null(condLabels)){
      condLabels <- labels
    } else {
      condLabels = condLabels[idxOrdered]
    }
  }
  inc <- (rge[2] - rge[1])/10
  # center plot
  par(mar=c(0, 0, 2, 2));
  image(tmp, main=main, axes=FALSE, col=colorScale, breaks=breaks);
  
  nRow <- dim(tmp)[1]
  if (plotLabels){
    if (!is.null(labels)){
      axis(2, at=(0:(nRow-1)/(nRow-1)), las=2, labels=labels, cex.axis=1.3)
      axis(1, at=(0:(nRow-1)/(nRow-1)), las=2, labels=condLabels, cex.axis=1.3)
    }
  }
  if (!is.null(condNumbers)){
    for (i in 2:length(condNumbers)){
      if (condNumbers[i] != condNumbers[i-1]){
        abline(h=(i-1.5)/(nRow-1))
        abline(v=(i-1.5)/(nRow-1))
      }
    }
  }
  
  # right plot
  par(mar=c(0, 2, 2, 4));
  image(t(as.matrix((1:256))), axes=FALSE, frame.plot=TRUE, col=colorScale)
  axis(4, at=(0:9)/9, las=2, labels=signif(rge[1] + inc * 0:9, 3), cex.axis=1.5)
}

##' @title Shrinks values and plots a histogram
##' @description Shrinks values and plots a histogram.
##' @param x a vector of values to shrink and plot a histogram from.
##' @param range two values specifying the range to shrink \code{x} to.
##' @param step a value specifying the step distance between the breaks for \code{hist()}.
##' @template addargs-template
##' @templateVar fun hist()
##' @template roxygen-template
##' @seealso \code{\link[graphics]{hist}}
##' @return Returns the histogram.
##' @examples
##' intHist(1:10)
intHist = function(x, range=c(round(min(x, na.rm=TRUE))-0.5, round(max(x, na.rm=TRUE))+0.5), step=1, ...){
  x = shrinkToRange(x, range)
  return(hist(x, breaks=seq(range[1], range[2]+step-1, by = step), ...))
}

##' @title Plots a heatmap
##' @description Plots a heatmap and adds a color key.
##' @param x the data to plot.
##' @template colors-template
##' @param lim two integers used for \code{breaks} argument passed to \code{heatmap.2()}.
##' @param cexCol an integer passed to \code{heatmap.2()}.
##' @param labRow a character vector, possibly modified, then passed to \code{heatmap.2()}.
##' @param margins an integer, possibly modified, then passed to \code{heatmap.2()}.
##' @param dendrogram a character passed to \code{heatmap.2()}.
##' @param Rowv a logical passed to \code{heatmap.2()}.
##' @param Colv a logical passed to \code{heatmap.2()}.
##' @param labCol a character vector passed to \code{heatmap.2()}.
##' @param key a logical passed to \code{heatmap.2()}. Will be set to FALSE if there is only one unique numeric in \code{x}.
##' @template addargs-template
##' @templateVar fun heatmap.2()
##' @template roxygen-template
##' @seealso \code{\link[gplots]{heatmap.2}}
##' @examples
##' ezHeatmap(x=matrix(1:100,10))
ezHeatmap = function(x, lim=c(-4, 4), colors=getBlueRedScale(),
                     dendrogram="none", margins=c(8,6), cexCol=1.1,
                     Rowv=TRUE, Colv=TRUE, labCol=colnames(x), labRow=rownames(x),
                     key=TRUE, ...){
  if (!is.matrix(x)){
    x = as.matrix(x)
  }
  requireNamespace("gplots", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  if (length(unique(as.numeric((x)))) == 1){
    key=FALSE
  }
  heatmap.2(x,
            breaks=seq(from=lim[1], to=lim[2], length.out=257), col=colors, na.color="black",
            Rowv=Rowv, Colv=Colv,
            dendrogram=dendrogram, density.info="none", trace="none", labCol=labCol, labRow=labRow, cexCol=cexCol,
            margins=margins, key=key, ...)
}

## see http://rpubs.com/gaston/dendrograms
##' @title Gets the color cluster labels
##' @description Gets the color cluster labels using \code{dendrapply()}.
##' @param hcd an object of the class dendrogram.
##' @template colors-template
##' @template roxygen-template
##' @seealso \code{\link[stats]{dendrapply}}
##' @return Returns the modified dendrogram.
##' @examples
##' hc = hclust(dist(USArrests), "ave")
##' dend = as.dendrogram(hc)
##' colorClusterLabels(dend, rainbow(6))
colorClusterLabels = function(hcd, colors) {
  dendrapply(hcd, colorNode, cols=colors)
}

##' @describeIn colorClusterLabels The functions used in \code{dendrapply()}.
colorNode = function(n, cols=NULL){
  if (is.leaf(n)) {
    a <- attributes(n)
    attr(n, "nodePar") <- list(lab.col = cols[a$label])
  }
  return(n)
}



########################################################################
## finished?
createDendogramReport <- function(x, annot, genes = row.names(x), multipalette = F, addLegend = F, cex.legend = 1, paletteList = NULL, ...) {
  # Description
  # Wrapper function for plotDendroAndColors -> plot(dendro)
  
  requireNamespace("WGCNA", quietly = T)
  requireNamespace("plyr", quietly = T)
  requireNamespace("pvclust", quietly = T)
  requireNamespace("RColorBrewer", quietly = T)
  requireNamespace("wesanderson", quietly = T)
  
  # Setup different default parameters for plotDendroAndColors arguments if not specified in function call
  # if(!exists("cex.colorLabels")) cex.colorLabels = 1
  # if(!exists("cex.dendroLabels")) cex.dendroLabels = 0.7
  if(!exists("colorHeight")) colorHeight = 0.15
  
  res <- list()
  
  # Subset of genes
  x <- x[genes, ]
  #print(dim(x))
  
  
  
  # Colors for annotation of dendograms
  if(is.null(paletteList)) {
    if(!multipalette) paletteList <- list("grenYll" = c('#4db6ac','#aed581','#dce775','#ffd54f'))
    if(multipalette)  paletteList <- list("Royal1" = wes_palette("Royal1"), 
                                          "Moonrise1" = wes_palette("Moonrise1"),
                                          'Moonrise2' = wes_palette("Moonrise2"),
                                          "Chevalier" = wes_palette("Chevalier"),
                                          "Zissou" = wes_palette("Zissou"),
                                          "Cavalcanti" = wes_palette("Cavalcanti"))
  }
  colList = list()
  for (j in 1:ncol(annot)) {
    gtab <- unique(annot[, j])
    colJ = length(paletteList) - (j %% length(paletteList))
    cols <- colorRampPalette(paletteList[[colJ]])(length(gtab))
    names(cols) = gtab
    colList[[colnames(annot)[j]]] = cols
  }
  
  colAnnot = annot
  for (nm in names(colAnnot)){
    colAnnot[[nm]] = colList[[nm]][annot[[nm]]]
  }
  
  # Dendograms
  d = as.dist(1-cor(x, use="complete.obs"));
  hc = hclust(d, method="ward.D2")
  hcd = as.dendrogram(hc, hang=-0.1)
  
  if(!addLegend) {
    op <- par(no.readonly=TRUE)
    plotDendroAndColors(hc, colAnnot, autoColorHeight = F, ...)
    par(op)
  }
  if(addLegend) {
    opar <- par(no.readonly=TRUE)
    parMar0 <- par()$mar
    layout(matrix(c(1:4), 2, 2), heights = c(1 - colorHeight, colorHeight), widths = c(1 - 0.25, 0.25))
    plotDendroAndColors(hc, colAnnot, 
                        autoColorHeight = F,
                        marAll = c(1, 5, 3, 0),
                        setLayout = FALSE, ...)
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    lNames = gsub ("\\.", " ", names(unlist(colList)))
    legend("center", legend=lNames, fill = unlist(colList), bty = "n", cex = cex.legend)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    par(mar=parMar0)
  }
  
  #   return(res)
  
}


## REFAC, but function is currently unused.
ezArrayImage = function(mat, file=NULL, colorRange=c(-3,3), xScale=1, yScale=1, colors=getBlueRedScale(256)){
  
  mat[mat > colorRange[2]] = colorRange[2]
  mat[mat < colorRange[1]] = colorRange[1]
  
  if (!is.null(file)){
    png(file, height=nrow(mat)*yScale, width=ncol(mat)*xScale)
    on.exit(dev.off())
  }
  par(mar=c(0,0,0,0))
  image(1:ncol(mat), 1:nrow(mat), t(mat), zlim=colorRange, axes=FALSE, frame=FALSE, col=colors,
        xlim=c(0, ncol(mat)), ylim=c(1, nrow(mat)+1))
}


## REFAC, but function is currently unused.
ezCdfPlot = function(x, itemName="gene", scoreName="expression", percentage=FALSE,
                     file=NULL, height=600, width=600, col=getSampleColors(colnames(x), colorNames = colnames(x))){
  if (!is.null(file)){
    switch(sub(".*\\.", "", file),
           pdf=pdf(file, height=height/100, width=width/100),
           png=png(file, height=height, width=width),
           stop("unsupported file: ", file))
    on.exit(dev.off())
  }
  xlim = c(1, max(x))
  xlab = paste(itemName, scoreName)
  if (percentage){
    ylab = paste0("percentage of ", itemName, "s")
    ylim = c(1, 100)
  } else {
    ylab = paste0("number of ", itemName, "s")  
    ylim = c(1, nrow(x))
  }
  
  plot(1, -100, type="l", log="x", xlim=xlim, ylim=ylim,
       xlab=xlab, ylab=ylab, main="Cumulative Distribution")
  for (sm in colnames(x)){
    cts = x[ ,sm]
    cts = cts[!is.na(cts) & cts > 0]
    if (percentage){
      yValue = (1:length(cts)) / length(cts) * 100
    } else {
      yValue = 1:length(cts)			
    }
    lines(sort(cts), yValue, col=sampleColors[sm])
  }
  legend("bottomright", colnames(x), col=sampleColors[colnames(x)], cex=1.2, pt.cex=1.5, pch=20, bty="o", pt.bg="white")
}


## REFAC, but function is currently unused.
goGroupBarPlot = function(xSub){
  ggr = xSub[xSub$Level %in% c("level 2", "unknown"), ]
  pngFile = ezValidFilename(paste0(name, ".png"), replace="-")
  png(pngFile, height=800, width=600)
  par(mar=c(4,16,4,4)) 
  bp = barplot(t(as.matrix(ggr[ ,c("Count", "Size")])), horiz=TRUE, beside=T, 
               legend.text =c("Count", "Size"),
               col=c("blue", "gray"),
               xlab="#Genes", las=2, cex.axis=0.8, cex.names=0.8,
               names.arg=ggr$Term,
               main=name)
  dev.off()
}


## still used?
normalized.distr <- function(data, bins, dmax=1, dmin=-1, ylim= NULL, add=FALSE, col="black", 
                             xlab="value",ylab="freq", pch=1) {
  if (is.null(dmax)) dmax <- max(data)
  if (is.null(dmin)) dmin <- min(data)
  
  breaks <- 0:bins
  breaks <- breaks*(dmax-dmin)/bins
  breaks <- breaks+dmin
  
  x <- vector(length=bins)
  y <- vector(length=bins)
  x[1:bins] <- 0
  y[1:bins] <- 0
  
  for (i in 1:bins)
  {
    if (i==1) idx <- which (data>=breaks[i] & data<=breaks[i+1])
    if (i>1) idx <- which (data>breaks[i] & data<=breaks[i+1])
    x[i] <- (breaks[i]+breaks[i+1])/2
    y[i] <- length(idx)
  }
  
  y <- y/sum(y)
  
  if (!add) 
  {
    plot(x,y,type="l",xlim=c(dmin,dmax), ylim=ylim, col=col, xlab=xlab, ylab=ylab)
    points(x,y,col=col, cex=0.5, pch=pch)
  }
  else 
  {
    lines(x,y, col=col)  
    points(x,y,col=col, cex=0.5, pch=pch)
  }
}


## still used?
ezProfilePlot <- function(x, err=NULL, colors=rainbow(nrow(x)), xaxs="i", yaxs="i", xlim=c(0, ncol(x)+10), ylim=NULL, log="", type="l",
                          names=rownames(x), legendPos="topright", lty=rep(1, nrow(x)), main="", plotXLabels=TRUE, xlab="", ylab="", ...){
  
  if (is.null(ylim)){
    ylim = range(x, na.rm=TRUE)
    if (!is.null(err)){
      ylim = range(x+err, x-err, na.rm=TRUE)
    }
  }
  
  x <- as.matrix(x)
  for(i in 1:nrow(x)){
    if (i==1){
      plot(x[i,], col=colors[i],
           xlim=xlim, ylim=ylim, xaxs=xaxs, yaxs=yaxs,
           type=type, axes=FALSE, frame=TRUE,
           xlab=xlab, ylab=ylab, log=log, lty=lty[i], main=main, ...)
    } else {
      lines(x[i,], col=colors[i], lty=lty[i], type=type, ...)
    }
    if (!is.null(err)){
      points(1:ncol(x), x[i, ], pch=16, col=colors[i])
      arrows(1:ncol(x), x[i, ] - err[i, ], 1:ncol(x), x[i, ] + err[i, ], col=colors[i],
             angle=90, code=3, length=0.1, ...)
    }
  }
  if (!is.null(legendPos)){
    legend(legendPos, names, col=colors, bty="n", cex=0.5, pt.bg="white", lty=lty, lwd=3, seg.len=4 )
  }
  if (plotXLabels){
    axis(1, at=1:ncol(x), labels=colnames(x), las=2, cex.axis=0.8)
  }
  axis(2)
  return()
}


## still used?
getBinColors = function(binNames, colorSet=c("darkorange", "gray70", "gray50", "gray70", "cyan")){
  colors = colorRampPalette(colorSet)(length(binNames))
  names(colors) = binNames
  return(colors)
}
