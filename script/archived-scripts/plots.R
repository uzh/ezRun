
ezVolcano = function(log2Ratio, pValue, yType="p-value",
                     xlim=NULL, ylim=NULL, isPresent=NULL,
                     types=NULL, pch=16, colors=rainbow(ncol(types)), legendPos="bottomright",
                     cex.main=1.0, cex=1, ...){
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
  return(list(x=log2Ratio, y=yValues))
}
